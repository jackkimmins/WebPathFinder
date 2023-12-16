#include <emscripten.h>
#include <emscripten/bind.h>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <sstream>

class City {
public:
    std::string name;
    int index;

    City(std::string name, int index) : name(std::move(name)), index(index) {}
    City() = default;
};

struct Route {
    std::vector<City> cities;
    double distance = 0.0;

    explicit Route(std::vector<City> cities) : cities(std::move(cities)) {}

    // Calculate the total distance of the route based on a given distance matrix
    void CalculateDistance(const std::vector<std::vector<double>>& distanceMatrix) {
        distance = 0.0;
        for (size_t i = 0, n = cities.size(); i < n - 1; ++i) {
            distance += distanceMatrix[cities[i].index][cities[i + 1].index];
        }
        // Add distance from last city to the first city to complete the loop
        distance += distanceMatrix[cities.back().index][cities.front().index];
    }

    // String representation of the route
    std::string ToString() const {
        std::ostringstream routeStr;
        for (size_t i = 0, n = cities.size(); i < n; ++i) {
            routeStr << cities[i].name << (i < n - 1 ? " -> " : "");
        }
        return routeStr.str();
    }

    bool operator<(const Route& other) const {
        return distance < other.distance;
    }
};

class Population {
public:
    std::vector<Route> routes;

    Population(size_t size) : routes() {
        routes.reserve(size);
    }
};

class TSPSolver {
private:
    bool isComplete = false;                            // Flag to check if the solution is complete
    std::vector<std::vector<double>> distanceMatrix;    // Matrix to store distances between cities
    std::vector<City> cities;                           // Vector of cities
    std::mt19937 gen;                                   // Mersenne Twister random number generator
    Population population;                              // Population of routes
    std::uniform_real_distribution<double> dist;        // Distribution for real numbers
    std::uniform_int_distribution<> intDist;            // Distribution for integers

    size_t currentGeneration = 0;                       // Current generation count in the genetic algorithm

    // Constants for the genetic algorithm
    const size_t NUM_GENERATIONS = 1000;
    const double MUTATION_RATE = 0.05;

    // Generate a random real number within a range
    double GetRandomNumber(double min, double max) {
        return dist(gen, decltype(dist)::param_type{min, max});
    }

    // Generate a random index within a range
    int GetRandomIndex(size_t size) {
        return intDist(gen, decltype(intDist)::param_type{0, static_cast<int>(size - 1)});
    }

    // Generate a route using the Nearest Neighbour algorithm
    Route NearestNeighbourRoute(int startCityIndex) {
        std::vector<City> nnRoute;
        nnRoute.reserve(cities.size());

        std::unordered_set<int> visitedCities;
        visitedCities.insert(startCityIndex);
        nnRoute.push_back(cities[startCityIndex]);

        int currentCityIndex = startCityIndex;
        // Loop to construct the route
        while (nnRoute.size() < cities.size()) {
            double minDistance = std::numeric_limits<double>::max();
            int nearestCityIndex = -1;

            // Find the nearest unvisited city
            for (int i = 0; i < cities.size(); ++i) {
                if (visitedCities.find(i) == visitedCities.end() && distanceMatrix[currentCityIndex][i] < minDistance) {
                    minDistance = distanceMatrix[currentCityIndex][i];
                    nearestCityIndex = i;
                }
            }

            nnRoute.push_back(cities[nearestCityIndex]);
            visitedCities.insert(nearestCityIndex);
            currentCityIndex = nearestCityIndex;
        }

        return Route(nnRoute);
    }

    // Calculate the total inverse distance of all routes in a population
    double GetTotalInverseDistance(const Population& population) {
        double totalInverseDistance = 0.0;
        for (const auto& route : population.routes) {
            if (route.distance > 0) {
                totalInverseDistance += 1.0 / route.distance;
            }
        }
        return totalInverseDistance;
    }

    // Select a route from the population using the Roulette Wheel Selection method
    Route RouletteWheelSelection(const Population& population) {
        double totalInverseDistance = GetTotalInverseDistance(population);
        double slice = GetRandomNumber(0, totalInverseDistance);
        double currentSum = 0.0;

        for (const auto& route : population.routes) {
            if (route.distance > 0) {
                currentSum += 1.0 / route.distance;
                if (currentSum >= slice) {
                    return route;
                }
            }
        }
        // In case of rounding errors, return a random route
        return population.routes[GetRandomIndex(population.routes.size())];
    }

    // Perform genetic crossover between two parent routes
    Route Crossover(const Route& parent1, const Route& parent2) {
        int start = GetRandomIndex(parent1.cities.size());
        int end = GetRandomNumber(start, parent1.cities.size());

        std::vector<City> childCities(parent1.cities.size());
        std::unordered_set<int> included;

        // Copy a segment from parent1 to child
        for (int i = start; i < end; ++i) {
            childCities[i] = parent1.cities[i];
            included.insert(parent1.cities[i].index);
        }

        // Fill the rest of child with cities from parent2
        size_t curIndex = end % parent1.cities.size();
        for (const auto& city : parent2.cities) {
            if (included.find(city.index) == included.end()) {
                childCities[curIndex] = city;
                curIndex = (curIndex + 1) % parent1.cities.size();
            }
        }

        return Route(childCities);
    }

    // Calculate the diversity of a population
    double CalculateDiversity(const Population& population) {
        double meanDistance = 0.0;
        // Calculate mean distance of all routes
        for (const auto& route : population.routes) {
            meanDistance += route.distance;
        }
        meanDistance /= population.routes.size();

        double variance = 0.0;
        // Calculate variance of distances
        for (const auto& route : population.routes) {
            double diff = route.distance - meanDistance;
            variance += diff * diff;
        }
        variance /= population.routes.size();
        return variance;
    }

    // Mutate a route
    void Mutate(Route& route) {
        double diversity = CalculateDiversity(population);
        double adaptiveMutationRate = MUTATION_RATE;

        // Increase mutation rate if diversity is low
        if (diversity < 0.6) {
            adaptiveMutationRate *= 2;
        }

        // Mutate each city in the route
        for (size_t i = 0; i < route.cities.size(); ++i) {
            if (GetRandomNumber(0, 1) < adaptiveMutationRate) {
                size_t swapIndex = GetRandomIndex(route.cities.size());
                std::swap(route.cities[i], route.cities[swapIndex]);
            }
        }
        route.CalculateDistance(distanceMatrix);
    }

    // Initialise the population with nearest neighbour routes
    void InitialisePopulation(size_t populationSize) {
        population.routes.clear();

        for (size_t i = 0; i < populationSize; ++i) {
            int startCityIndex = GetRandomIndex(cities.size());
            Route nnRoute = NearestNeighbourRoute(startCityIndex);
            nnRoute.CalculateDistance(distanceMatrix);
            population.routes.push_back(nnRoute);
        }
    }

    // Apply the 3-opt optimisation algorithm on a route
    void ThreeOpt(Route& route) {
        bool improvement = true;
        while (improvement) {
            improvement = false;
            for (size_t i = 0; i < route.cities.size() - 2; ++i) {
                for (size_t j = i + 1; j < route.cities.size() - 1; ++j) {
                    for (size_t k = j + 1; k < route.cities.size(); ++k) {
                        double currentDistance = distanceMatrix[route.cities[i].index][route.cities[i + 1].index] +
                                                distanceMatrix[route.cities[j].index][route.cities[j + 1].index] +
                                                distanceMatrix[route.cities[k].index][route.cities[(k + 1) % route.cities.size()].index];

                        double newDistance = distanceMatrix[route.cities[i].index][route.cities[j].index] +
                                            distanceMatrix[route.cities[i + 1].index][route.cities[k].index] +
                                            distanceMatrix[route.cities[j + 1].index][route.cities[(k + 1) % route.cities.size()].index];

                        if (newDistance < currentDistance) {
                            std::reverse(route.cities.begin() + i + 1, route.cities.begin() + j + 1);
                            std::reverse(route.cities.begin() + j + 1, route.cities.begin() + k + 1);
                            route.CalculateDistance(distanceMatrix);
                            improvement = true;
                        }
                    }
                }
            }
        }
    }

    void InitialiseFromCSV(const std::string& csvData) {
        std::istringstream ss(csvData);
        std::string line;

        // Read the header line to get city names
        if (!std::getline(ss, line)) {
            std::cerr << "Error reading CSV data." << std::endl;
            return;
        }

        std::vector<City> cityNames;
        std::string cityName;
        size_t index = 0;

        // Read city names from the header, skipping the first empty cell
        std::istringstream headerStream(line);
        std::getline(headerStream, cityName, ','); // Skip first cell
        while (std::getline(headerStream, cityName, ',')) {
            cityNames.emplace_back(cityName, index++);
        }

        // Read the distances
        std::vector<std::vector<double>> localDistanceMatrix;
        double distance;
        char comma;

        while (std::getline(ss, line)) {
            std::istringstream lineStream(line);
            std::vector<double> distances;

            std::getline(lineStream, cityName, ','); // Read and ignore city name

            while (lineStream >> distance) {
                distances.push_back(distance);
                lineStream >> comma; // Read and ignore comma
            }

            localDistanceMatrix.push_back(std::move(distances));
        }

        this->cities = std::move(cityNames);
        this->distanceMatrix = std::move(localDistanceMatrix);
    }

public:
    TSPSolver(const std::vector<std::vector<double>>& matrix, const std::vector<City>& inputCities, size_t populationSize) 
    : distanceMatrix(matrix), cities(inputCities), population(populationSize),
      dist(0.0, 1.0), intDist(0, static_cast<int>(inputCities.size() - 1)) {
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        gen.seed(seed);
        InitialisePopulation(populationSize);
    }

    TSPSolver(const std::string& csvData, size_t populationSize) 
    : population(populationSize), dist(0.0, 1.0) {
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        gen.seed(seed);
        InitialiseFromCSV(csvData);
        InitialisePopulation(populationSize);
    }

    bool Run() {
        size_t numElites = 0.05 * population.routes.size();

        for (currentGeneration = 0; currentGeneration < NUM_GENERATIONS; ++currentGeneration) {
            Population newPopulation(population.routes.size());

            // Sort the current population by route distance
            std::sort(population.routes.begin(), population.routes.end());

            // Carry over elite routes
            for (size_t i = 0; i < numElites; ++i) {
                newPopulation.routes.push_back(population.routes[i]);
            }

            // Fill the rest of the new population
            while (newPopulation.routes.size() < population.routes.size()) {
                Route parent1 = RouletteWheelSelection(population);
                Route parent2 = RouletteWheelSelection(population);
                Route child = Crossover(parent1, parent2);
                Mutate(child);
                ThreeOpt(child); // Apply 3-opt optimisation
                newPopulation.routes.push_back(std::move(child));
            }

            population = std::move(newPopulation);
        }

        isComplete = true;
        return true;
    }

    std::string GetBestRoute() {
        if (!isComplete) return "Not complete";

        const Route& bestRoute = *std::min_element(population.routes.begin(), population.routes.end());

        std::stringstream ss;
        for (const auto& city : bestRoute.cities) {
            std::string name = city.name;
            ss << name << std::endl;
        }
        return ss.str();
    }

    int GetBestRouteLength() {
        if (!isComplete) return -1;

        const Route& bestRoute = *std::min_element(population.routes.begin(), population.routes.end());
        return bestRoute.distance;
    }
};

EMSCRIPTEN_BINDINGS(tsp_class) {
    emscripten::class_<TSPSolver>("TSPSolver")
    .constructor<const std::string&, size_t>()
    .function("run", &TSPSolver::Run)
    .function("GetBestRoute", &TSPSolver::GetBestRoute)
    .function("GetBestRouteLength", &TSPSolver::GetBestRouteLength)
    ;
};

int main()
{
    return 0;
}