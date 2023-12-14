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
    int index = -1;

    City(const std::string& name, int index) : name(name), index(index) {}
    City() = default;
};

class Route {
public:
    std::vector<City> cities;
    double distance = 0.0;

    Route(const std::vector<City>& cities) : cities(cities) {}

    void calculateDistance(const std::vector<std::vector<double>>& distanceMatrix) {
        distance = 0;
        for (size_t i = 0; i < cities.size() - 1; ++i) {
            distance += distanceMatrix[cities[i].index][cities[i + 1].index];
        }
        distance += distanceMatrix[cities.back().index][cities.front().index];
    }

    std::string toString() const {
        std::ostringstream routeStr;
        for (size_t i = 0; i < cities.size(); ++i) {
            routeStr << cities[i].name << (i < cities.size() - 1 ? " -> " : "");
        }
        return routeStr.str();
    }

    bool operator<(const Route& other) const {
        return this->distance < other.distance;
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
    bool isComplete = false;
    std::vector<std::vector<double>> distanceMatrix;
    std::vector<City> cities;
    std::mt19937 gen;
    Population population;
    std::uniform_real_distribution<double> dist;
    std::uniform_int_distribution<> intDist;

    const size_t NUM_GENERATIONS = 1000;
    const double MUTATION_RATE = 0.05;
    const size_t TOURNAMENT_SELECTION_SIZE = 5;

    double getRandomNumber(double min, double max) {
        return dist(gen, decltype(dist)::param_type{min, max});
    }

    int getRandomIndex(size_t size) {
        return intDist(gen, decltype(intDist)::param_type{0, static_cast<int>(size - 1)});
    }

    Route tournamentSelection(const Population& population) {
        Population tournament(TOURNAMENT_SELECTION_SIZE);
        for (size_t i = 0; i < TOURNAMENT_SELECTION_SIZE; ++i) {
            tournament.routes.push_back(population.routes[getRandomIndex(population.routes.size())]);
        }
        return *std::min_element(tournament.routes.begin(), tournament.routes.end());
    }

    Route crossover(const Route& parent1, const Route& parent2) {
        int start = getRandomIndex(parent1.cities.size());
        int end = getRandomNumber(start, parent1.cities.size());

        std::vector<City> childCities(parent1.cities.size());
        std::unordered_set<int> included;

        for (int i = start; i < end; ++i) {
            childCities[i] = parent1.cities[i];
            included.insert(parent1.cities[i].index);
        }

        size_t curIndex = end % parent1.cities.size();
        for (const auto& city : parent2.cities) {
            if (included.find(city.index) == included.end()) {
                childCities[curIndex] = city;
                curIndex = (curIndex + 1) % parent1.cities.size();
            }
        }

        return Route(childCities);
    }

    void mutate(Route& route) {
        for (size_t i = 0; i < route.cities.size(); ++i) {
            if (getRandomNumber(0, 1) < MUTATION_RATE) {
                size_t swapIndex = getRandomIndex(route.cities.size());
                std::swap(route.cities[i], route.cities[swapIndex]);
            }
        }
        route.calculateDistance(distanceMatrix);
    }

    // Population Initialization
    void initializePopulation(size_t populationSize) {
        for (size_t i = 0; i < populationSize; ++i) {
            std::shuffle(cities.begin(), cities.end(), gen);
            Route newRoute(cities);
            newRoute.calculateDistance(distanceMatrix);
            population.routes.push_back(newRoute);
        }
    }

void initializeFromCSV(const std::string& csvData) {
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
        initializePopulation(populationSize);
    }

    TSPSolver(const std::string& csvData, size_t populationSize) 
    : population(populationSize), dist(0.0, 1.0) {
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        gen.seed(seed);
        initializeFromCSV(csvData);
        initializePopulation(populationSize);
    }

    bool run() {
        for (size_t gen = 0; gen < NUM_GENERATIONS; ++gen) {
            Population newPopulation(population.routes.size());

            for (size_t i = 0; i < population.routes.size(); ++i) {
                Route parent1 = tournamentSelection(population);
                Route parent2 = tournamentSelection(population);
                Route child = crossover(parent1, parent2);
                mutate(child);
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
            // ss << city.name << std::endl;

            //Remove the quotes
            std::string name = city.name;
            // name.erase(std::remove(name.begin(), name.end(), ''), name.end());
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
    .function("run", &TSPSolver::run)
    .function("GetBestRoute", &TSPSolver::GetBestRoute)
    .function("GetBestRouteLength", &TSPSolver::GetBestRouteLength)
    ;
};

int main()
{
    return 0;
}