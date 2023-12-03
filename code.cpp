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

    City(const std::string& name, int index) : name(name), index(index) {}
};

class Route {
public:
    std::vector<City> cities;
    double distance;

    Route(const std::vector<City>& cities) : cities(cities), distance(0) {}

    void calculateDistance(const std::vector<std::vector<double>>& distanceMatrix) {
        distance = 0;
        for (size_t i = 0; i < cities.size() - 1; ++i) {
            distance += distanceMatrix[cities[i].index][cities[i + 1].index];
        }
        distance += distanceMatrix[cities.back().index][cities.front().index];
    }

    std::string toString() const {
        std::string routeStr;
        for (const auto& city : cities) {
            routeStr += city.name + " -> ";
        }
        routeStr += cities.front().name; // Loop back to the start
        return routeStr;
    }
};

class Population {
public:
    std::vector<Route> routes;

    Population(size_t size) {
        routes.reserve(size);
    }

    void addRoute(const Route& route) {
        routes.push_back(route);
    }
};

class TSPSolver {
private:
    std::vector<std::vector<double>> distanceMatrix;
    std::vector<City> cities; // Added to store cities
    Population population;
    std::mt19937 gen;
    int maxIterations = 1000; // Example value for maximum iterations
    int tournamentSize = 5; // Example value for tournament size

    void initializePopulation() {
        for (size_t i = 0; i < population.routes.capacity(); ++i) {
            std::vector<City> shuffledCities = cities;
            std::shuffle(shuffledCities.begin(), shuffledCities.end(), gen);
            Route newRoute(shuffledCities);
            newRoute.calculateDistance(distanceMatrix);
            population.addRoute(newRoute);
        }
    }

    Route& tournamentSelection() {
        std::uniform_int_distribution<size_t> dist(0, population.routes.size() - 1);
        Route& best = population.routes[dist(gen)];
        for (int i = 0; i < tournamentSize; ++i) {
            Route& rival = population.routes[dist(gen)];
            if (rival.distance < best.distance) {
                best = rival;
            }
        }
        return best;
    }

    void crossover(const Route& parent1, const Route& parent2, Route& child) {
        std::uniform_int_distribution<size_t> dist(0, parent1.cities.size() - 1);
        size_t start = dist(gen);
        size_t end = dist(gen);

        if (start > end) std::swap(start, end);

        std::unordered_set<int> included;
        for (size_t i = start; i <= end; ++i) {
            child.cities[i] = parent1.cities[i];
            included.insert(parent1.cities[i].index);
        }

        size_t currentIndex = (end + 1) % parent1.cities.size();
        for (const auto& city : parent2.cities) {
            if (included.find(city.index) == included.end()) {
                child.cities[currentIndex] = city;
                currentIndex = (currentIndex + 1) % parent1.cities.size();
            }
        }
        child.calculateDistance(distanceMatrix);
    }

    void mutate(Route& route) {
        std::uniform_int_distribution<size_t> dist(0, route.cities.size() - 1);
        size_t i = dist(gen);
        size_t j = dist(gen);
        std::swap(route.cities[i], route.cities[j]);
        route.calculateDistance(distanceMatrix);
    }

    void twoOpt(Route& route) {
        bool improvement = true;
        while (improvement) {
            improvement = false;
            for (size_t i = 0; i < route.cities.size() - 1; ++i) {
                for (size_t k = i + 1; k < route.cities.size(); ++k) {
                    double delta = - distanceMatrix[route.cities[i].index][route.cities[(i + 1) % route.cities.size()].index]
                                   - distanceMatrix[route.cities[k].index][route.cities[(k + 1) % route.cities.size()].index]
                                   + distanceMatrix[route.cities[i].index][route.cities[k].index]
                                   + distanceMatrix[route.cities[(i + 1) % route.cities.size()].index][route.cities[(k + 1) % route.cities.size()].index];

                    if (delta < 0) {
                        std::reverse(route.cities.begin() + i + 1, route.cities.begin() + k + 1);
                        route.calculateDistance(distanceMatrix);
                        improvement = true;
                    }
                }
            }
        }
    }

    void initializeFromCSV(const std::string& csvData) {
        std::stringstream ss(csvData);
        std::string line;

        // First, read the header line to get city names
        if (!std::getline(ss, line)) {
            std::cerr << "Error reading CSV data." << std::endl;
            return;
        }

        std::stringstream headerStream(line);
        std::string cityName;

        // Skip the first empty cell in the header
        std::getline(headerStream, cityName, ',');

        std::vector<City> cityNames;
        int index = 0;
        while (std::getline(headerStream, cityName, ',')) {
            cityNames.emplace_back(cityName, index++);
        }

        // Now read the distances
        std::vector<std::vector<double>> distanceMatrix;
        while (std::getline(ss, line)) {
            std::stringstream lineStream(line);
            std::vector<double> distances;

            // Read the city name (but don't use it)
            std::getline(lineStream, cityName, ',');

            double distance;
            while (lineStream >> distance) {
                distances.push_back(distance);

                // If the next character is a comma, ignore it
                if (lineStream.peek() == ',') {
                    lineStream.ignore();
                }
            }

            distanceMatrix.push_back(distances);
        }

        // Initialize the solver with the read data
        this->cities = cityNames;
        this->distanceMatrix = distanceMatrix;
        initializePopulation(); // Assuming this reinitializes the population
    }

public:
    TSPSolver(const std::vector<std::vector<double>>& matrix, const std::vector<City>& inputCities, size_t populationSize) : distanceMatrix(matrix), cities(inputCities), population(populationSize) {
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        gen = std::mt19937(seed);
        initializePopulation();
    }

    TSPSolver(const std::string& csvData, size_t populationSize) : population(populationSize) {
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        gen = std::mt19937(seed);
        initializeFromCSV(csvData);
    }

    std::string run() {
        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            // Perform crossover and mutation on selected individuals
            for (size_t i = 0; i < population.routes.size(); ++i) {
                Route& parent1 = tournamentSelection();
                Route& parent2 = tournamentSelection();
                Route child(cities);
                crossover(parent1, parent2, child);
                mutate(child);
                twoOpt(child); // Apply 2-opt local search
                population.routes[i] = child;
            }
            // Sort population by fitness (route distance)
            std::sort(population.routes.begin(), population.routes.end(), [](const Route& a, const Route& b) {
                return a.distance < b.distance;
            });
        }
        // Return the best route found
        // return population.routes[0];

        //Return as csv
        std::stringstream ss;
        for (const auto& city : population.routes[0].cities) {
            // ss << city.name << std::endl;

            //Remove the quotes
            std::string name = city.name;
            name.erase(std::remove(name.begin(), name.end(), '\"'), name.end());
            ss << name << std::endl;
        }
        return ss.str();
    }
};

EMSCRIPTEN_BINDINGS(tsp_class) {
    emscripten::class_<TSPSolver>("TSPSolver")
    .constructor<const std::string&, size_t>()
    .function("run", &TSPSolver::run)
    ;
};

int main()
{
    // Create distance matrix
    // std::vector<std::vector<double>> distanceMatrix = {
    //     { 0, 1, 2, 3, 4 },
    //     { 1, 0, 1, 2, 3 },
    //     { 2, 1, 0, 1, 2 },
    //     { 3, 2, 1, 0, 1 },
    //     { 4, 3, 2, 1, 0 }
    // };

    // // Create cities
    // std::vector<City> cities = {
    //     City("A", 0),
    //     City("B", 1),
    //     City("C", 2),
    //     City("D", 3),
    //     City("E", 4)
    // };

    // std::string csvData = " ,Taunton,St Ives,Bristol,London,Middlesbrough,Weston-super-Mare\nTaunton,0,229,76,245,491,46\nSt Ives,228,0,301,456,716,271\nBristol,75,301,0,191,425,37\nLondon,244,456,190,0,404,227\nMiddlesbrough,490,716,424,403,0,452\nWeston-super-Mare,46,272,37,227,452,0";

    // TSPSolver solver(csvData, 100);
    // Route finalRoute = solver.run();
    // std::cout << "Best Route: " << finalRoute.toString() << std::endl;



    return 0;
}