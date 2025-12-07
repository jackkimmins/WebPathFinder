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
#include <limits>
#include <numeric>

class City {
public:
    std::string name;
    int index;

    City(std::string name, int index) : name(std::move(name)), index(index) {}
    City() = default;
};

class TSPSolver {
private:
    bool isComplete = false;
    std::vector<std::vector<double>> distanceMatrix;
    std::vector<City> cities;
    std::vector<int> bestPath;
    double bestDistance = std::numeric_limits<double>::max();
    std::mt19937 gen;

    // Threshold for using exact algorithms (can handle up to ~12 cities with brute force quickly)
    static const size_t EXACT_THRESHOLD = 12;
    // Threshold for branch and bound (can handle up to ~18 cities reasonably)
    static const size_t BRANCH_BOUND_THRESHOLD = 18;

    // Calculate path distance (open path - no return to start)
    double CalculatePathDistance(const std::vector<int>& path) const {
        double dist = 0.0;
        for (size_t i = 0; i < path.size() - 1; ++i) {
            dist += distanceMatrix[path[i]][path[i + 1]];
        }
        return dist;
    }

    // Calculate tour distance (closed loop - returns to start)
    double CalculateTourDistance(const std::vector<int>& path) const {
        double dist = CalculatePathDistance(path);
        if (!path.empty()) {
            dist += distanceMatrix[path.back()][path.front()];
        }
        return dist;
    }

    // ==================== EXACT ALGORITHMS ====================

    // Brute force - guaranteed optimal for small instances
    void SolveBruteForce() {
        std::vector<int> path(cities.size());
        std::iota(path.begin(), path.end(), 0);

        bestDistance = std::numeric_limits<double>::max();
        
        do {
            double dist = CalculateTourDistance(path);
            if (dist < bestDistance) {
                bestDistance = dist;
                bestPath = path;
            }
        } while (std::next_permutation(path.begin(), path.end()));
    }

    // Branch and Bound with pruning - optimal for medium instances
    void BranchAndBound() {
        size_t n = cities.size();
        std::vector<int> currentPath;
        std::vector<bool> visited(n, false);
        
        // Start with a good upper bound using nearest neighbor heuristic
        std::vector<int> nnPath = NearestNeighborPath(0);
        bestDistance = CalculateTourDistance(nnPath);
        bestPath = nnPath;
        
        // Try starting from each city and apply 2-opt to improve
        for (size_t start = 0; start < n; ++start) {
            std::vector<int> path = NearestNeighborPath(start);
            TwoOptImprove(path);
            double dist = CalculateTourDistance(path);
            if (dist < bestDistance) {
                bestDistance = dist;
                bestPath = path;
            }
        }
        
        // Now do branch and bound with this upper bound
        currentPath.push_back(0);
        visited[0] = true;
        BranchAndBoundRecursive(currentPath, visited, 0.0);
    }

    void BranchAndBoundRecursive(std::vector<int>& currentPath, std::vector<bool>& visited, double currentDist) {
        size_t n = cities.size();
        
        if (currentPath.size() == n) {
            // Complete tour - add return to start
            double totalDist = currentDist + distanceMatrix[currentPath.back()][currentPath.front()];
            if (totalDist < bestDistance) {
                bestDistance = totalDist;
                bestPath = currentPath;
            }
            return;
        }

        // Calculate lower bound for remaining path
        double lowerBound = currentDist + CalculateLowerBound(visited, currentPath.back());
        if (lowerBound >= bestDistance) {
            return; // Prune this branch
        }

        // Try adding each unvisited city
        for (size_t i = 0; i < n; ++i) {
            if (!visited[i]) {
                double newDist = currentDist + distanceMatrix[currentPath.back()][i];
                
                // Early pruning
                if (newDist >= bestDistance) continue;
                
                visited[i] = true;
                currentPath.push_back(i);
                
                BranchAndBoundRecursive(currentPath, visited, newDist);
                
                currentPath.pop_back();
                visited[i] = false;
            }
        }
    }

    // Calculate lower bound using minimum edge heuristic
    double CalculateLowerBound(const std::vector<bool>& visited, int lastCity) const {
        double bound = 0.0;
        size_t n = cities.size();
        
        // Minimum edge from last city to any unvisited city
        double minFromLast = std::numeric_limits<double>::max();
        for (size_t i = 0; i < n; ++i) {
            if (!visited[i]) {
                minFromLast = std::min(minFromLast, distanceMatrix[lastCity][i]);
            }
        }
        if (minFromLast < std::numeric_limits<double>::max()) {
            bound += minFromLast;
        }

        // For each unvisited city, add minimum outgoing edge
        for (size_t i = 0; i < n; ++i) {
            if (!visited[i]) {
                double minEdge = std::numeric_limits<double>::max();
                for (size_t j = 0; j < n; ++j) {
                    if (i != j && (!visited[j] || j == 0)) { // Can go to unvisited or back to start
                        minEdge = std::min(minEdge, distanceMatrix[i][j]);
                    }
                }
                if (minEdge < std::numeric_limits<double>::max()) {
                    bound += minEdge;
                }
            }
        }

        return bound;
    }

    // ==================== HEURISTIC ALGORITHMS ====================

    // Nearest Neighbor heuristic
    std::vector<int> NearestNeighborPath(int startCity) const {
        size_t n = cities.size();
        std::vector<int> path;
        path.reserve(n);
        std::vector<bool> visited(n, false);
        
        path.push_back(startCity);
        visited[startCity] = true;
        
        while (path.size() < n) {
            int current = path.back();
            double minDist = std::numeric_limits<double>::max();
            int nearest = -1;
            
            for (size_t i = 0; i < n; ++i) {
                if (!visited[i] && distanceMatrix[current][i] < minDist) {
                    minDist = distanceMatrix[current][i];
                    nearest = i;
                }
            }
            
            path.push_back(nearest);
            visited[nearest] = true;
        }
        
        return path;
    }

    // 2-opt improvement - very effective local search
    bool TwoOptImprove(std::vector<int>& path) const {
        size_t n = path.size();
        bool improved = true;
        bool anyImprovement = false;
        
        while (improved) {
            improved = false;
            for (size_t i = 0; i < n - 1; ++i) {
                for (size_t j = i + 2; j < n; ++j) {
                    // Calculate improvement from reversing segment [i+1, j]
                    size_t nextJ = (j + 1) % n;
                    
                    double oldDist = distanceMatrix[path[i]][path[i + 1]] +
                                     distanceMatrix[path[j]][path[nextJ]];
                    double newDist = distanceMatrix[path[i]][path[j]] +
                                     distanceMatrix[path[i + 1]][path[nextJ]];
                    
                    if (newDist < oldDist - 1e-10) {
                        std::reverse(path.begin() + i + 1, path.begin() + j + 1);
                        improved = true;
                        anyImprovement = true;
                    }
                }
            }
        }
        return anyImprovement;
    }

    // Or-opt: move a sequence of 1, 2, or 3 consecutive cities
    bool OrOptImprove(std::vector<int>& path) const {
        size_t n = path.size();
        bool improved = true;
        bool anyImprovement = false;
        
        while (improved) {
            improved = false;
            // Try moving segments of length 1, 2, and 3
            for (int segLen = 1; segLen <= 3 && segLen <= (int)n - 2; ++segLen) {
                for (size_t i = 0; i < n - segLen; ++i) {
                    for (size_t j = 0; j < n; ++j) {
                        if (j >= i && j <= i + segLen) continue;
                        
                        // Calculate cost of removing segment [i, i+segLen-1] and inserting after j
                        size_t prevI = (i + n - 1) % n;
                        size_t nextSeg = (i + segLen) % n;
                        size_t nextJ = (j + 1) % n;
                        
                        double oldCost = distanceMatrix[path[prevI]][path[i]] +
                                         distanceMatrix[path[i + segLen - 1]][path[nextSeg]] +
                                         distanceMatrix[path[j]][path[nextJ]];
                        
                        double newCost = distanceMatrix[path[prevI]][path[nextSeg]] +
                                         distanceMatrix[path[j]][path[i]] +
                                         distanceMatrix[path[i + segLen - 1]][path[nextJ]];
                        
                        if (newCost < oldCost - 1e-10) {
                            // Perform the move
                            std::vector<int> segment(path.begin() + i, path.begin() + i + segLen);
                            path.erase(path.begin() + i, path.begin() + i + segLen);
                            
                            size_t insertPos = (j > i) ? j - segLen + 1 : j + 1;
                            path.insert(path.begin() + insertPos, segment.begin(), segment.end());
                            
                            improved = true;
                            anyImprovement = true;
                            break;
                        }
                    }
                    if (improved) break;
                }
                if (improved) break;
            }
        }
        return anyImprovement;
    }

    // Lin-Kernighan style improvement with multiple restarts
    void SolveWithLKHeuristic() {
        size_t n = cities.size();
        bestDistance = std::numeric_limits<double>::max();
        
        // Generate multiple starting solutions and improve each
        for (size_t start = 0; start < n; ++start) {
            std::vector<int> path = NearestNeighborPath(start);
            
            // Apply iterative improvement
            bool improved = true;
            while (improved) {
                improved = TwoOptImprove(path);
                improved = OrOptImprove(path) || improved;
            }
            
            double dist = CalculateTourDistance(path);
            if (dist < bestDistance) {
                bestDistance = dist;
                bestPath = path;
            }
        }
        
        // Also try random starting permutations
        size_t numRandomStarts = std::min((size_t)50, n * n);
        std::vector<int> randomPath(n);
        std::iota(randomPath.begin(), randomPath.end(), 0);
        
        for (size_t r = 0; r < numRandomStarts; ++r) {
            std::shuffle(randomPath.begin(), randomPath.end(), gen);
            std::vector<int> path = randomPath;
            
            bool improved = true;
            while (improved) {
                improved = TwoOptImprove(path);
                improved = OrOptImprove(path) || improved;
            }
            
            double dist = CalculateTourDistance(path);
            if (dist < bestDistance) {
                bestDistance = dist;
                bestPath = path;
            }
        }
    }

    void InitialiseFromCSV(const std::string& csvData) {
        std::istringstream ss(csvData);
        std::string line;

        if (!std::getline(ss, line)) {
            std::cerr << "Error reading CSV data." << std::endl;
            return;
        }

        std::vector<City> cityNames;
        std::string cityName;
        size_t index = 0;

        std::istringstream headerStream(line);
        std::getline(headerStream, cityName, ',');
        while (std::getline(headerStream, cityName, ',')) {
            cityNames.emplace_back(cityName, index++);
        }

        std::vector<std::vector<double>> localDistanceMatrix;
        double distance;
        char comma;

        while (std::getline(ss, line)) {
            std::istringstream lineStream(line);
            std::vector<double> distances;

            std::getline(lineStream, cityName, ',');

            while (lineStream >> distance) {
                distances.push_back(distance);
                lineStream >> comma;
            }

            localDistanceMatrix.push_back(std::move(distances));
        }

        this->cities = std::move(cityNames);
        this->distanceMatrix = std::move(localDistanceMatrix);
    }

public:
    TSPSolver(const std::vector<std::vector<double>>& matrix, const std::vector<City>& inputCities, size_t /* populationSize - ignored */) 
    : distanceMatrix(matrix), cities(inputCities) {
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        gen.seed(seed);
    }

    TSPSolver(const std::string& csvData, size_t /* populationSize - ignored */) {
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        gen.seed(seed);
        InitialiseFromCSV(csvData);
    }

    bool Run() {
        size_t n = cities.size();
        
        if (n <= 1) {
            bestPath.clear();
            if (n == 1) bestPath.push_back(0);
            bestDistance = 0;
            isComplete = true;
            return true;
        }
        
        if (n <= EXACT_THRESHOLD) {
            // Use brute force for small instances - guaranteed optimal
            SolveBruteForce();
        } else if (n <= BRANCH_BOUND_THRESHOLD) {
            // Use branch and bound for medium instances - guaranteed optimal but slower
            BranchAndBound();
        } else {
            // Use LK-style heuristic for larger instances - very good but not guaranteed optimal
            SolveWithLKHeuristic();
        }

        isComplete = true;
        return true;
    }

    std::string GetBestRoute() {
        if (!isComplete) return "Not complete";

        std::stringstream ss;
        for (int idx : bestPath) {
            ss << cities[idx].name << std::endl;
        }
        return ss.str();
    }

    int GetBestRouteLength() {
        if (!isComplete) return -1;
        return static_cast<int>(bestDistance);
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