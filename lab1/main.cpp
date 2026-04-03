#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include "engine.h"
#include "data_processing.h"

// g++ -std=c++23 main.cpp engine.cpp data_processing.cpp -o main.exe

int main() {
    bool relaxation_experiment = false;
    bool overrelaxation_experiment = true;
    bool grid_refinement_experiment = false;

    const int MAX_ITER = 300;
    const double TOLERANCE = 1e-6;

    // ----------------------------------------------------
    // RELAXATION
    // experiment over the grid size for n in [2, 3, 4, 5]
    // ----------------------------------------------------
    if (relaxation_experiment) {
        std::cout << "Starting relaxation experiment..." << std::endl;

        std::vector<int> n_values = {2, 3, 4, 5, 6};
        for (int n : n_values) {
            int N = (1 << n) + 1;  // 2^n + 1
            auto [grid, delta_x] = create_grid(N);
            auto start = std::chrono::high_resolution_clock::now();
            gauss_seidel("relaxation", 'r', grid, N, delta_x, 1.0, MAX_ITER, TOLERANCE);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "\nn = " << n << " | Elapsed time: " << std::fixed << std::setprecision(6) << elapsed.count() << " seconds" << std::endl;
        }
        std::vector<std::string> cases;
        for (int n : n_values) {
            cases.push_back(std::to_string(1 << n));
        }

        std::cout << "Relaxation experiment finished." << std::endl;
    }

    // ----------------------------------------------------
    // OVERRELAXATION
    // experiment for n = 4 over the weights w in [0.5, 1.0, 1.5, 1.9]
    // ----------------------------------------------------
    if (overrelaxation_experiment) {
        std::cout << "\nStarting overrelaxation experiment..." << std::endl;

        std::vector<double> weights = {0.5, 1.0, 1.5, 1.9};
        int n = 5;
        int N = (1 << n) + 1;
        auto [grid, delta_x] = create_grid(N);
        for (double w : weights) {
            auto start = std::chrono::high_resolution_clock::now();
            gauss_seidel("overrelaxation", 'o', grid, N, delta_x, w, MAX_ITER, TOLERANCE);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "\nw = " << w << " | Elapsed time: " << std::fixed << std::setprecision(6) << elapsed.count() << " seconds" << std::endl;
        }
    }


    return 0;
}
