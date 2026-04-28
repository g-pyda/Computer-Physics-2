#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include <tuple>
#include "engine.h"
#include "data_processing.h"

// g++ -std=c++23 main.cpp engine.cpp data_processing.cpp -o main.exe

int main() {
    bool relaxation_experiment = true;
    bool overrelaxation_experiment = true;
    bool grid_refinement_experiment = true;

    const int MAX_ITER = 1000;
    const double TOLERANCE = 1e-6;

    // ----------------------------------------------------
    // RELAXATION
    // experiment over the grid size for n in [2, 3, 4, 5, 6, 7]
    // ----------------------------------------------------
    if (relaxation_experiment) {
        std::cout << "Starting relaxation experiment..." << std::endl;

        std::vector<int> n_values = {2, 3, 4, 5, 6, 7};
        for (int n : n_values) {
            int N = (1 << n) + 1;  // 2^n + 1
            auto [grid, delta_x] = create_grid(N);
            auto start = std::chrono::high_resolution_clock::now();
            gauss_seidel("relaxation", 'r', grid, N, delta_x, 1.0, MAX_ITER, TOLERANCE);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "\nn = " << n << " | Elapsed time: " << std::fixed << std::setprecision(6) << elapsed.count() << " seconds" << std::endl;
            std::cout << "==========================================================" << std::endl;
        }
        std::vector<std::string> cases;
        for (int n : n_values) {
            cases.push_back(std::to_string(1 << n));
        }

        std::cout << "Relaxation experiment finished." << std::endl;
        std::cout << "==========================================================" << std::endl;
    }

    // ----------------------------------------------------
    // OVERRELAXATION - w
    // experiment for n = 5 over the weights w in range(0.5, 2.0) with step 0.1
    //
    // OVERRELAXATION - n
    // experiment for n = 1.5 over the sizes N in range(2^n + 1) for n in [2, 3, 4, 5, 6, 7]
    // ----------------------------------------------------
    if (overrelaxation_experiment) {
        std::cout << "\nStarting overrelaxation experiment..." << std::endl;

        std::vector<double> weights;
        for (double w = 0.5; w <= 2.0; w += 0.1) {
            weights.push_back(w);
        }
        int n = 5;
        int N = (1 << n) + 1;
        for (double w : weights) {
            auto [grid, delta_x] = create_grid(N);
            auto start = std::chrono::high_resolution_clock::now();
            gauss_seidel("overrelaxation-w", 'o', grid, N, delta_x, w, MAX_ITER, TOLERANCE);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "\nw = " << w << " | Elapsed time: " << std::fixed << std::setprecision(6) << elapsed.count() << " seconds" << std::endl;
            std::cout << "==========================================================" << std::endl;
        }

        std::vector<double> n_values = {2, 3, 4, 5, 6, 7};
        float best_weight = 1.5;
        for (int n : n_values) {
            int N = (1 << n) + 1;
            auto [grid, delta_x] = create_grid(N);
            auto start = std::chrono::high_resolution_clock::now();
            gauss_seidel("overrelaxation-n", 'r', grid, N, delta_x, best_weight, MAX_ITER, TOLERANCE);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "\nn = " << n << " | Elapsed time: " << std::fixed << std::setprecision(6) << elapsed.count() << " seconds" << std::endl;
            std::cout << "==========================================================" << std::endl;
        }

        std::cout << "Overrelaxation experiment finished." << std::endl;
        std::cout << "==========================================================" << std::endl;
    }

    // ----------------------------------------------------
    // GRID REFINEMENT
    // experiment for n = 3 over the grid refinement levels in range(0,5) with step 1
    // ----------------------------------------------------
    if (grid_refinement_experiment) {
        std::cout << "\nStarting grid refinement experiment..." << std::endl;

        int n = 3;
        int N = (1 << n) + 1;
        auto [grid, delta_x] = create_grid(N);
        for (int level = 0; level <= 5; ++level) {
            if (level > 0) {
                // refine the grid by doubling the number of points in each dimension
                int new_N = 2 * (N - 1) + 1;
                std::tie(grid, delta_x) = refine_grid(N, grid);
                N = new_N;
            }
            auto start = std::chrono::high_resolution_clock::now();
            gauss_seidel("grid_refinement", 'r',grid, N, delta_x, 1.5, MAX_ITER, TOLERANCE);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "\nRefinement level: n = " << N << " | Elapsed time: " << std::fixed << std::setprecision(6) << elapsed.count() << " seconds" << std::endl;
            std::cout << "==========================================================" << std::endl;
        }

        std::cout << "Grid refinement experiment finished." << std::endl;
        std::cout << "==========================================================" << std::endl;
    }


    return 0;
}
