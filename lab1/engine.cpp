#define _USE_MATH_DEFINES
#include "engine.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "data_processing.h"

// ------------------------------------------------------
// GRID CREATION AND MANIPULATION
// ------------------------------------------------------

std::pair<std::vector<std::vector<double>>, double> create_grid(int N) {
    std::vector<std::vector<double>> grid(N, std::vector<double>(N, 0.0));
    double delta_x = 3.0 / (N - 1);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == 0 || i == N - 1 || j == 0 || j == N - 1) {
                grid[i][j] = 0.0;
            } else {
                grid[i][j] = poisson_equation(i * delta_x, j * delta_x);
            }
        }
    }
    return {grid, delta_x};
}


std::pair<std::vector<std::vector<double>>, double> refine_grid(
    // allocating the new grid
    int prev_N, const std::vector<std::vector<double>>& coarse_grid) {
    int N = 2 * (prev_N - 1) + 1;
    std::vector<std::vector<double>> grid(N, std::vector<double>(N, 0.0));
    double delta_x = 3.0 / (N - 1);
    
    // filling with exact indices from previous grid
    std::vector<int> coarse_indices(prev_N);
    for (int i = 0; i < prev_N; ++i) {
        coarse_indices[i] = i;
    }
    for (int i : coarse_indices) {
        for (int j : coarse_indices) {
            grid[i*2][j*2] = coarse_grid[i][j];
        }
    }

    // filling the nodes INLINE with the exact values with their mean (2-point mean)
    for (int i = 0; i < N; ++i) {
        for (int j = (i + 1)%2; j < N; ++j) {
            if (i % 2 == 0) {
                grid[i][j] = 0.5 * (grid[i][j - 1] + grid[i][j + 1]);
            } else {
                grid[i][j] = 0.5 * (grid[i - 1][j] + grid[i + 1][j]);
            }
        }
    }

    // filling the nodes with the mean of their 4 neighbors (4-point mean)
    for (int i = 1; i < N - 1; i += 2) {
        for (int j = 1; j < N - 1; j += 2) {
            grid[i][j] = 0.25 * (grid[i - 1][j] + grid[i + 1][j] + grid[i][j - 1] + grid[i][j + 1]);
        }
    }

    return {grid, delta_x};
}

double poisson_equation(double x, double y) {
    return 2 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}

double residuum(const std::vector<std::vector<double>>& u, int i, int j, int N, double delta_x) {
    double res = (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]) / (delta_x * delta_x) +
                 poisson_equation(i * delta_x, j * delta_x);
    return res;
}

double energy_functional(const std::vector<std::vector<double>>& u, int N, double delta_x) {
    double energy = 0.0;
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            energy += 0.5 * std::pow((u[i + 1][j] - u[i - 1][j]) / (2 * delta_x), 2) +
                      0.5 * std::pow((u[i][j + 1] - u[i][j - 1]) / (2 * delta_x), 2) +
                      poisson_equation(i * delta_x, j * delta_x) * u[i][j];
        }
    }
    return energy * delta_x * delta_x;
}

void gauss_seidel(const std::string& dirname, const char study, std::vector<std::vector<double>>& u, int N, double delta_x, double weight, int max_iterations, double tolerance) {
    std::string study_dir;
    switch (study) {
        case 'r':
            study_dir = std::to_string(N);
            break;
        case 'o':
            study_dir = std::to_string(weight);
            break;
        default:
            study_dir = std::to_string(N);
            return;
    }
    cleanup_directory("./data/" + dirname + "/" + study_dir);
    const double energy_exact = -std::pow(3 * M_PI / 2, 2);
    float previous_energy = energy_functional(u, N, delta_x);
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        if (iteration % (max_iterations / 10) == 0) {
            std::cout << "\rIteration: " << iteration << "/" << max_iterations << std::flush;
        }
        // Update all grid points
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                u[i][j] = weight * 0.25 * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - delta_x * delta_x * poisson_equation(i * delta_x, j * delta_x)) +
                          (1 - weight) * u[i][j];
            }
        }
        // Calculate energy once per iteration
        double energy = energy_functional(u, N, delta_x);
        double diff = std::abs(energy - previous_energy);
        previous_energy = energy;
        if (diff < tolerance && iteration > 0) {
            std::cout << "\nConvergence achieved after " << iteration << " iterations." << std::endl;
            // saving the last residual
            for (int i = 1; i < N - 1; ++i) {
                for (int j = 1; j < N - 1; ++j) {
                    double res = residuum(u, i, j, N, delta_x);
                    write_to_file("./data/" + dirname + "/" + study_dir + "/residual_" + std::to_string(iteration) + ".txt", std::to_string(res) + " ");
                }
                write_to_file("./data/" + dirname + "/" + study_dir + "/residual_" + std::to_string(iteration) + ".txt", "\n");
            }
            break;
        }
        // saving residual every 10th iteration
        if (iteration % (max_iterations / 10) == 0) {
            for (int i = 1; i < N - 1; ++i) {
                for (int j = 1; j < N - 1; ++j) {
                    double res = residuum(u, i, j, N, delta_x);
                    write_to_file("./data/" + dirname + "/" + study_dir + "/residual_" + std::to_string(iteration) + ".txt", std::to_string(res) + " ");
                }
                write_to_file("./data/" + dirname + "/" + study_dir + "/residual_" + std::to_string(iteration) + ".txt", "\n");
            }
        }
        // saving energy functional
        write_to_file("./data/" + dirname + "/" + study_dir + "/energy_functional.txt", std::to_string(iteration) + " " + std::to_string(energy) + " " + std::to_string(energy_exact) + "\n");
        if (diff < tolerance) {
            break;
        }
    }
    save_grid_to_file("./data/" + dirname + "/" + study_dir + "/solution_grid.txt", u, N);
}
