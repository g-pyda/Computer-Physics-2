#define _USE_MATH_DEFINES
#include "engine.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <filesystem>
#include <lapacke.h>

// ========== STIFFNESS MATRICES (From Lab 2) ========== //

const std::vector<std::vector<double>> bilinear_stiffness_matrix = {
    { 0.666667, -0.166667, -0.333333, -0.166667},
    {-0.166667,  0.666667, -0.166667, -0.333333},
    {-0.333333, -0.166667,  0.666667, -0.166667},
    {-0.166667, -0.333333, -0.166667,  0.666667}
};

const std::vector<std::vector<double>> biparabolic_stiffness_matrix = {
    { 0.622222, -0.033333, -0.033333, -0.022222, -0.200000,  0.111111, -0.200000,  0.111111, -0.355556},
    {-0.033333,  0.622222, -0.022222, -0.033333, -0.200000, -0.200000,  0.111111,  0.111111, -0.355556},
    {-0.033333, -0.022222,  0.622222, -0.033333,  0.111111,  0.111111, -0.200000, -0.200000, -0.355556},
    {-0.022222, -0.033333, -0.033333,  0.622222,  0.111111, -0.200000,  0.111111, -0.200000, -0.355556},
    {-0.200000, -0.200000,  0.111111,  0.111111,  1.955556, -0.355556, -0.355556,  0.000000, -1.066667},
    { 0.111111, -0.200000,  0.111111, -0.200000, -0.355556,  1.955556,  0.000000, -0.355556, -1.066667},
    {-0.200000,  0.111111, -0.200000,  0.111111, -0.355556,  0.000000,  1.955556, -0.355556, -1.066667},
    { 0.111111,  0.111111, -0.200000, -0.200000,  0.000000, -0.355556, -0.355556,  1.955556, -1.066667},
    {-0.355556, -0.355556, -0.355556, -0.355556, -1.066667, -1.066667, -1.066667, -1.066667,  5.688889}
};

// ========== WAVE EQUATION EXPERIMENT ========== //

void run_experiment(const std::string& expName, double omega) {
    const double L = 5.0; 
    const double v = 1.0; 
    const int N = 41; 
    const double dx = L / (N - 1);
    
    const double dt = 0.05; 
    const int timeSteps = 2000;
    
    const int numNodes = N * N;
    const int centerNode = (N / 2) * N + (N / 2);

    std::string baseDir = "data/" + expName;
    std::filesystem::create_directories(baseDir);

    std::vector<double> S(numNodes * numNodes, 0.0);
    std::vector<double> O(numNodes * numNodes, 0.0);
    
    // Analytically computed Overlap matrix for bilinear elements
    double h2_36 = (dx * dx) / 36.0;
    const std::vector<std::vector<double>> bilinear_overlap_matrix = {
        { 4.0 * h2_36, 2.0 * h2_36, 1.0 * h2_36, 2.0 * h2_36 },
        { 2.0 * h2_36, 4.0 * h2_36, 2.0 * h2_36, 1.0 * h2_36 },
        { 1.0 * h2_36, 2.0 * h2_36, 4.0 * h2_36, 2.0 * h2_36 },
        { 2.0 * h2_36, 1.0 * h2_36, 2.0 * h2_36, 4.0 * h2_36 }
    };

    // Global matrix assembly for S and O
    for (int ej = 0; ej < N - 1; ++ej) {
        for (int ei = 0; ei < N - 1; ++ei) {
            int base = ej * N + ei;
            
            // CCW order mapping: Bottom-Left, Bottom-Right, Top-Right, Top-Left
            std::vector<int> nlg = { base, base + 1, base + N + 1, base + N };
            
            for (int i1 = 0; i1 < 4; ++i1) {
                for (int i2 = 0; i2 < 4; ++i2) {
                    int row = nlg[i1];
                    int col = nlg[i2];
                    S[row * numNodes + col] += bilinear_stiffness_matrix[i1][i2];
                    O[row * numNodes + col] += bilinear_overlap_matrix[i1][i2];
                }
            }
        }
    }

    // Apply Dirichlet Boundary Conditions directly to matrices
    for (int i = 0; i < numNodes; ++i) {
        int grid_x = i % N;
        int grid_y = i / N;
        bool is_boundary = (grid_x == 0 || grid_x == N - 1 || grid_y == 0 || grid_y == N - 1);
        bool is_center = (i == centerNode);

        if (is_boundary || is_center) {
            for (int j = 0; j < numNodes; ++j) {
                O[i * numNodes + j] = 0.0;
                S[i * numNodes + j] = 0.0; 
            }
            O[i * numNodes + i] = 1.0; 
        }
    }

    // LAPACK initialization & factorization
    std::vector<lapack_int> ipiv(numNodes);
    std::vector<double> O_factored = O;
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, numNodes, numNodes, O_factored.data(), numNodes, ipiv.data());

    std::vector<double> c_prev(numNodes, 0.0); 
    std::vector<double> c_curr(numNodes, 0.0); 
    std::vector<double> c_next(numNodes, 0.0);

    // Time loop
    for (int step = 1; step <= timeSteps; ++step) {
        std::cout << "\r  Time Step: " << std::setw(4) << step << "/" << timeSteps << std::flush;
        double t_next = step * dt;
        std::vector<double> rhs(numNodes, 0.0);

        for (int i = 0; i < numNodes; ++i) {
            int grid_x = i % N;
            int grid_y = i / N;
            bool is_boundary = (grid_x == 0 || grid_x == N - 1 || grid_y == 0 || grid_y == N - 1);

            if (is_boundary) {
                rhs[i] = 0.0; 
            } else if (i == centerNode) {
                rhs[i] = std::sin(omega * t_next); 
            } else {
                double s_term = 0.0;
                double o_curr_term = 0.0;
                double o_prev_term = 0.0;

                for (int j = 0; j < numNodes; ++j) {
                    s_term += S[i * numNodes + j] * c_curr[j];
                    o_curr_term += O[i * numNodes + j] * c_curr[j];
                    o_prev_term += O[i * numNodes + j] * c_prev[j];
                }
                
                rhs[i] = (-v * v * dt * dt * s_term + 2.0 * o_curr_term) - o_prev_term; 
            }
        }

        c_next = rhs;
        LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', numNodes, 1, O_factored.data(), numNodes, ipiv.data(), c_next.data(), 1);

        // Export data
        if (step % 10 == 0) {
            std::string filename = baseDir + "/step_" + std::to_string(step) + ".dat";
            std::ofstream out(filename);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    out << i * dx << " " << j * dx << " " << c_next[j * N + i] << "\n";
                }
                out << "\n";
            }
        }

        c_prev = c_curr;
        c_curr = c_next;
    }
}

void run_convergence_test(const std::string& expName, double omega) {
    const double L = 5.0; 
    const double v = 1.0; 
    const int N = 21; 
    const int numNodes = N * N;
    const int centerNode = (N / 2) * N + (N / 2);
    const double totalTime = 30.0;
    const double checkpointInterval = 0.1;

    std::string baseDir = "data/" + expName;
    std::filesystem::create_directories(baseDir);

    std::vector<double> S(numNodes * numNodes, 0.0);
    std::vector<double> O(numNodes * numNodes, 0.0);
    
    double dx = L / (N - 1);
    double h2_36 = (dx * dx) / 36.0;
    const std::vector<std::vector<double>> bilinear_overlap_matrix = {
        { 4.0 * h2_36, 2.0 * h2_36, 1.0 * h2_36, 2.0 * h2_36 },
        { 2.0 * h2_36, 4.0 * h2_36, 2.0 * h2_36, 1.0 * h2_36 },
        { 1.0 * h2_36, 2.0 * h2_36, 4.0 * h2_36, 2.0 * h2_36 },
        { 2.0 * h2_36, 1.0 * h2_36, 2.0 * h2_36, 4.0 * h2_36 }
    };

    // Global matrix assembly
    for (int ej = 0; ej < N - 1; ++ej) {
        for (int ei = 0; ei < N - 1; ++ei) {
            int base = ej * N + ei;
            std::vector<int> nlg = { base, base + 1, base + N + 1, base + N };
            for (int i1 = 0; i1 < 4; ++i1) {
                for (int i2 = 0; i2 < 4; ++i2) {
                    int row = nlg[i1];
                    int col = nlg[i2];
                    S[row * numNodes + col] += bilinear_stiffness_matrix[i1][i2];
                    O[row * numNodes + col] += bilinear_overlap_matrix[i1][i2];
                }
            }
        }
    }

    // Boundary Conditions
    for (int i = 0; i < numNodes; ++i) {
        int grid_x = i % N;
        int grid_y = i / N;
        if (grid_x == 0 || grid_x == N - 1 || grid_y == 0 || grid_y == N - 1 || i == centerNode) {
            for (int j = 0; j < numNodes; ++j) {
                O[i * numNodes + j] = 0.0;
                S[i * numNodes + j] = 0.0; 
            }
            O[i * numNodes + i] = 1.0; 
        }
    }

    // Matrix Factorization (Done once for all runs since O doesn't change)
    std::vector<lapack_int> ipiv(numNodes);
    std::vector<double> O_factored = O;
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, numNodes, numNodes, O_factored.data(), numNodes, ipiv.data());

    // Helper lambda to run a simulation and either store or compare results
    auto runSim = [&](double dt, std::vector<std::vector<double>>& refData, bool isReference) {
        std::vector<double> c_prev(numNodes, 0.0); 
        std::vector<double> c_curr(numNodes, 0.0); 
        std::vector<double> c_next(numNodes, 0.0);
        
        int timeSteps = std::round(totalTime / dt);
        double nextCheckpoint = checkpointInterval;
        int refIdx = 0;
        
        std::ofstream outFile;
        if (!isReference) {
            std::string filename = baseDir + "/error_dt_" + std::to_string(dt) + ".dat";
            outFile.open(filename);
            outFile << "time error\n";
        }

        for (int step = 1; step <= timeSteps; ++step) {
            double t_next = step * dt;
            std::vector<double> rhs(numNodes, 0.0);

            for (int i = 0; i < numNodes; ++i) {
                int grid_x = i % N;
                int grid_y = i / N;
                if (grid_x == 0 || grid_x == N - 1 || grid_y == 0 || grid_y == N - 1) {
                    rhs[i] = 0.0; 
                } else if (i == centerNode) {
                    rhs[i] = std::sin(omega * t_next); 
                } else {
                    double s_term = 0.0, o_curr_term = 0.0, o_prev_term = 0.0;
                    for (int j = 0; j < numNodes; ++j) {
                        s_term += S[i * numNodes + j] * c_curr[j];
                        o_curr_term += O[i * numNodes + j] * c_curr[j];
                        o_prev_term += O[i * numNodes + j] * c_prev[j];
                    }
                    rhs[i] = (-v * v * dt * dt * s_term + 2.0 * o_curr_term) - o_prev_term; 
                }
            }

            c_next = rhs;
            LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', numNodes, 1, O_factored.data(), numNodes, ipiv.data(), c_next.data(), 1);

            // Checkpoint saving / Error calculation
            if (t_next >= nextCheckpoint - 1e-9) {
                if (isReference) {
                    refData.push_back(c_next);
                } else {
                    double mse = 0.0;
                    for (int i = 0; i < numNodes; ++i) {
                        double diff = c_next[i] - refData[refIdx][i];
                        mse += diff * diff;
                    }
                    double rmse = std::sqrt(mse / numNodes);
                    outFile << t_next << " " << rmse << "\n";
                    refIdx++;
                }
                nextCheckpoint += checkpointInterval;
            }

            c_prev = c_curr;
            c_curr = c_next;
        }
    };

    // 1. Generate Reference Solution
    std::cout << "  Generating reference solution (dt = 0.0005)..." << std::endl;
    std::vector<std::vector<double>> reference_states;
    runSim(0.0005, reference_states, true);

    // 2. Test specific dt values
    std::vector<double> test_dts = {0.1, 0.05, 0.025, 0.01, 0.005};
    for (double test_dt : test_dts) {
        std::cout << "  Testing dt = " << test_dt << "..." << std::endl;
        runSim(test_dt, reference_states, false);
    }
}