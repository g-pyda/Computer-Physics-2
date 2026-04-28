#include "engine.h"
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <lapacke.h>

// ========== STIFFNESS MATRICES ========== //

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

// ========== SHAPE FUNCTIONS ========== //

double g_func(int i, double xi1, double xi2) {
    auto f1 = [](double xi) { return (1.0 - xi) / 2.0; };
    auto f2 = [](double xi) { return (1.0 + xi) / 2.0; };
    if (i == 1) return f1(xi1) * f1(xi2);
    if (i == 2) return f2(xi1) * f1(xi2);
    if (i == 3) return f1(xi1) * f2(xi2);
    if (i == 4) return f2(xi1) * f2(xi2);
    return 0;
}

double h_func(int i, double xi1, double xi2) {
    auto q1 = [](double xi) { return xi * (xi - 1.0) / 2.0; };
    auto q2 = [](double xi) { return (1.0 - xi) * (1.0 + xi); };
    auto q3 = [](double xi) { return xi * (xi + 1.0) / 2.0; };
    switch(i) {
        case 1: return q1(xi1) * q1(xi2); case 2: return q3(xi1) * q1(xi2);
        case 3: return q1(xi1) * q3(xi2); case 4: return q3(xi1) * q3(xi2);
        case 5: return q2(xi1) * q1(xi2); case 6: return q3(xi1) * q2(xi2);
        case 7: return q1(xi1) * q2(xi2); case 8: return q2(xi1) * q3(xi2);
        case 9: return q2(xi1) * q2(xi2);
        default: return 0;
    }
}

// ========== EXACT SOLUTION AND BC ========== //

double get_exact_potential(double x, double y) {
    double xc = -1.0, yc = 0.0;
    double r = std::sqrt(std::pow(x - xc, 2) + std::pow(y - yc, 2));
    return -1.0 / (2.0 * M_PI) * std::log(r);
}

// ========== LINEAR SYSTEM SOLVER (LAPACKE) ========== //

void solve_and_save(const std::vector<std::vector<double>>& S, std::vector<double>& b, const std::string& path, int total_nodes) {
    int n = total_nodes;
    int nrhs = 1;
    int lda = n;
    int ldb = 1;

    // creating a flat buffer for LAPACK
    std::vector<double> a_flat(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a_flat[i * n + j] = S[i][j];
        }
    }

    std::vector<lapack_int> ipiv(n);

    // solving Sc = b
    lapack_int info = LAPACKE_dgesv(
        LAPACK_ROW_MAJOR, 
        n, 
        nrhs, 
        a_flat.data(), 
        lda, 
        ipiv.data(), 
        b.data(), 
        nrhs
    );

    if (info > 0) {
        std::cerr << "error: matrix is singular at U(" << info << "," << info << ")" << std::endl;
        return;
    } else if (info < 0) {
        std::cerr << "error: illegal argument at position " << -info << " in lapack call" << std::endl;
        return;
    }

    // saving coefficients to solution.txt
    std::ofstream f_sol(path + "/solution.txt");
    if (!f_sol.is_open()) {
        std::cerr << "error: could not open solution file for writing" << std::endl;
        return;
    }

    for (int i = 0; i < n; i++) {
        f_sol << i + 1 << " " << std::fixed << std::setprecision(10) << b[i] << "\n";
    }
    f_sol.close();
}

// ========== MESH AND GLOBAL ASSEMBLY ========== //

void run_simulation(int N, double L, bool biparabolic) {
    // calculating grid dimensions
    int nodes_per_side = biparabolic ? (2 * N + 1) : (N + 1);
    int total_nodes = nodes_per_side * nodes_per_side;
    double h_step = L / (nodes_per_side - 1);

    // generating global node coordinates
    std::vector<Node> p(total_nodes);
    for (int j = 0; j < nodes_per_side; ++j) {
        for (int i = 0; i < nodes_per_side; ++i) {
            p[j * nodes_per_side + i] = { i * h_step, j * h_step };
        }
    }

    // building mapping table nlg
    std::vector<std::vector<int>> nlg;
    int step = biparabolic ? 2 : 1;
    for (int ej = 0; ej < N; ++ej) {
        for (int ei = 0; ei < N; ++ei) {
            int base = ej * step * nodes_per_side + ei * step;
            if (!biparabolic) {
                nlg.push_back({ base, base + 1, base + nodes_per_side, base + nodes_per_side + 1 });
            } else {
                nlg.push_back({
                    base, base + 2, base + 2 * nodes_per_side, base + 2 * nodes_per_side + 2,
                    base + 1, base + nodes_per_side + 2, base + nodes_per_side, base + 2 * nodes_per_side + 1,
                    base + nodes_per_side + 1
                });
            }
        }
    }

    // assembling global stiffness matrix
    std::vector<std::vector<double>> S(total_nodes, std::vector<double>(total_nodes, 0.0));
    std::vector<double> b(total_nodes, 0.0);
    const auto& local_S = biparabolic ? biparabolic_stiffness_matrix : bilinear_stiffness_matrix;

    for (int k = 0; k < (int)nlg.size(); ++k) {
        for (int i1 = 0; i1 < (int)local_S.size(); ++i1) {
            for (int i2 = 0; i2 < (int)local_S.size(); ++i2) {
                S[nlg[k][i1]][nlg[k][i2]] += local_S[i1][i2];
            }
        }
    }

    // applying dirichlet boundary conditions
    double eps = 1e-9;
    for (int i = 0; i < total_nodes; ++i) {
        if (p[i].x < eps || p[i].x > L - eps || p[i].y < eps || p[i].y > L - eps) {
            std::fill(S[i].begin(), S[i].end(), 0.0);
            S[i][i] = 1.0;
            b[i] = get_exact_potential(p[i].x, p[i].y);
        }
    }

    // finalizing simulation and solving
    std::string folder = biparabolic ? "biparabolic" : "bilinear";
    std::string path = "./data/" + folder + "/" + std::to_string(N);
    std::filesystem::create_directories(path);

    // solving Sc = b and saving solution.txt
    solve_and_save(S, b, path, total_nodes);

    // exporting mapping and nodes
    std::ofstream f_nodes(path + "/nodes.txt");
    for (int i = 0; i < total_nodes; i++) f_nodes << i + 1 << " " << p[i].x << " " << p[i].y << "\n";

    std::ofstream f_nlg(path + "/nlg.txt");
    for (int k = 0; k < (int)nlg.size(); ++k) {
        f_nlg << k + 1 << " ";
        for (int node : nlg[k]) f_nlg << node + 1 << " ";
        f_nlg << "\n";
    }
}