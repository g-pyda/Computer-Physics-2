#include "engine.h"
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <lapacke.h>

// ========== STIFFNESS MATRICES AND COMPONENTS ========== //

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

const std::vector<double> w = {
    (18.0 - std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 - std::sqrt(30.0)) / 36.0
};

const std::vector<double> gamma_pts = {
    -std::sqrt(3.0/7.0 + 2.0/7.0*std::sqrt(6.0/5.0)),
    -std::sqrt(3.0/7.0 - 2.0/7.0*std::sqrt(6.0/5.0)),
     std::sqrt(3.0/7.0 - 2.0/7.0*std::sqrt(6.0/5.0)),
     std::sqrt(3.0/7.0 + 2.0/7.0*std::sqrt(6.0/5.0))
};
// const std::vector<double> gamma = {
//     -1*sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)),
//     sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)),
//     sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0)),
//     -1*sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0))
// }


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

void solve_and_save(
    std::vector<std::vector<double>>& H,
    std::vector<std::vector<double>>& O,
    std::vector<double>& E,
    const std::string& path,
    int total_nodes
) {
    int n = total_nodes;
    int nrhs = 1;
    int lda = n;
    int ldb = n;

    // creating flat buffers for LAPACK
    std::vector<double> H_flat(n * n);
    std::vector<double> O_flat(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            H_flat[i * n + j] = H[i][j];
            O_flat[i * n + j] = O[i][j];
        }
    }

    // solving generalized symmetric-definite eigenproblem (eigenvalues in E, eigenvectors in H_flat)
    lapack_int info = LAPACKE_dsygv(
        LAPACK_ROW_MAJOR, 
        1,
        'V',
        'U',
        n,
        H_flat.data(),
        lda,
        O_flat.data(),
        ldb,
        E.data()
    );

    if (info > 0) {
        std::cerr << "error: matrix is singular at U(" << info << "," << info << ")" << std::endl;
        return;
    } else if (info < 0) {
        std::cerr << "error: illegal argument at position " << -info << " in lapack call" << std::endl;
        return;
    }

    // saving eigenvalues to eigenvalues.txt
    std::ofstream f_sol(path + "/eigenvalues.txt");
    if (!f_sol.is_open()) {
        std::cerr << "error: could not open solution file for writing" << std::endl;
        return;
    }

    for (int i = 0; i < n; i++) {
        f_sol << i + 1 << " " << std::fixed << std::setprecision(10) << E[i] << "\n";
    }
    f_sol.close();

    // saving eigenvectors to eigenvectors.txt
    std::ofstream f_vec(path + "/eigenvectors.txt");
    if (!f_vec.is_open()) {
        std::cerr << "error: could not open eigenvectors file for writing" << std::endl;
        return;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            f_vec << std::fixed << std::setprecision(10) << H_flat[i*n + j] << " ";
        f_vec << "\n";
    }
    f_vec.close();
}

// ========== MESH AND GLOBAL ASSEMBLY ========== //

void run_simulation(int N, double L, bool biparabolic) {
    double a = L / N;

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

    // assembling global stiffness matrix H and overlap matric O
    std::vector<std::vector<double>> H(total_nodes, std::vector<double>(total_nodes, 0.0));
    std::vector<std::vector<double>> O(total_nodes, std::vector<double>(total_nodes, 0.0));
    std::vector<double> E(total_nodes, 0.0);
    std::vector<double> c(total_nodes, 0.0);
    const auto& local_S = biparabolic ? biparabolic_stiffness_matrix : bilinear_stiffness_matrix;

    for (int k = 0; k < (int)nlg.size(); ++k) {
        int m = (int)local_S.size();
        // assemble stiffness (H) as before
        for (int i1 = 0; i1 < m; ++i1) {
            for (int i2 = 0; i2 < m; ++i2) {
                H[nlg[k][i1]][nlg[k][i2]] += 0.5 * local_S[i1][i2];
            }
        }

        // assemble overlap (mass) matrix using 2D Gauss quadrature on reference square [-1,1]^2
        std::vector<std::vector<double>> local_O(m, std::vector<double>(m, 0.0));
        for (int gi = 0; gi < 4; ++gi) {
            for (int gj = 0; gj < 4; ++gj) {
                double wt = w[gi] * w[gj];
                for (int i1 = 0; i1 < m; ++i1) {
                    double phi_i = (m == 4) ? g_func(i1 + 1, gamma_pts[gi], gamma_pts[gj]) : h_func(i1 + 1, gamma_pts[gi], gamma_pts[gj]);
                    for (int i2 = 0; i2 < m; ++i2) {
                        double phi_j = (m == 4) ? g_func(i2 + 1, gamma_pts[gi], gamma_pts[gj]) : h_func(i2 + 1, gamma_pts[gi], gamma_pts[gj]);
                        local_O[i1][i2] += wt * phi_i * phi_j;
                    }
                }
            }
        }
        // account for Jacobian of mapping from reference to physical element: |J| = (a/2)*(a/2) = a*a/4
        double jac = a * a / 4.0;
        for (int i1 = 0; i1 < m; ++i1) {
            for (int i2 = 0; i2 < m; ++i2) {
                O[nlg[k][i1]][nlg[k][i2]] += jac * local_O[i1][i2];
            }
        }
    }

    // applying boundary conditions
    double eps = 1e-9;
    for (int i = 0; i < total_nodes; ++i) {
        if (p[i].x < eps || p[i].x > L - eps || p[i].y < eps || p[i].y > L - eps) {
            for (int j = 0; j < total_nodes; ++j)
                O[i][j] = O[j][i] = H[i][j] = H[j][i] = 0.0;
            
            H[i][i] = -1.0;
            O[i][i] = 1.0;
        }
    }

    // finalizing simulation and solving
    std::string folder = biparabolic ? "biparabolic" : "bilinear";
    std::string path = "./data/" + folder + "/" + std::to_string(N);
    std::filesystem::create_directories(path);

    // solving Hc = EOc and saving solution.txt
    solve_and_save(H, O, E, path, total_nodes);

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