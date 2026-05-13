#include "engine.h"
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <lapacke.h>
#include <complex>
#include <algorithm>

typedef std::complex<double> complexd;

const double h_bar = 1.0;
const double dt = 0.25;      
const double T_max = 500.0;
const int N = 16;
const double L = 5.0;

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

// Solve generalized eigenvalue problem returning eigenvectors safely
std::vector<std::vector<double>> solve_eigen(
    const std::vector<std::vector<double>>& H,
    const std::vector<std::vector<double>>& O,
    std::vector<double>& E,
    int total_nodes
) {
    int n = total_nodes;
    std::vector<double> H_flat(n * n);
    std::vector<double> O_flat(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            H_flat[i * n + j] = H[i][j];
            O_flat[i * n + j] = O[i][j];
        }
    }

    lapack_int info = LAPACKE_dsygv(
        LAPACK_ROW_MAJOR, 1, 'V', 'U', n,
        H_flat.data(), n, O_flat.data(), n, E.data()
    );

    std::vector<std::vector<double>> eigenvecs(n, std::vector<double>(n, 0.0));
    if (info == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                eigenvecs[i][j] = H_flat[i * n + j]; 
            }
        }
    }
    return eigenvecs;
}

// Assemble external potential matrix V
void assemble_potential_matrix(
    int total_nodes,
    int N_elem, double L_total,
    const std::vector<std::vector<int>>& nlg,
    std::vector<std::vector<double>>& V
) {
    double a = L_total / N_elem;
    double jac = a * a / 4.0;

    for (int k = 0; k < (int)nlg.size(); ++k) {
        for (int gi = 0; gi < 4; ++gi) {
            for (int gj = 0; gj < 4; ++gj) {
                double wt = w[gi] * w[gj];
                double xi1 = gamma_pts[gi];
                double xi2 = gamma_pts[gj];

                double phys_x = (k % N_elem) * a + (1.0 + xi1) * a / 2.0;

                for (int i1 = 0; i1 < 9; ++i1) {
                    int row = nlg[k][i1];
                    double phi_i = h_func(i1 + 1, xi1, xi2);
                    for (int i2 = 0; i2 < 9; ++i2) {
                        int col = nlg[k][i2];
                        double phi_j = h_func(i2 + 1, xi1, xi2);
                        V[row][col] += jac * wt * phys_x * phi_i * phi_j;
                    }
                }
            }
        }
    }
}

// Perform Crank-Nicolson time evolution
void time_evolution(int total_nodes, std::vector<double>& c_initial, 
                    std::vector<std::vector<double>>& O,
                    std::vector<std::vector<double>>& H0, 
                    std::vector<std::vector<double>>& V,
                    const std::vector<Node>& nodes_coords,
                    std::string path,
                    double omega,
                    double eF,
                    std::vector<double>& psi0,
                    std::vector<double>& psi1,
                    std::vector<double>& psi2
                ) {
    
    int n = total_nodes;
    double C_const = dt / (2.0 * h_bar); 
    complexd i_comp(0.0, 1.0); 
    
    // Apply Dirichlet boundary conditions to stabilize Crank-Nicolson scheme
    double eps = 1e-9;
    for (int i = 0; i < n; ++i) {
        bool is_boundary = (nodes_coords[i].x < eps || nodes_coords[i].x > L - eps || 
                            nodes_coords[i].y < eps || nodes_coords[i].y > L - eps);
        if (is_boundary) {
            for (int j = 0; j < n; ++j) {
                H0[i][j] = 0.0; H0[j][i] = 0.0;
                V[i][j] = 0.0;  V[j][i] = 0.0;
                O[i][j] = 0.0;  O[j][i] = 0.0;
            }
            O[i][i] = 1.0; 
        }
    }

    std::vector<complexd> A_flat(n * n, 0.0);
    std::vector<complexd> B_flat(n, 0.0);
    std::vector<complexd> c(n);
    for (int i = 0; i < n; ++i) c[i] = complexd(c_initial[i], 0.0);
    std::vector<int> ipiv(n);

    std::ofstream f_out(path + "projections.txt");
    
    for (double t = 0; t <= T_max; t += dt) {
        std::fill(A_flat.begin(), A_flat.end(), complexd(0.0, 0.0));
        std::fill(B_flat.begin(), B_flat.end(), complexd(0.0, 0.0));
        
        double st = std::sin(omega * t);
        double st_dt = std::sin(omega * (t + dt));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                complexd h_t = complexd(H0[i][j] + eF * V[i][j] * st, 0.0);
                complexd h_t_dt = complexd(H0[i][j] + eF * V[i][j] * st_dt, 0.0);

                // LHS matrix A
                A_flat[i * n + j] = complexd(O[i][j], 0.0) + i_comp * C_const * h_t_dt;

                // RHS vector B
                complexd m_right_ij = complexd(O[i][j], 0.0) - i_comp * C_const * h_t;
                B_flat[i] += m_right_ij * c[j];
            }
        }
        
        lapack_int info = LAPACKE_zgesv(
            LAPACK_ROW_MAJOR, n, 1,
            (lapack_complex_double*)A_flat.data(), n,
            ipiv.data(),
            (lapack_complex_double*)B_flat.data(), 1
        );

        if (info != 0) break;

        c = B_flat;

        // Compute projections using complex conjugate
        complexd p0 = 0.0, p1 = 0.0, p2 = 0.0, norm_sq = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                complexd o_ij = complexd(O[i][j], 0.0);
                norm_sq += std::conj(c[i]) * o_ij * c[j]; 
                p0 += complexd(psi0[i], 0.0) * o_ij * c[j]; 
                p1 += complexd(psi1[i], 0.0) * o_ij * c[j];
                p2 += complexd(psi2[i], 0.0) * o_ij * c[j];
            }
        }

        // Calculate leakage as probability of being in states higher than p0, p1, p2
        double current_norm = std::real(norm_sq);
        double prob_sum = std::norm(p0) + std::norm(p1) + std::norm(p2);
        double leakage = current_norm - prob_sum;

        f_out << t << " " << std::norm(p0) << " " << std::norm(p1) << " " 
              << std::norm(p2) << " " << leakage << "\n";
    }
    f_out.close();
}

double get_ground_state(
    std::vector<std::vector<double>>& H,
    std::vector<std::vector<double>>& O,
    std::vector<std::vector<int>>& nlg,
    std::vector<double>& c,
    std::vector<double>& psi0,
    std::vector<double>& psi1,
    std::vector<double>& psi2,
    int N_elem, double L_total, bool biparabolic
) {
    double a = L_total / N_elem;
    int nodes_per_side = biparabolic ? (2 * N_elem + 1) : (N_elem + 1);
    int total_nodes = nodes_per_side * nodes_per_side;
    double h_step = L_total / (nodes_per_side - 1);

    std::vector<Node> p(total_nodes);
    for (int j = 0; j < nodes_per_side; ++j) {
        for (int i = 0; i < nodes_per_side; ++i) {
            p[j * nodes_per_side + i] = { i * h_step, j * h_step };
        }
    }

    int step = biparabolic ? 2 : 1;
    for (int ej = 0; ej < N_elem; ++ej) {
        for (int ei = 0; ei < N_elem; ++ei) {
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

    std::vector<double> E(total_nodes, 0.0);
    const auto& local_S = biparabolic ? biparabolic_stiffness_matrix : bilinear_stiffness_matrix;
    int m = (int)local_S.size();

    for (int k = 0; k < (int)nlg.size(); ++k) {
        for (int i1 = 0; i1 < m; ++i1) {
            for (int i2 = 0; i2 < m; ++i2) {
                H[nlg[k][i1]][nlg[k][i2]] += 0.5 * local_S[i1][i2];
            }
        }
        double jac = a * a / 4.0;
        for (int gi = 0; gi < 4; ++gi) {
            for (int gj = 0; gj < 4; ++gj) {
                double wt = w[gi] * w[gj];
                for (int i1 = 0; i1 < m; ++i1) {
                    double phi_i = (m == 4) ? g_func(i1 + 1, gamma_pts[gi], gamma_pts[gj]) : h_func(i1 + 1, gamma_pts[gi], gamma_pts[gj]);
                    for (int i2 = 0; i2 < m; ++i2) {
                        double phi_j = (m == 4) ? g_func(i2 + 1, gamma_pts[gi], gamma_pts[gj]) : h_func(i2 + 1, gamma_pts[gi], gamma_pts[gj]);
                        O[nlg[k][i1]][nlg[k][i2]] += jac * wt * phi_i * phi_j;
                    }
                }
            }
        }
    }

    double b_eps = 1e-9;
    for (int i = 0; i < total_nodes; ++i) {
        if (p[i].x < b_eps || p[i].x > L_total - b_eps || p[i].y < b_eps || p[i].y > L_total - b_eps) {
            for (int j = 0; j < total_nodes; ++j) O[i][j] = O[j][i] = H[i][j] = H[j][i] = 0.0;
            H[i][i] = 1e12; 
            O[i][i] = 1.0;
        }
    }

    std::vector<std::vector<double>> evecs = solve_eigen(H, O, E, total_nodes);
    
    // Sort eigenvalues to identify physical states
    std::vector<std::pair<double, int>> indexed_E;
    for(int i=0; i<total_nodes; ++i) indexed_E.push_back({E[i], i});
    std::sort(indexed_E.begin(), indexed_E.end());

    int idx0 = indexed_E[0].second;
    int idx1 = indexed_E[1].second;
    int idx2 = indexed_E[2].second;

    for (int i = 0; i < total_nodes; ++i) {
        psi0[i] = evecs[i][idx0];
        psi1[i] = evecs[i][idx1];
        psi2[i] = evecs[i][idx2];
    }

    // Normalize initial state w.r.t. Overlap matrix
    auto normalize = [&](std::vector<double>& psi) {
        double norm2 = 0.0;
        for (int i = 0; i < total_nodes; ++i) {
            for (int j = 0; j < total_nodes; ++j) {
                norm2 += psi[i] * O[i][j] * psi[j];
            }
        }
        double factor = 1.0 / std::sqrt(norm2);
        for (int i = 0; i < total_nodes; ++i) psi[i] *= factor;
    };

    normalize(psi0);
    normalize(psi1);
    normalize(psi2);
    
    for (int i = 0; i < total_nodes; ++i) c[i] = psi0[i];

    // Clear boundary penalties from Hamiltonian before time evolution
    for (int i = 0; i < total_nodes; ++i) {
        if (p[i].x < b_eps || p[i].x > L_total - b_eps || p[i].y < b_eps || p[i].y > L_total - b_eps) {
            H[i][i] = 0.0; 
        }
    }

    return indexed_E[1].first - indexed_E[0].first;
}

// Main execution function
void run_simulation(std::string path, bool half_omega, double eF) {
    std::string folder = "./data/" + path + "/";
    std::filesystem::create_directories(folder);

    int nodes_per_side = 2 * N + 1;
    int total_nodes = nodes_per_side * nodes_per_side;
    double h_step = L / (nodes_per_side - 1);

    std::vector<Node> nodes_coords(total_nodes);
    for (int j = 0; j < nodes_per_side; ++j) {
        for (int i = 0; i < nodes_per_side; ++i) {
            nodes_coords[j * nodes_per_side + i] = { i * h_step, j * h_step };
        }
    }

    std::vector<std::vector<double>> H0(total_nodes, std::vector<double>(total_nodes, 0.0));
    std::vector<std::vector<double>> O(total_nodes, std::vector<double>(total_nodes, 0.0));
    std::vector<std::vector<double>> V(total_nodes, std::vector<double>(total_nodes, 0.0));
    std::vector<double> c(total_nodes, 0.0);
    std::vector<std::vector<int>> nlg;
    
    std::vector<double> psi0(total_nodes), psi1(total_nodes), psi2(total_nodes);

    double dE = get_ground_state(H0, O, nlg, c, psi0, psi1, psi2, N, L, true);

    double omega_res = dE / h_bar;
    if (half_omega) omega_res /= 2.0;

    assemble_potential_matrix(total_nodes, N, L, nlg, V);
    time_evolution(total_nodes, c, O, H0, V, nodes_coords, folder, omega_res, eF, psi0, psi1, psi2);
}