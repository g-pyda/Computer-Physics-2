#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <lapacke.h>

// Gauss quadrature weights and points for numerical integration
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

// Bilinear shape functions defined on the reference square [-1, 1]^2
double g_func(int i, double xi1, double xi2) {
    auto f1 = [](double xi) { return (1.0 - xi) / 2.0; };
    auto f2 = [](double xi) { return (1.0 + xi) / 2.0; };
    if (i == 1) return f1(xi1) * f1(xi2);
    if (i == 2) return f2(xi1) * f1(xi2);
    if (i == 3) return f1(xi1) * f2(xi2);
    if (i == 4) return f2(xi1) * f2(xi2);
    return 0;
}

// Numerical derivative using finite differences: (a/2) * [g_i(p + delta) - g_i(p - delta)] / (2*delta)
double dg_dxi(int i, double xi1, double xi2, int direction) {
    double delta = 1e-5;
    if (direction == 1) {
        return (g_func(i, xi1 + delta, xi2) - g_func(i, xi1 - delta, xi2)) / (2.0 * delta);
    } else {
        return (g_func(i, xi1, xi2 + delta) - g_func(i, xi1, xi2 - delta)) / (2.0 * delta);
    }
}

// Main function encapsulating the full simulation logic
void run_simulation(double D, double vx, double vy, int N, const std::string& mode_name) {
    std::cout << "Starting simulation: " << mode_name << "\n";
    
    double L = 5.0; 
    double dt = 0.02; 
    double a = L / N;
    int total_nodes = N * N; // For periodic boundary conditions, total nodes = total elements
    
    // Element-to-global node mapping table (nlg) implementing periodic boundary conditions
    std::vector<std::vector<int>> nlg;
    for (int ej = 0; ej < N; ++ej) {
        for (int ei = 0; ei < N; ++ei) {
            int n1 = ej * N + ei;
            int n2 = ej * N + (ei + 1) % N;
            int n3 = ((ej + 1) % N) * N + ei;
            int n4 = ((ej + 1) % N) * N + (ei + 1) % N;
            nlg.push_back({n1, n2, n3, n4});
        }
    }
    
    // Matrix initialization (flat arrays for LAPACK compatibility)
    std::vector<double> O_flat(total_nodes * total_nodes, 0.0);
    std::vector<double> R_flat(total_nodes * total_nodes, 0.0); // Combined operator matrix: R = -C - D*S

    // Global matrix assembly using 2D Gauss quadrature
    for (int k = 0; k < (int)nlg.size(); ++k) {
        for (int gi = 0; gi < 4; ++gi) {
            for (int gj = 0; gj < 4; ++gj) {
                double wt = w[gi] * w[gj];
                double xi1 = gamma_pts[gi];
                double xi2 = gamma_pts[gj];
                
                for (int i1 = 0; i1 < 4; ++i1) {
                    double phi_i = g_func(i1 + 1, xi1, xi2);
                    double dxi1_i = dg_dxi(i1 + 1, xi1, xi2, 1);
                    double dxi2_i = dg_dxi(i1 + 1, xi1, xi2, 2);
                    
                    for (int i2 = 0; i2 < 4; ++i2) {
                        double phi_j = g_func(i2 + 1, xi1, xi2);
                        double dxi1_j = dg_dxi(i2 + 1, xi1, xi2, 1);
                        double dxi2_j = dg_dxi(i2 + 1, xi1, xi2, 2);
                        
                        // Overlap (mass) matrix component O
                        double o_val = wt * phi_i * phi_j * (a * a / 4.0);
                        
                        // Stiffness (diffusion) matrix component S
                        double s_val = wt * (dxi1_i * dxi1_j + dxi2_i * dxi2_j); 
                        
                        // Advection matrix component C
                        double c_val = wt * phi_i * (vx * dxi1_j + vy * dxi2_j) * (a / 2.0); 
                        
                        int row = nlg[k][i1];
                        int col = nlg[k][i2];
                        int idx = row * total_nodes + col;
                        
                        O_flat[idx] += o_val;
                        R_flat[idx] += (-c_val - D * s_val);
                    }
                }
            }
        }
    }
    
    // Construct LHS and RHS matrices for the Crank-Nicolson scheme
    std::vector<double> M_LHS(total_nodes * total_nodes, 0.0);
    std::vector<double> M_RHS(total_nodes * total_nodes, 0.0);
    for (int i = 0; i < total_nodes * total_nodes; ++i) {
        M_LHS[i] = O_flat[i] - (dt / 2.0) * R_flat[i]; 
        M_RHS[i] = O_flat[i] + (dt / 2.0) * R_flat[i]; 
    }
    
    // Set up the Gaussian initial condition profile
    std::vector<double> c_vec(total_nodes, 0.0);
    for (int ej = 0; ej < N; ++ej) {
        for (int ei = 0; ei < N; ++ei) {
            double x = ei * a;
            double y = ej * a;
            c_vec[ej * N + ei] = std::exp(-2.0 * std::pow(x - L/2.0, 2) - 2.0 * std::pow(y - L/2.0, 2));
        }
    }
    
    int num_steps = 500; // Total time T = 10.0 allows the packet to return twice under v=1
    std::vector<int> ipiv(total_nodes);
    std::ofstream out_file("./data/results_" + mode_name + ".txt");
    
    // Time-stepping loop
    for (int step = 1; step <= num_steps; ++step) {
        // Calculate Right-Hand Side vector: B = M_RHS * c_vec
        std::vector<double> B(total_nodes, 0.0);
        for (int i = 0; i < total_nodes; ++i) {
            for (int j = 0; j < total_nodes; ++j) {
                B[i] += M_RHS[i * total_nodes + j] * c_vec[j];
            }
        }
        
        // LAPACK overwrites the LHS matrix, so a fresh copy is required for each step
        std::vector<double> A_copy = M_LHS;
        
        // Solve the linear system using LAPACKE
        lapack_int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, total_nodes, 1, A_copy.data(), total_nodes, ipiv.data(), B.data(), 1);
        if (info != 0) {
            std::cerr << "LAPACK error: info = " << info << "\n";
            return;
        }
        
        // Update the solution vector
        c_vec = B; 
        
        // Compute statistics (min, max, and center of mass) to verify lab deliverables
        double c_min = c_vec[0], c_max = c_vec[0];
        double center_x = 0, center_y = 0, mass = 0;
        for (int ej = 0; ej < N; ++ej) {
            for (int ei = 0; ei < N; ++ei) {
                double val = c_vec[ej * N + ei];
                if (val < c_min) c_min = val;
                if (val > c_max) c_max = val;
                
                center_x += (ei * a) * val;
                center_y += (ej * a) * val;
                mass += val;
            }
        }
        if (mass > 0) {
            center_x /= mass;
            center_y /= mass;
        }
        
        // Log output periodically and save to file
        if (step % 50 == 0 || step == num_steps) {
            std::cout << "Step: " << step << ", t: " << step * dt 
                      << " | max_u: " << c_max << " min_u: " << c_min 
                      << " | center: (" << center_x << ", " << center_y << ")\n";
        }
        if (step % 5 == 0 || step == num_steps) {
            out_file << step * dt << " " << c_min << " " << c_max << " " << center_x << " " << center_y << "\n";
        }
    }
    out_file.close();
    std::cout << "Finished " << mode_name << "\n\n";
}