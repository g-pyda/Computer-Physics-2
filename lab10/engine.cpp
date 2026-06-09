#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <filesystem>
#include <iomanip>



// ============ CONSTANTS ============
const double L = 1.0;
const int N = 51;
const double delta_x = 1;
const double D_center = 1.0;
const double D_edges_1 = 0.0001;
const double D_edges_2 = 0.005;
const double s = 0.01;
const double delta_t = 0.125;
const double max_t = 250.0;

// =========== FLUX FUNCTION =============
double D_comp(std::vector<std::vector<double>>& D, int i, int j, int direction) {
    switch (direction) {
        case 0: // flux to the left
            return 0.5*(D[i][j] + D[i-1][j]);
        case 1: // flux to the bottom
            return 0.5*(D[i][j] + D[i][j-1]);
        case 2: // flux from the right
            return 0.5*(D[i+1][j] + D[i][j]);
        case 3: // flux from the top
            return 0.5*(D[i][j+1] + D[i][j]);
    }
}

double F(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& D, int i, int j, int direction) {
    double D_calc = -1.0*D_comp(D, i, j, direction);
    switch (direction) {
        case 0: // flux to the left
            return D_calc*(u[i][j] - u[i-1][j])/delta_x;
        case 1: // flux to the bottom
            return D_calc*(u[i][j] - u[i][j-1])/delta_x;
        case 2: // flux from the right
            return D_calc*(u[i+1][j] - u[i][j])/delta_x;
        case 3: // flux from the top
            return D_calc*(u[i][j+1] - u[i][j])/delta_x;
    }
}

// ================ SOURCE ================
double S(int i, int j) {
    if (i == (N) && j == (N)) return s;
    return 0;
}

// ======= SNAPSHOTING =========
void save_snapshot(const std::vector<std::vector<double>>& u, double t) {
    std::ofstream snapshot("./data/"+std::to_string(t)+".txt");
    for (int i = 0; i < 2*N+1; ++i) {
        for (int j = 0; j < 2*N+1; ++j) {
            snapshot << u[i][j] << " ";
        }
        snapshot << std::endl;
    }
}

// ========== SIMULATION CODE =============
void run_simulation() {
    // grid init
    std::vector<std::vector<double>> u (2*N+1, std::vector<double>(2*N+1, 0.0));
    std::vector<std::vector<double>> u_new (2*N+1, std::vector<double>(2*N+1, 0.0));
    std::vector<std::vector<double>> D (2*N+1, std::vector<double>(2*N+1, D_center));
    for (int i = 0; i < 2*N+1; ++i) {
        for (int j = 0; j < 5; ++j) {
            D[i][j] = D_edges_1;
            D[i][2*N-j] = D_edges_1;
            D[j][i] = D_edges_1;
            D[2*N-j][i] = D_edges_1;
        }
        if (i < 5 || i > 2*N-5) continue;
        for (int j = 5; j < 10; ++j) {
            D[i][j] = D_edges_2;
            D[i][2*N-j] = D_edges_2;
            D[j][i] = D_edges_2;
            D[2*N-j][i] = D_edges_2;
        }
    }
    // time evolution
    std::filesystem::create_directories("./data/");
    double t = 0.0;
    save_snapshot(u, t);
    while (t < max_t) {
        t += delta_t;
        std::cout << "\rTimestep " << std::setprecision(3) << t << "/" << std::setprecision(3) << max_t << std::flush;
        for (int i = 1; i < 2*N; ++i) {
            for (int j = 1; j < 2*N; ++j) {
                u_new[i][j] = u[i][j] + delta_t * (
                    (F(u, D, i, j, 0) - F(u, D, i, j, 2)) / delta_x
                    + (F(u, D, i, j, 1) - F(u, D, i, j, 3)) / delta_x
                    + S(i, j)
                );
            }
        }
        u = u_new;
        // saving the snapshot
        save_snapshot(u, t);
    }
    std::cout << "\nSimulation completed!" << std::endl;
    
    // calculating the integral
    double integral = 0.0;
    for (int i = 0; i < 2*N+1; ++i) {
        for (int j = 0; j < 2*N+1; ++j) {
            integral += u[i][j] * delta_x * delta_x;
        }
    }
    std::cout << "Integral: " << integral << " | st: " << s*max_t << std::endl;
}