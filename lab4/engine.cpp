#include "engine.h"
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <fstream>
#include <filesystem>

// ============= CONSTANTS ================= //

const double M = 1.0;
const double H_c = 1.0;
const double Wx = 1.0;
const double Wy = 2.0;

const int N = 16;
const double DELTA = 4.0/N;
const double L = DELTA*N;

const int MAX_ITER = 1e3;
const double MAX_TOL = 1e-5;

// =========== PHYSICS OPERATIONS =============== //

double get_exact_energy(double nx, double ny) {
    return H_c*(Wx*(nx + 0.5) + Wy*(ny + 0.5));
}

double potential(double x, double y) {
    return 0.5 * M * (std::pow(Wx*x, 2) + std::pow(Wy*y, 2));
}

double operate_hamiltonian(
    std::vector<std::vector<double>>& wavef,
    int i,
    int j
) {
    double result = 0.0;
    result += (wavef[i+1][j] + wavef[i-1][j] - 2*wavef[i][j]) / std::pow(DELTA, 2);
    result += (wavef[i][j+1] + wavef[i][j-1] - 2*wavef[i][j]) / std::pow(DELTA, 2);
    result *= -0.5 * std::pow(H_c, 2) * M;

    return result + potential(DELTA*i - 0.5*L, DELTA*j - 0.5*L)*wavef[i][j];
}

double mean_energy(
    std::vector<std::vector<double>>& wavef
) {
    double energy = 0.0;
    for (int i = 1; i < N-1; ++i)
        for (int j = 1; j < N-1; ++j)
            energy += wavef[i][j]*operate_hamiltonian(wavef, i, j)*std::pow(DELTA,2);

    return energy;
}

// ========== LINEAR SYSTEM SOLVER (LAPACKE) ========== //

// void solve_and_save(
//     std::vector<std::vector<double>>& H,
//     std::vector<std::vector<double>>& O,
//     std::vector<double>& E,
//     const std::string& path,
//     int total_nodes
// ) {
//     int n = total_nodes;
//     int nrhs = 1;
//     int lda = n;
//     int ldb = n;

//     // creating flat buffers for LAPACK
//     std::vector<double> H_flat(n * n);
//     std::vector<double> O_flat(n * n);
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             H_flat[i * n + j] = H[i][j];
//             O_flat[i * n + j] = O[i][j];
//         }
//     }

//     // solving generalized symmetric-definite eigenproblem (eigenvalues in E, eigenvectors in H_flat)
//     lapack_int info = LAPACKE_dsygv(
//         LAPACK_ROW_MAJOR, 
//         1,
//         'V',
//         'U',
//         n,
//         H_flat.data(),
//         lda,
//         O_flat.data(),
//         ldb,
//         E.data()
//     );

//     if (info > 0) {
//         std::cerr << "error: matrix is singular at U(" << info << "," << info << ")" << std::endl;
//         return;
//     } else if (info < 0) {
//         std::cerr << "error: illegal argument at position " << -info << " in lapack call" << std::endl;
//         return;
//     }

//     // saving eigenvalues to eigenvalues.txt
//     std::ofstream f_sol(path + "/eigenvalues.txt");
//     if (!f_sol.is_open()) {
//         std::cerr << "error: could not open solution file for writing" << std::endl;
//         return;
//     }

//     for (int i = 0; i < n; i++) {
//         f_sol << i + 1 << " " << std::fixed << std::setprecision(10) << E[i] << "\n";
//     }
//     f_sol.close();

//     // saving eigenvectors to eigenvectors.txt
//     std::ofstream f_vec(path + "/eigenvectors.txt");
//     if (!f_vec.is_open()) {
//         std::cerr << "error: could not open eigenvectors file for writing" << std::endl;
//         return;
//     }

//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j)
//             f_vec << std::fixed << std::setprecision(10) << H_flat[i*n + j] << " ";
//         f_vec << "\n";
//     }
//     f_vec.close();
// }

// ========== MAIN SIMULATION FROM GROUND STATE ========== //

void run_simulation(double alpha, int excited_states, std::string path) {
    // initializing file directory
    std::filesystem::create_directories(path);

    // initializing the mesh with random values (except BC)
    std::vector<std::vector<double>> wavefunction(N, std::vector<double>(N, 0.0));
    for (int i = 1; i < N-1; ++i) 
        for (int j = 1; j < N-1; ++j)
            wavefunction[i][j] = double(rand())/(0.5*RAND_MAX) - 1.0;
    
    double prev_en = 0.0;
    double diff = 0.0;
    std::ofstream f_energy (path + "/grnd_energy.out");
    for (int e = 0; e < MAX_ITER; ++e) {
        if ((e+1)%10 == 0)
            std::cout << "Iteration " << e+1 << "\tenergy diff: " << diff << "\n";
        // iterative application of Hamiltonian on the wavefunction
        std::vector<std::vector<double>> hamil (N, std::vector<double>(N, 0.0));
        for (int i = 1; i < N-1; ++i) 
            for (int j = 1; j < N-1; ++j) 
                hamil[i][j] =  operate_hamiltonian(wavefunction, i, j);
                
        for (int i = 1; i < N-1; ++i) 
            for (int j = 1; j < N-1; ++j) 
                wavefunction[i][j] = wavefunction[i][j] - alpha*hamil[i][j];

        // wave normalization
        double I = 0.0;
        for (int i = 1; i < N-1; ++i) 
            for (int j = 1; j < N-1; ++j) 
                I += std::pow(wavefunction[i][j]*DELTA, 2);

        I = std::sqrt(I);
        for (int i = 1; i < N-1; ++i) 
            for (int j = 1; j < N-1; ++j) 
                wavefunction[i][j] /= I;
            
        // saving the energy to a file
        double current_en = mean_energy(wavefunction);
        f_energy << std::to_string(e+1) << "\t" << std::setprecision(6) << current_en << " " << std::setprecision(6) << get_exact_energy(0, 0) << "\n";
        
        // convergence check
        diff = std::abs(current_en - prev_en);
        if (e > 0) 
            if (diff < MAX_TOL)
                break;
        
        prev_en = current_en;
    }
    f_energy.close();
}