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

// ========== GROUND STATE GENERATION ========= //

void get_ground_state(
    std::vector<std::vector<double>>& wavefunction,
    std::string path,
    double alpha
) {
    double prev_en = 0.0;
    double diff = 0.0;
    std::ofstream f_energy (path + "/grnd_energy.out");
    for (int e = 0; e < MAX_ITER; ++e) {
        if ((e+1)%10 == 0)
            std::cout << "Iteration " << e+1 << "\tenergy diff: " << diff << "\r";
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
    std::cout << "\n";
    f_energy.close();
}

// ========== EXCITED STATES GENERATION ========== //

double c_k(
    std::vector<std::vector<std::vector<double>>>& excitations,
    std::vector<std::vector<double>>& wavefunction,
    int k
) {
    double c = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            c += excitations[k][i][j]*wavefunction[i][j]*std::pow(DELTA, 2);
        }
    }
    return c;
}

void get_excited_state(
    std::vector<std::vector<std::vector<double>>>& excitations,
    std::vector<std::vector<double>>& wavefunction,
    std::string path,
    double alpha,
    int excitation_level
) {
    double prev_en = 0.0;
    double diff = 0.0;
    std::ofstream f_energy (path + "/" + std::to_string(excitation_level) + "_energy.out");
    for (int e = 0; e < MAX_ITER; ++e) {
        if ((e+1)%10 == 0)
            std::cout << "Iteration " << e+1 << "\tenergy diff: " << diff << "\r";

        // iterative application of Hamiltonian on the wavefunction
        std::vector<std::vector<double>> hamil (N, std::vector<double>(N, 0.0));
        for (int i = 1; i < N-1; ++i) 
            for (int j = 1; j < N-1; ++j) 
                hamil[i][j] =  operate_hamiltonian(wavefunction, i, j);
                
        for (int i = 1; i < N-1; ++i) 
            for (int j = 1; j < N-1; ++j) 
                wavefunction[i][j] = wavefunction[i][j] - alpha*hamil[i][j];
        
        // orthogonalization
        std::vector<double> c(excitation_level, 0.0);
        for (int k = 0; k < excitation_level; ++k) 
            c[k] = c_k(excitations, wavefunction, k);
        
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < excitation_level; ++k) {
                    wavefunction[i][j] = wavefunction[i][j] - c[k]*excitations[k][i][j];
                }
            }
        }

        // wave normalization
        double I = 0.0;
        for (int i = 0; i < N; ++i) 
            for (int j = 0; j < N; ++j) 
                I += std::pow(wavefunction[i][j]*DELTA, 2);

        I = std::sqrt(I);
        for (int i = 0; i < N; ++i) 
            for (int j = 0; j < N; ++j) 
                wavefunction[i][j] /= I;
            
        // saving the energy to a file
        double current_en = mean_energy(wavefunction);

        int nx = excitation_level, ny = 0;
        if (excitation_level == 3) {
            nx = 0;
            ny = 1;
        } else if (excitation_level == 4) nx -= 1;

        f_energy << std::to_string(e+1) << "\t" << std::setprecision(6) << current_en << " " << std::setprecision(6) << get_exact_energy(nx, ny) << "\n";
        
        // convergence check
        diff = std::abs(current_en - prev_en);
        if (e > 0) 
            if (diff < MAX_TOL)
                break;
        
        prev_en = current_en;
    }
    std::cout << "\n";
    f_energy.close();
}

void save_wavefunction(std::vector<std::vector<double>>& wavefunction, std::string path) {
    std::ofstream f_wf (path + "_wf.out");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            f_wf << wavefunction[i][j] << " ";
        }
        f_wf << "\n";
    }
    f_wf.close();
}

// ========== MAIN SIMULATION ========== //

void run_simulation(double alpha, int K, std::string path, std::vector<std::vector<std::vector<double>>>& excitations, bool save_excitation) {
    // initializing file directory
        std::filesystem::create_directories(path);

        // initializing the mesh with random values (except BC)
        std::vector<std::vector<double>> wavefunction(N, std::vector<double>(N, 0.0));
        for (int i = 1; i < N-1; ++i) 
            for (int j = 1; j < N-1; ++j)
                wavefunction[i][j] = double(rand())/(0.5*RAND_MAX) - 1.0;

    if (K == 0) {
        // generation of ground state 
        get_ground_state(wavefunction, path, alpha);
    } else {
        // generating excitation state based on previous ones
        get_excited_state(excitations, wavefunction, path, alpha, K);
        
    }
    if (save_excitation) {
        excitations.push_back(wavefunction);
        save_wavefunction(wavefunction, path + "/" + std::to_string(K));
    }
}