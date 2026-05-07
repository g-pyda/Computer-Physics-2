#include "engine.h"
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <complex>

using namespace std::complex_literals;

// ============= CONSTANTS ================= //

const double M = 1.0;
const double H_c = 1.0;
const double Wx = 1.0;
const double Wy = 1.001;

const int N = 16;
const double DELTA = 4.0/N;
const double L = DELTA*N;

const double DELTA_T = 1e-3;
const double MAX_T = 4*M_PI;
const int MAX_ITER_CN = 5;
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

std::complex<double> operate_hamiltonian(
    std::vector<std::vector<std::complex<double>>>& wavef,
    int i,
    int j
) {
    std::complex<double> result = 0.0;

    double inv_dx2 = 1.0 / (DELTA * DELTA);

    result += (wavef[i+1][j] + wavef[i-1][j] - 2.0 * wavef[i][j]) * inv_dx2;
    result += (wavef[i][j+1] + wavef[i][j-1] - 2.0 * wavef[i][j]) * inv_dx2;

    result *= -0.5 * H_c * H_c / M;

    return result + potential(DELTA*i - 0.5*L, DELTA*j - 0.5*L) * wavef[i][j];
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


double mean_energy(
    std::vector<std::vector<std::complex<double>>>& wavef
) {
    std::complex<double> energy = 0.0;

    double dA = DELTA * DELTA;

    for (int i = 1; i < N-1; ++i)
        for (int j = 1; j < N-1; ++j)
            energy += std::conj(wavef[i][j]) 
                    * operate_hamiltonian(wavef, i, j) 
                    * dA;

    return energy.real();
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

// ========== WAVEFUNCTION GENERATION ========== //

void generate_excitation(
    double alpha,
    int K,
    std::string path,
    std::vector<std::vector<std::vector<double>>>& excitations,
    bool save_excitation) 
    {
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
    }
}

void get_initial_condition(
    std::vector<std::vector<std::vector<double>>>& excitations,
    std::vector<std::vector<std::complex<double>>>& init,
    bool proper
) {
    init.assign(N, std::vector<std::complex<double>>(N, 0.0 + 0i));
    double norm = 1.0 / std::sqrt(3);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (proper == true)
                init[i][j] = std::complex<double>(norm*(excitations[0][i][j] + excitations[1][i][j]), norm*excitations[2][i][j]);
            else init[i][j] = excitations[2][i][j] + 0i;
        }
    }
}

void run_simulation(
    std::vector<std::vector<std::complex<double>>>& init,
    std::string path
) {
    // initializing file directory
    std::filesystem::create_directories(path);

    // mid-step holders initialization
    std::vector<std::vector<std::complex<double>>> wavef_t (
        N, std::vector<std::complex<double>> (N, 0.0)
    );
    std::vector<std::vector<std::complex<double>>> wavef_t_dt (
        N, std::vector<std::complex<double>> (N, 0.0)
    );
    std::vector<std::vector<std::complex<double>>> wavef_new (
        N, std::vector<std::complex<double>> (N, 0.0)
    );

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            wavef_t[i][j] = init[i][j];
        }
    }

    // file setup
    std::ofstream f_deliverables (path + "deliverables.txt");

    // simulation for all time steps
    std::complex<double> wavef_norm = DELTA_T / std::complex<double>(0.0, (2*H_c));
    for (double t = 0; t < MAX_T; t += DELTA_T) {
        // initialization of starting point
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                wavef_t_dt[i][j] = wavef_t[i][j];
            }
        }

        // Crank-Nicolson scheme
        for (int iter = 0; iter < MAX_ITER_CN; ++iter) {
            for (int i = 1; i < N-1; ++i) {
                for (int j = 1; j < N-1; ++j) {
                    wavef_new[i][j] = wavef_t[i][j] + wavef_norm*(
                        operate_hamiltonian(wavef_t, i, j) + operate_hamiltonian(wavef_t_dt, i, j)
                    );
                }
            }

            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    wavef_t_dt[i][j] = wavef_new[i][j];
                }
            }
        }

        // normalization
        std::complex<double> I = 0.0;
        for (int i = 1; i < N-1; ++i) 
            for (int j = 1; j < N-1; ++j) 
                I += std::norm(wavef_t_dt[i][j])*std::pow(DELTA, 2);

        I = std::sqrt(I);
        for (int i = 1; i < N-1; ++i) 
            for (int j = 1; j < N-1; ++j) 
                wavef_t_dt[i][j] /= I;

        // saving the deliverables
        double avg_x = 0.0, avg_y = 0.0;
        double dA = std::pow(DELTA, 2); // Area element dx*dy

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double x_coord = DELTA * i - 0.5 * L;
                double y_coord = DELTA * j - 0.5 * L;
        
                double prob_density = std::norm(wavef_t_dt[i][j]);
        
                avg_x += prob_density * x_coord * dA;
                avg_y += prob_density * y_coord * dA;
            }
        }

        f_deliverables << t << "\t"
            << std::setprecision(6) << std::norm(I) << "\t"
            << std::setprecision(6) << mean_energy(wavef_t_dt) << "\t"
            << std::setprecision(6) << avg_x << "\t"
            << std::setprecision(6) << avg_y << "\t"
            << std::endl;

        // shift of time
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                wavef_t[i][j] = wavef_t_dt[i][j];
            }
        }
    }
    f_deliverables.close();
}
