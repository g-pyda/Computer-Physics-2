#include "engine.h"
#include <cmath>
#include <iostream>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -O3

int main() {
    const double L = 5.0;

    // std::cout << "Starting Convergence Test (omega = pi/L)..." << std::endl;
    // run_convergence_test("convergence_omega_pi_L", M_PI / L);

    // std::cout << "Starting simulation for omega = pi/L..." << std::endl;
    // run_experiment("omega_pi_L", M_PI / L);

    // std::cout << "\nStarting simulation for omega = 2*pi/L..." << std::endl;
    // run_experiment("omega_2pi_L", 2.0 * M_PI / L);

    // std::cout << "\nStarting simulation for omega = pi/(2*L)..." << std::endl;
    // run_experiment("omega_half_pi_L", M_PI / (2.0 * L));

    std::cout << "\nStarting simulation for omega = sqrt(2.0)*pi/L..." << std::endl;
    run_experiment("omega_sqrt2_pi_L", std::sqrt(2.0) * M_PI / L);

    std::cout << "\nAll simulations completed. Data saved in 'data' directory." << std::endl;

    return 0;
}