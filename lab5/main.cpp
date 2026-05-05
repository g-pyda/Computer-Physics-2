#include <iostream>
#include <vector>
#include <string.h>
#include "engine.h"
#include <time.h>
#include <iomanip>
#include <complex>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -lblas
// or
// g++ -I ~/lapacke/lapacke/include/ -L ~/lapacke/lapacke/ main.cpp engine.cpp -llapacke -llapack -lblas -lm

int main() {
    std::cout << "Initializing the experiments..." << std::endl;
    std::vector<std::vector<std::vector<double>>> excitations;

    // INITIAL CONDITION GENERATON

    double alpha_crit = 0.03125;
    double alpha_optimal = 0.6*alpha_crit;
    int max_excitation = 3;
    for (int k = 0; k < max_excitation; ++k) {
        std::cout << "Generating excitation for K = " << k << "/" << max_excitation-1 << std::endl;
        generate_excitation(alpha_optimal, k, "data/excitation", excitations, true);
    }

    // INITIAL CONDITION WITH 3 EIGENSTATES

    std::cout << "Running the experiment with proper setup" << std::endl;
    std::vector<std::vector<std::complex<double>>> proper_init;
    get_initial_condition(excitations, proper_init, true);
    run_simulation(proper_init, "./data/proper_init/");


    std::cout << "Experiments completed successfully." << std::endl;

    return 0;
}