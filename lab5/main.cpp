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
    std::vector<std::vector<std::complex<double>>> init;

    // INITIAL CONDITION GENERATON

    double alpha_crit = 0.0078125;
    double alpha_optimal = 0.6*alpha_crit;
    int max_excitation = 3;
    for (int k = 0; k < max_excitation; ++k) {
        std::cout << "Generating excitation for K = " << k << "/" << max_excitation-1 << std::endl;
        generate_excitation(alpha_optimal, k, "data/excitation", excitations, true);
    }

    // MULTIPLE EIGENSTATES

    std::cout << "Running the experiment with multi-eigenvalue setup" << std::endl;
    
    get_initial_condition(excitations, init, 0);
    run_simulation(init, "./data/mult/");

    // MULTIPLE EIGENSTATES - REVERSE ROTATION

    std::cout << "Running the experiment with reverse rotation" << std::endl;
    get_initial_condition(excitations, init, 1);
    run_simulation(init, "./data/mult_rev/");

    // OSCILLATION - X DIRECTION

    std::cout << "Running the experiment with x-direcetd oscillation" << std::endl;
    get_initial_condition(excitations, init, 2);
    run_simulation(init, "./data/osc_x/");

    // OSCILLATION - Y DIRECTION

    std::cout << "Running the experiment with y-directed oscilation" << std::endl;
    get_initial_condition(excitations, init, 3);
    run_simulation(init, "./data/osc_y/");

    // ONE EIGENSTATE

    std::cout << "Running the experiment with single eigenstate" << std::endl;
    get_initial_condition(excitations, init, 4);
    run_simulation(init, "./data/sing/");


    std::cout << "Experiments completed successfully." << std::endl;

    return 0;
}