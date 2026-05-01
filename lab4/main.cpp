#include <iostream>
#include <vector>
#include <string.h>
#include "engine.h"
#include <time.h>
#include <iomanip>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -lblas
// or
// g++ -I ~/lapacke/lapacke/include/ -L ~/lapacke/lapacke/ main.cpp engine.cpp -llapacke -llapack -lblas -lm

int main() {
    std::cout << "Starting experiments..." << std::endl;
    std::vector<std::vector<std::vector<double>>> excitations;

    // GROUND STATE GENERATION ANALYSIS

    double alpha_crit = 0.03125;
    std::vector<double> factors = {
        0.99, 0.9,                              // near critical values
        0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1  // stable values
    };
    for (double f : factors) {
        double alpha = alpha_crit * f;
        std::cout << "Running experiment with alpha = " << std::fixed << std::setprecision(6) << alpha << std::endl;
        run_simulation(alpha, 0, "data/grnd_alpha_"+std::to_string(alpha), excitations, false);
    }

    // EXCITED STATES GENERATION ANALYSIS
    double alpha_optimal = 0.6*alpha_crit;
    int max_excitation = 5;
    for (int k = 0; k < max_excitation; ++k) {
        std::cout << "Running excitation experiment for K = " << k << "/" << max_excitation-1 << std::endl;
        run_simulation(alpha_optimal, k, "data/excitation", excitations, true);
    }

    std::cout << "Experiments completed successfully." << std::endl;

    return 0;
}