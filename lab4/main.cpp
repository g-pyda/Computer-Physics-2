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

    double alpha_crit = 0.03125;
    std::vector<double> alphas = {
        0.03125, 0.031, 0.030, // values near critical
        0.028125
    };
    for (double alpha : alphas) {
        std::cout << "Running experiment with alpha = " << std::fixed << std::setprecision(6) << alpha << std::endl;
        run_simulation(alpha, 0, "data/grnd_alpha_"+std::to_string(alpha));
    }

    std::cout << "Experiments completed successfully." << std::endl;

    return 0;
}