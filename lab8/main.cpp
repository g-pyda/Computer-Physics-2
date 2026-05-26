#include "engine.h"
#include <cmath>
#include <iostream>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -lblas -O3
// or
// g++ -I ~/lapacke/lapacke/include/ -L ~/lapacke/lapacke/ main.cpp engine.cpp -llapacke -llapack -lblas -lm

int main() {
    //experiment over different omega values
    for (double factor = 0.5; factor <= 5.0; factor += 0.05) {
        double omega = factor * M_PI;
        std::cout << "\nStarting simulation for omega = " << factor << "*pi..." << std::endl;
        run_simulation("omega/omega_" + std::to_string(factor) + "_pi", omega, 1.0);
    }

    //experiment over different damping values
    for (double d = 0.0; d <= 2.0; d += 0.1) {
        std::cout << "\nStarting simulation for damping = " << d << "..." << std::endl;
        run_simulation("damping/damping_" + std::to_string(d), M_PI, d);        
    }

    // hybrid experiment over omega and damping
    for (double factor = 0.5; factor <= 5.0; factor += 0.1) {
        double omega = factor * M_PI;
        for (double d = 0.0; d <= 2.0; d += 0.2) {
            std::cout << "\nStarting simulation for omega = " << factor << "*pi and damping = " << d << "..." << std::endl;
            run_simulation("hybrid/omega_" + std::to_string(factor) + "_pi_damping_" + std::to_string(d), omega, d);
        }   
    }

    std::cout << "\nAll simulations completed. Data saved in 'data' directory." << std::endl;

    return 0;
}