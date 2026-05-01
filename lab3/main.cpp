#include <iostream>
#include <vector>
#include "engine.h"
#include <time.h>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -lblas
// or
// g++ -I ~/lapacke/lapacke/include/ -L ~/lapacke/lapacke/ main.cpp engine.cpp -llapacke -llapack -lblas -lm

int main() {
    double L = 5.0;
    std::cout << "Starting experiments..." << std::endl;

    std::vector<int> N_values = {4, 8, 16, 32};
    for (int N : N_values) {
        std::cout << "Running bilinear experiment with N = " << N << std::endl;
        run_simulation(N, L, false);
    
        std::cout << "Running biparabolic experiment with N = " << N << std::endl;
        run_simulation(N, L, true);
    }

    std::cout << "Experiments completed successfully." << std::endl;

    return 0;
}