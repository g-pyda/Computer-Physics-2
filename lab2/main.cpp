#include <iostream>
#include <vector>
#include "engine.h"
#include <time.h>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -lblas

int main() {
    double L = 5.0;
    std::cout << "Starting experiments..." << std::endl;

    std::vector<int> bilinear_N_values = {4, 8, 16, 32};
    for (int N : bilinear_N_values) {
        std::cout << "Running bilinear experiment with N = " << N << std::endl;
        run_simulation(N, L, false);
    
        std::cout << "Running biparabolic experiment with N = " << N << std::endl;
        run_simulation(N, L, true);
    }

    std::cout << "Experiments completed successfully." << std::endl;

    return 0;
}