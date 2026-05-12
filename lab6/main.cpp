#include <iostream>
#include <vector>
#include "engine.h"
#include <time.h>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -lblas
// or
// g++ -I ~/lapacke/lapacke/include/ -L ~/lapacke/lapacke/ main.cpp engine.cpp -llapacke -llapack -lblas -lm

int main() {
    std::cout << "Starting experiments..." << std::endl;

    std::cout << "Experiment for first two degenerate states with full omega" << std::endl;
    run_simulation("full_omega", false);

    std::cout << "Experiment for first two degenerate states with half omega" << std::endl;
    run_simulation("half_omega", true);

    std::cout << "Experiments completed successfully." << std::endl;

    return 0;
}