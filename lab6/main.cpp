#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "engine.h"
#include <time.h>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -lblas
// or
// g++ -I ~/lapacke/lapacke/include/ -L ~/lapacke/lapacke/ main.cpp engine.cpp -llapacke -llapack -lblas -lm

int main() {
    std::cout << "Starting experiments..." << std::endl;
    
    // experiment for full omega
    std::cout << "Experiment for first two degenerate states with full omega..." << std::endl;
    run_simulation("full_omega", false, 1.0/5.0);

    // experiment for half of omega
    std::cout << "Experiment for first two degenerate states with half omega..." << std::endl;
    run_simulation("half_omega", true, 1.0/5.0);

    // experiment over different eF values
    std::cout << "Experiment for different eF values..." << std::endl;
    std::vector<double> factors = {0.1, 0.5, 1.0, 2.0, 5.0};
    for(auto& factor : factors) {
        double eF = factor / 5.0;
        std::ostringstream stream;
        stream << std::fixed << std::setprecision(2) << eF;
        std::string s_value = stream.str();
        std::string directory = "ef/" + s_value;
        std::cout << "-> eF = " << eF << std::endl;
        run_simulation(directory, false, eF);
    }

    std::cout << "Experiments completed successfully." << std::endl;

    return 0;
}