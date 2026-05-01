#ifndef LAB2_ENGINE_H
#define LAB2_ENGINE_H

#include <vector>
#include <string>


void run_simulation(
    double alpha,
    int excited_states,
    std::string path,
    std::vector<std::vector<std::vector<double>>>& excitations,
    bool save_excitation
);

#endif