#ifndef LAB2_ENGINE_H
#define LAB2_ENGINE_H

#include <vector>
#include <string>
#include <complex>


void generate_excitation(
    double alpha,
    int excited_states,
    std::string path,
    std::vector<std::vector<std::vector<double>>>& excitations,
    bool save_excitation
);

void run_simulation(
    std::vector<std::vector<std::complex<double>>>& init,
    std::string path
);

void get_initial_condition(
    std::vector<std::vector<std::vector<double>>>& excitations,
    std::vector<std::vector<std::complex<double>>>& init,
    bool proper
);


#endif