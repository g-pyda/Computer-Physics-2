#ifndef ENGINE_H
#define ENGINE_H

#include <string>

void run_experiment(const std::string& expName, double omega);
void run_convergence_test(const std::string& expName, double omega);

#endif