#ifndef ENGINE_H
#define ENGINE_H

#include <vector>
#include <string>
#include <utility>

std::pair<std::vector<std::vector<double>>, double> create_grid(int N);
double poisson_equation(double x, double y);
double residuum(const std::vector<std::vector<double>>& u, int i, int j, int N, double delta_x);
double energy_functional(const std::vector<std::vector<double>>& u, int N, double delta_x);
void gauss_seidel(const std::string& dirname, const char study, std::vector<std::vector<double>>& u, int N, double delta_x, double weight = 1.0, int max_iterations = 1000, double tolerance = 1e-6);

#endif // ENGINE_H