#ifndef DATA_PROCESSING_H
#define DATA_PROCESSING_H

#include <string>
#include <vector>

void write_to_file(const std::string& filename, const std::string& content);
void save_grid_to_file(const std::string& filename, const std::vector<std::vector<double>>& grid, int N);
void cleanup_directory(const std::string& dirname);
void visualize_grids(const std::string& study, const std::vector<std::string>& cases, const std::string& output_filename);
void plot_energy_functional(int max_iterations, const std::string& study, const std::vector<std::string>& cases, const std::string& output_filename);

#endif // DATA_PROCESSING_H