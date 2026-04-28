#include "data_processing.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>

void write_to_file(const std::string& filename, const std::string& content) {
    std::filesystem::path dir = std::filesystem::path(filename).parent_path();
    if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
    }
    std::ofstream file(filename, std::ios::app);
    file << content;
    file.close();
}

void save_grid_to_file(const std::string& filename, const std::vector<std::vector<double>>& grid, int N) {
    std::filesystem::path dir = std::filesystem::path(filename).parent_path();
    if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
    }
    std::ofstream file(filename);
    for (int row = 0; row < N; ++row) {
        for (int col = 0; col < N; ++col) {
            file << std::fixed << std::setprecision(6) << grid[row][col] << " ";
        }
        file << "\n";
    }
    file.close();
}

void cleanup_directory(const std::string& dirname) {
    if (std::filesystem::exists(dirname)) {
        for (const auto& entry : std::filesystem::directory_iterator(dirname)) {
            if (std::filesystem::is_regular_file(entry)) {
                std::filesystem::remove(entry);
            }
        }
    }
}