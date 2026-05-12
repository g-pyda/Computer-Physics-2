#ifndef LAB2_ENGINE_H
#define LAB2_ENGINE_H

#include <vector>
#include <string>

// ========== STRUCTURES ========== //

struct Node {
    double x, y;
};

// ========== CORE FUNCTIONS ========== //

double get_exact_potential(double x, double y);

void run_simulation(int N, double L, bool biparabolic);

// ========== SHAPE FUNCTIONS ========== //

double g_func(int i, double xi1, double xi2);

double h_func(int i, double xi1, double xi2);

#endif