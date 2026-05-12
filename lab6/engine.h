#ifndef LAB2_ENGINE_H
#define LAB2_ENGINE_H

#include <vector>
#include <string>

// ========== STRUCTURES ========== //

struct Node {
    double x, y;
};

void run_simulation (
    std::string path,
    bool half_omega
);


#endif