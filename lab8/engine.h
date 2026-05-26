#ifndef LAB2_ENGINE_H
#define LAB2_ENGINE_H

#include <vector>
#include <string>


struct Node {
    double x, y;
};

void run_simulation (
    std::string path,
    double omega,
    double damp_const
);


#endif