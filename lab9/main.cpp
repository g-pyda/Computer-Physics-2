#include "engine.h"
#include <cmath>
#include <iostream>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -lblas -O3
// or
// g++ -I ~/lapacke/lapacke/include/ -L ~/lapacke/lapacke/ main.cpp engine.cpp -llapacke -llapack -lblas -lm

int main() {
    int N = 40;
    
    // Deliverable 1: Pure advection (D=0, vx=vy=1)
    run_simulation(0.0, 1.0, 1.0, N, "pure_advection");
    
    // Deliverable 2: Pure diffusion (D=0.1, vx=vy=0)
    run_simulation(0.1, 0.0, 0.0, N, "pure_diffusion");
    
    // Deliverable 3: Advection-diffusion (D=0.1, vx=vy=1)
    run_simulation(0.1, 1.0, 1.0, N, "advection_diffusion");
    
    return 0;
}