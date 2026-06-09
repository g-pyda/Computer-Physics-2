#include "engine.h"
#include <cmath>
#include <iostream>

// g++ main.cpp engine.cpp -o program -llapacke -llapack -lblas -O3
// or
// g++ -I ~/lapacke/lapacke/include/ -L ~/lapacke/lapacke/ main.cpp engine.cpp -llapacke -llapack -lblas -lm

int main() {
    run_simulation();
    
    return 0;
}