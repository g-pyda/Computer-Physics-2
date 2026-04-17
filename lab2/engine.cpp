#include <stdio.h>
#include <vector>

// ========== BILINEAR SHAPE FUNCTION ========== //
double _f1(double x) {
    return (1 - x) / 2.0;
}

double _f2(double x) {
    return (1 + x) / 2.0;
}

double _g(int i, double x1, double x2) {
    double result = 1.0;

    if (i % 2 == 1)
        result *= _f1(x1) * _f2(x2);
    else
        result *= _f1(x2) * _f2(x1);

    return result;
}

// ========== BIPARABOLIC SHAPE FUNCTION ========== //
double _q1(double x) {
    return x * (x-1) / 2.0;
}

double _q2(double x) {
    return (1-x) * (1+x);
}

double _q3(double x) {
    return x * (x+1) / 2.0;
}

double _h(int i, double x1, double x2) {
    double result = 1.0;
    // first multiplication
    switch (i) {
        case 1:
        case 3:
        case 7:
            result *= _q1(x1);
            break;
        case 2:
        case 4:
        case 6:
            result *= _q3(x1);
            break;
        case 5:
        case 8:
        case 9:
            result *= _q2(x1);
            break;
    }

    // second multiplication
    switch (i) {
        case 1:
        case 2:
        case 5:
            result *= _q1(x2);
            break;
        case 3:
        case 4:
        case 8:
            result *= _q3(x2);
            break;
        case 6:
        case 7:
        case 9:
            result *= _q2(x2);
            break;
    }

    return result;
}

// ========== EXPORTABLE ============ //

void process_node(
    std::vector<std::vector<double>>& x,
    std::vector<std::vector<double>>& y,
    int node_id,
    int N,
    double L,
    bool biparabolic
) {
    double jump = L / N;
    if (biparabolic)
        jump = L / (2 * N);
        
    int i = node_id / (N + 1);
    int j = node_id % (N + 1);

    if (biparabolic) {
        // creating the combiantions to iterate over
        // combinations[k] = (gl_i, gl_j, x1, xj)
        std::vector<std::tuple<int, int, double, double>> combinations = {
            std::tuple<int, int, double, double>(i+0, j+0, -1, -1),
            std::tuple<int, int, double, double>(i+2, j+0, 1, -1),
            std::tuple<int, int, double, double>(i+0, j+2, -1, 1),
            std::tuple<int, int, double, double>(i+2, j+2, 1, 1),
            std::tuple<int, int, double, double>(i+1, j+0, 0, -1),
            std::tuple<int, int, double, double>(i+2, j+1, 1, 0),
            std::tuple<int, int, double, double>(i+0, j+1, -1, 0),
            std::tuple<int, int, double, double>(i+1, j+2, 0, 1),
            std::tuple<int, int, double, double>(i+1, j+1, 0, 0)
        };

        for (int k = 0; k < 9; k++) {
            int i = std::get<0>(combinations[k]);
            int j = std::get<1>(combinations[k]);
            double x1 = std::get<2>(combinations[k]);
            double x2 = std::get<3>(combinations[k]);

            x[i][j] += jump * i * _h(k+1, x1, x2);
            x[i][j] += jump * j * _h(k+1, x1, x2);
        }
    } else {
        // node 1
        x[i][j] += jump*i*_g(1, -1, -1);
        y[i][j] += jump*j*_g(1, -1, -1);

        // node 2
        x[i+1][j] += jump*(i+1)*_g(2, 1, -1);
        y[i+1][j] += jump*j*_g(2, 1, -1);

        // node 3
        x[i][j+1] += jump*i*_g(3, -1, 1);
        y[i][j+1] += jump*(j+1)*_g(3, -1, 1);

        // node 4
        x[i+1][j+1] += jump*(i+1)*_g(4, 1, 1);
        y[i+1][j+1] += jump*(j+1)*_g(4, 1, 1);
    }
}

void run_experiment(int N, double L, bool biparabolic) {
    int n;
    if (biparabolic) n = 2*N + 1;
    else n = N + 1;

    std::vector<std::vector<double>> x (
        std::vector<double>(0.0, n), n
    );
    std::vector<std::vector<double>> y (
        std::vector<double>(0.0, n), n
    );

    int nr_nodes = N^2;
    
}

