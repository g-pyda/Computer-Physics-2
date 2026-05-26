#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <filesystem>
#include <iomanip>



// ============ CONSTANTS ============
const double L = 1.0;
const int N = 100;
const double vs = 1.0;
const double delta_x = L/N;

// ============ MAIN SIM ============

void run_simulation (
    std::string path,
    double omega,
    double damp_const
) {
    // parameters
    const double v = vs;
    const double d = damp_const; // damping constant 'd' in equation (2d * u_t)
    const double dt = delta_x/(2.0*v); // stability choice v*dt = dx/2
    const double t_final = 80.0;
    const int steps = static_cast<int>(std::ceil(t_final / dt));

    // grid and solution arrays
    std::vector<std::vector<double>> u_prev(N, std::vector<double>(N, 0.0));
    std::vector<std::vector<double>> u_cur(N, std::vector<double>(N, 0.0));
    std::vector<std::vector<double>> u_next(N, std::vector<double>(N, 0.0));

    int c = N/2; // center node

    // prepare output directory
    std::filesystem::path outdir = std::filesystem::path("data") / path;
    std::filesystem::create_directories(outdir);

    std::ofstream energy_file((outdir / "energy_time.txt").string());
    energy_file << std::setprecision(12);

    double energy_sum = 0.0;
    int energy_count = 0;

    // time-stepping loop
    for (int step = 0; step <= steps; ++step) {
        double t = step * dt;

        // compute u_next for interior points
        for (int i = 1; i < N-1; ++i) {
            for (int j = 1; j < N-1; ++j) {
                double lap = (u_cur[i+1][j] + u_cur[i-1][j] + u_cur[i][j+1] + u_cur[i][j-1] - 4.0*u_cur[i][j]) / (delta_x*delta_x);
                double coeff_cur = 2.0 - damp_const * dt;
                double coeff_prev = damp_const * dt - 1.0;
                double forcing = 0.0;
                if (i == c && j == c) forcing = std::sin(omega * t);
                u_next[i][j]    = u_cur[i][j] * coeff_cur 
                                + u_prev[i][j] * coeff_prev 
                                + (dt*dt * v * v) * lap
                                + (dt*dt) * forcing;
            }
        }

        // compute energy E(t) = 0.5 * "+" as discrete sum
        double E = 0.0;
        for (int i = 0; i < N-1; ++i) {
            for (int j = 0; j < N-1; ++j) {
                double u_t = (u_cur[i][j] - u_prev[i][j]) / dt;
                double u_x = (u_cur[i+1][j] - u_cur[i][j]) / (delta_x);
                double u_y = (u_cur[i][j+1] - u_cur[i][j]) / (delta_x);
                E += 0.5 * (u_t*u_t + u_x*u_x + u_y*u_y);
            }
        }
        E += 0.5 * std::pow((u_cur[N-1][N-1] - u_prev[N-1][N-1]) / dt, 2);
        for (int i = 0; i < N-1; ++i) {
            E += 0.5 * std::pow((u_cur[i][N-1] - u_prev[i][N-1]) / dt, 2);
            E += 0.5 * std::pow((u_cur[N-1][i] - u_prev[N-1][i]) / dt, 2);
        }
        
        E *= (delta_x * delta_x);

        // record energy
        energy_file << t << " " << E << "\n";

        // accumulate average over (10,80)
        if (t >= 10.0 && t <= t_final) {
            energy_sum += E;
            ++energy_count;
        }

        if (step % 1000 == 0) {
            std::cout << "t=" << t << " E=" << E << "\r" << std::flush;
        }

        // rotate time levels
        u_prev.swap(u_cur);
        u_cur.swap(u_next);
        // zero out u_next boundaries (maintain Dirichlet 0)
        for (int i = 0; i < N; ++i) {
            u_next[i][0] = 0.0; u_next[i][N-1] = 0.0;
            u_next[0][i] = 0.0; u_next[N-1][i] = 0.0;
        }
    }

    energy_file.close();

    // write average energy
    std::ofstream avgf((outdir / "average_energy.txt").string());
    double avgE = (energy_count>0) ? (energy_sum / energy_count) : 0.0;
    avgf << std::setprecision(12) << avgE << "\n";
    avgf.close();

    // save final snapshot
    std::ofstream sol((outdir / "final_grid.txt").string());
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            sol << u_cur[i][j] << (j+1==N?"":" ");
        }
        sol << "\n";
    }
    sol.close();
}