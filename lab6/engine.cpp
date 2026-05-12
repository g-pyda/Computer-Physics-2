#include "engine.h"
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <lapacke.h>
#include <complex>
#include <algorithm>

typedef std::complex<double> complexd;

// ======== SIMULATION CONSTANTS ========== //

const double    h_bar = 1.0;
const double    dt = 0.25;
const double    T_max = 500.0;
const int       N = 32;
const double    L = 5.0;
const double    eF = 1.0/L;

// ========== STIFFNESS MATRICES AND COMPONENTS ========== //

const std::vector<std::vector<double>> bilinear_stiffness_matrix = {
    { 0.666667, -0.166667, -0.333333, -0.166667},
    {-0.166667,  0.666667, -0.166667, -0.333333},
    {-0.333333, -0.166667,  0.666667, -0.166667},
    {-0.166667, -0.333333, -0.166667,  0.666667}
};

const std::vector<std::vector<double>> biparabolic_stiffness_matrix = {
    { 0.622222, -0.033333, -0.033333, -0.022222, -0.200000,  0.111111, -0.200000,  0.111111, -0.355556},
    {-0.033333,  0.622222, -0.022222, -0.033333, -0.200000, -0.200000,  0.111111,  0.111111, -0.355556},
    {-0.033333, -0.022222,  0.622222, -0.033333,  0.111111,  0.111111, -0.200000, -0.200000, -0.355556},
    {-0.022222, -0.033333, -0.033333,  0.622222,  0.111111, -0.200000,  0.111111, -0.200000, -0.355556},
    {-0.200000, -0.200000,  0.111111,  0.111111,  1.955556, -0.355556, -0.355556,  0.000000, -1.066667},
    { 0.111111, -0.200000,  0.111111, -0.200000, -0.355556,  1.955556,  0.000000, -0.355556, -1.066667},
    {-0.200000,  0.111111, -0.200000,  0.111111, -0.355556,  0.000000,  1.955556, -0.355556, -1.066667},
    { 0.111111,  0.111111, -0.200000, -0.200000,  0.000000, -0.355556, -0.355556,  1.955556, -1.066667},
    {-0.355556, -0.355556, -0.355556, -0.355556, -1.066667, -1.066667, -1.066667, -1.066667,  5.688889}
};

const std::vector<double> w = {
    (18.0 - std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 - std::sqrt(30.0)) / 36.0
};

const std::vector<double> gamma_pts = {
    -std::sqrt(3.0/7.0 + 2.0/7.0*std::sqrt(6.0/5.0)),
    -std::sqrt(3.0/7.0 - 2.0/7.0*std::sqrt(6.0/5.0)),
     std::sqrt(3.0/7.0 - 2.0/7.0*std::sqrt(6.0/5.0)),
     std::sqrt(3.0/7.0 + 2.0/7.0*std::sqrt(6.0/5.0))
};
// const std::vector<double> gamma = {
//     -1*sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)),
//     sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)),
//     sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0)),
//     -1*sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0))
// }


// ========== SHAPE FUNCTIONS ========== //

double g_func(int i, double xi1, double xi2) {
    auto f1 = [](double xi) { return (1.0 - xi) / 2.0; };
    auto f2 = [](double xi) { return (1.0 + xi) / 2.0; };
    if (i == 1) return f1(xi1) * f1(xi2);
    if (i == 2) return f2(xi1) * f1(xi2);
    if (i == 3) return f1(xi1) * f2(xi2);
    if (i == 4) return f2(xi1) * f2(xi2);
    return 0;
}

double h_func(int i, double xi1, double xi2) {
    auto q1 = [](double xi) { return xi * (xi - 1.0) / 2.0; };
    auto q2 = [](double xi) { return (1.0 - xi) * (1.0 + xi); };
    auto q3 = [](double xi) { return xi * (xi + 1.0) / 2.0; };
    switch(i) {
        case 1: return q1(xi1) * q1(xi2); case 2: return q3(xi1) * q1(xi2);
        case 3: return q1(xi1) * q3(xi2); case 4: return q3(xi1) * q3(xi2);
        case 5: return q2(xi1) * q1(xi2); case 6: return q3(xi1) * q2(xi2);
        case 7: return q1(xi1) * q2(xi2); case 8: return q2(xi1) * q3(xi2);
        case 9: return q2(xi1) * q2(xi2);
        default: return 0;
    }
}

// ========== EXACT SOLUTION AND BC ========== //

double get_exact_potential(double x, double y) {
    double xc = -1.0, yc = 0.0;
    double r = std::sqrt(std::pow(x - xc, 2) + std::pow(y - yc, 2));
    return -1.0 / (2.0 * M_PI) * std::log(r);
}

// ========== LINEAR SYSTEM SOLVER (LAPACKE) ========== //

void solve(
    std::vector<std::vector<double>>& H,
    std::vector<std::vector<double>>& O,
    std::vector<double>& E,
    int total_nodes
) {
    int n = total_nodes;
    int nrhs = 1;
    int lda = n;
    int ldb = n;

    // creating flat buffers for LAPACK
    std::vector<double> H_flat(n * n);
    std::vector<double> O_flat(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            H_flat[i * n + j] = H[i][j];
            O_flat[i * n + j] = O[i][j];
        }
    }

    // solving generalized symmetric-definite eigenproblem (eigenvalues in E, eigenvectors in H_flat)
    lapack_int info = LAPACKE_dsygv(
        LAPACK_ROW_MAJOR, 
        1,
        'V',
        'U',
        n,
        H_flat.data(),
        lda,
        O_flat.data(),
        ldb,
        E.data()
    );

    if (info > 0) {
        std::cerr << "error: matrix is singular at U(" << info << "," << info << ")" << std::endl;
        return;
    } else if (info < 0) {
        std::cerr << "error: illegal argument at position " << -info << " in lapack call" << std::endl;
        return;
    }

    // shifting the H_flat to original H
    for (int i = 0; i < n*n; i++)
        H[i/n][i%n] = H_flat[i];
}

// ========== POTENTIAL ASSEMBLY ========= //

void assemble_potential_matrix(
    int total_nodes,
    int N, double L,
    const std::vector<std::vector<int>>& nlg,
    std::vector<std::vector<double>>& V
) {
    double a = L / N;
    double jac = (a * a / 4.0);
    for (int k = 0; k < nlg.size(); ++k) {
        for (int hi = 0; hi < 4; ++hi) {
            for (int hj = 0; hj < 4; ++hj) {
                double wt = w[hi] * w[hj];
                // Calculation of physical x
                double xi1 = gamma_pts[hi];
                double xi2 = gamma_pts[hj];
                
                // x = x_start_el + (1 + xi1)*a/2
                double phys_x = (k % N) * a + (1.0 + xi1) * a / 2.0; 

                for (int i1 = 0; i1 < 9; ++i1) {
                    double phi_i = h_func(i1 + 1, xi1, xi2);
                    for (int i2 = 0; i2 < 9; ++i2) {
                        double phi_j = h_func(i2 + 1, xi1, xi2);
                        // V_ij = integral( eF * x * phi_i * phi_j )
                        V[nlg[k][i1]][nlg[k][i2]] += jac * wt * phys_x * phi_i * phi_j;
                    }
                }
            }
        }
    }
}

// ============= TIME EVOLUTION =============== //

void time_evolution(int total_nodes, std::vector<double>& c_initial, 
                    std::vector<std::vector<double>>& O,
                    std::vector<std::vector<double>>& H0, 
                    std::vector<std::vector<double>>& V,
                    std::string path,
                    double omega,
                    std::vector<double> psi0,
                    std::vector<double> psi1,
                    std::vector<double> psi2
                ) {
    
    int n = total_nodes;
    complexd form_const = complexd(0, -0.5*dt/h_bar);
    std::vector<complexd> c(n);
    std::vector<complexd> A_flat(n*n, 0.0);
    std::vector<complexd> B_flat(n, 0.0);
    std::vector<complexd> H_flat_t(n*n, 0.0);
    std::vector<complexd> H_flat_t_dt(n*n, 0.0);
    std::vector<int> ipiv(n);
    for(int i=0; i<n; ++i) c[i] = c_initial[i];

    std::ofstream f_out(path + "projections.txt");
    for (double t = 0; t < T_max; t += dt) {
        //  H(t) = H0 + eF * V * sin(omega * t) assembly
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                H_flat_t[i*n+j] = H0[i][j] + eF*V[i][j]*std::sin(omega*t);
                H_flat_t_dt[i*n+j] = H0[i][j] + eF*V[i][j]*std::sin(omega*(t+dt));
            }
        }
        // A and B matrices for equations: 
        // [O - (dt/2hi)H(t+dt)]c(t+dt) = [O + (dt/2hi)H(t)]c(t) 
        //          A           c(t+dt) = B
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A_flat[i*n+j] = O[i][j] - form_const*H_flat_t_dt[i*n+j];
                B_flat[i] += (O[i][j] + form_const*H_flat_t[i*n+j]) * c[j];
            }
        }
        
        // solution of the equation
        lapack_int info = LAPACKE_zgesv(
            LAPACK_ROW_MAJOR, 
            n,
            1,
            (lapack_complex_double*)A_flat.data(), 
            n,
            ipiv.data(), 
            (lapack_complex_double*)B_flat.data(), 
            1
        );

        if (info != 0) {
            std::cerr << "LAPACK error: " << info << std::endl;
            break;
        }

        // C actualization (solution is saved at B_flat)
        c = B_flat;

        // Monitoring: saving the norm and p_i
        complexd p0 = 0.0, p1 = 0.0, p2 = 0.0, norm_sq = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                norm_sq += std::conj(c[i]) * O[i][j] * c[j]; // N(t) = c* O c 
                p0 += psi0[i] * O[i][j] * c[j];  // p_i(t) = c_i^T O c(t) 
                p1 += psi1[i] * O[i][j] * c[j];
                p2 += psi2[i] * O[i][j] * c[j];
            }
        }

        // Format: t | norm | |p0|^2 | |p1|^2  |p2|^2 
        f_out << t << " " << std::real(norm_sq) << " " 
              << std::norm(p0) << " " << std::norm(p1) << " " << std::norm(p2) << "\n";
    }
    f_out.close();
    
}


// ========== MESH AND GLOBAL ASSEMBLY ========== //

double get_ground_state(
    std::vector<std::vector<double>>& H,
    std::vector<std::vector<double>>& O,
    std::vector<std::vector<int>>& nlg,
    std::vector<double>& c,
    int N, double L, bool biparabolic
) {
    double a = L / N;

    // calculating grid dimensions
    int nodes_per_side = biparabolic ? (2 * N + 1) : (N + 1);
    int total_nodes = nodes_per_side * nodes_per_side;
    double h_step = L / (nodes_per_side - 1);

    // generating global node coordinates
    std::vector<Node> p(total_nodes);
    for (int j = 0; j < nodes_per_side; ++j) {
        for (int i = 0; i < nodes_per_side; ++i) {
            p[j * nodes_per_side + i] = { i * h_step, j * h_step };
        }
    }

    // building mapping table nlg
    int step = biparabolic ? 2 : 1;
    for (int ej = 0; ej < N; ++ej) {
        for (int ei = 0; ei < N; ++ei) {
            int base = ej * step * nodes_per_side + ei * step;
            if (!biparabolic) {
                nlg.push_back({ base, base + 1, base + nodes_per_side, base + nodes_per_side + 1 });
            } else {
                nlg.push_back({
                    base, base + 2, base + 2 * nodes_per_side, base + 2 * nodes_per_side + 2,
                    base + 1, base + nodes_per_side + 2, base + nodes_per_side, base + 2 * nodes_per_side + 1,
                    base + nodes_per_side + 1
                });
            }
        }
    }

    // assembling global stiffness matrix H and overlap matric O
    std::vector<double> E(total_nodes, 0.0);
    const auto& local_S = biparabolic ? biparabolic_stiffness_matrix : bilinear_stiffness_matrix;

    for (int k = 0; k < (int)nlg.size(); ++k) {
        int m = (int)local_S.size();
        // assemble stiffness (H) as before
        for (int i1 = 0; i1 < m; ++i1) {
            for (int i2 = 0; i2 < m; ++i2) {
                H[nlg[k][i1]][nlg[k][i2]] += 0.5 * local_S[i1][i2];
            }
        }

        // assemble overlap (mass) matrix using 2D Gauss quadrature on reference square [-1,1]^2
        std::vector<std::vector<double>> local_O(m, std::vector<double>(m, 0.0));
        for (int gi = 0; gi < 4; ++gi) {
            for (int gj = 0; gj < 4; ++gj) {
                double wt = w[gi] * w[gj];
                for (int i1 = 0; i1 < m; ++i1) {
                    double phi_i = (m == 4) ? g_func(i1 + 1, gamma_pts[gi], gamma_pts[gj]) : h_func(i1 + 1, gamma_pts[gi], gamma_pts[gj]);
                    for (int i2 = 0; i2 < m; ++i2) {
                        double phi_j = (m == 4) ? g_func(i2 + 1, gamma_pts[gi], gamma_pts[gj]) : h_func(i2 + 1, gamma_pts[gi], gamma_pts[gj]);
                        local_O[i1][i2] += wt * phi_i * phi_j;
                    }
                }
            }
        }
        // account for Jacobian of mapping from reference to physical element: |J| = (a/2)*(a/2) = a*a/4
        double jac = a * a / 4.0;
        for (int i1 = 0; i1 < m; ++i1) {
            for (int i2 = 0; i2 < m; ++i2) {
                O[nlg[k][i1]][nlg[k][i2]] += jac * local_O[i1][i2];
            }
        }
    }

    // applying boundary conditions
    double eps = 1e-9;
    for (int i = 0; i < total_nodes; ++i) {
        if (p[i].x < eps || p[i].x > L - eps || p[i].y < eps || p[i].y > L - eps) {
            for (int j = 0; j < total_nodes; ++j)
                O[i][j] = O[j][i] = H[i][j] = H[j][i] = 0.0;
            
            H[i][i] = -1.0;
            O[i][i] = 1.0;
        }
    }

    // solving Hc = EOc
    solve(H, O, E, total_nodes);

    // returning the difference between two first energies
    return E[1] - E[0];
}

// ========== MAIN SIMULATION ========== //

void run_simulation (
    std::string path,
    bool half_omega
) {
    std::string folder = "./data/" + path + "/";
    std::filesystem::create_directories(folder);

    int total_nodes = (2 * N + 1) * (2 * N + 1);
    std::vector<std::vector<double>> H(total_nodes, std::vector<double>(total_nodes, 0.0));
    std::vector<std::vector<double>> O(total_nodes, std::vector<double>(total_nodes, 0.0));
    std::vector<std::vector<double>> V(total_nodes, std::vector<double>(total_nodes, 0.0));
    std::vector<double> c(total_nodes, 0.0);
    std::vector<std::vector<int>> nlg;

    // getting the ground state (biparabolic approach)
    double omega = get_ground_state(H, O, nlg, c, N, L, true);
    std::vector<double> psi0(total_nodes);
    std::vector<double> psi1(total_nodes);
    std::vector<double> psi2(total_nodes);

    for (int i = 0; i < total_nodes; ++i) {
        psi0[i] = H[i][0];
        psi1[i] = H[i][1];
        psi2[i] = H[i][2];
    }

    // setting up the omega
    omega /= h_bar;
    if (half_omega) omega /= 2.0;

    // assembling the potential matrix
    assemble_potential_matrix(total_nodes, N, L, nlg, V);

    // performing the time evolution
    time_evolution(total_nodes, c, O, H, V, path, omega, psi0, psi1, psi2);
}

/*
 * INSTRUKCJA IMPLEMENTACJI EWOLUCJI CZASOWEJ (CRANK-NICOLSON)
 * KROK PO KROKU:
 *
 * 1. ZMIANA STRUKTURY DANYCH W SOLVE_AND_SAVE:
 * - Obecnie solve_and_save zapisuje wyniki do pliku. Musisz ją zmodyfikować tak, 
 * aby zwracała (np. przez referencję) wektory własne (H_flat po wykonaniu LAPACKE_dsygv).
 * - Będą one potrzebne jako:
 * a) Stan początkowy c(t=0) - jest nim wektor własny dla najniższej energii[cite: 10].
 * b) Baza do obliczania rzutów p_i(t) = <Psi_i | Psi(t)>[cite: 35].
 *
 * 2. PARAMETRYZACJA REZONANSU (w run_simulation):
 * - Po rozwiązaniu problemu własnego (H0), odczytaj dwie najniższe energie E0 i E1.
 * - Oblicz pulsację rezonansową: omega = (E1 - E0) / h_bar[cite: 34].
 * - Ustaw amplitudę pola eF = 0.5 / L[cite: 34].
 *
 * 3. ASAMBLACJA MACIERZY POTENCJAŁU V:
 * - Użyj funkcji assemble_potential_matrix, którą już masz, ale upewnij się, 
 * że obsługuje ona przypadek biparaboliczny (obecnie w kodzie masz 'for i1 < 4').
 * - Macierz V reprezentuje operator eFx (bez członu sin(wt))[cite: 13, 32].
 *
 * 4. IMPLEMENTACJA PĘTLI CZASU (w time_evolution):
 * - Przygotuj wektor zespolony c typu std::complex<double> o rozmiarze total_nodes.
 * - Zainicjalizuj c[i] wartościami pierwszego wektora własnego (ground state).
 * - W każdej iteracji (t += dt):
 *
 * A. ZBUDUJ MACIERZ UKŁADU (Lewa strona Crank-Nicolson):
 * M_left = S - (dt / (2 * i * h_bar)) * H(t + dt) [cite: 31]
 * gdzie H(t + dt) = H0 + V * sin(omega * (t + dt))[cite: 6, 32].
 * Uwaga: S w równaniu (5) to Twoja macierz Overlap O[cite: 31, 35].
 *
 * B. OBLICZ WEKTOR PRAWEJ STRONY (Prawa strona Crank-Nicolson):
 * M_right = S + (dt / (2 * i * h_bar)) * H(t) [cite: 31]
 * RHS = M_right * c(t)
 *
 * C. ROZWIĄŻ UKŁAD RÓWNAŃ:
 * M_left * c(t + dt) = RHS
 * Użyj funkcji LAPACKE_zgesv (dedykowana dla układów liniowych zespolonych).
 *
 * 5. MONITOROWANIE I ZAPIS WYNIKÓW:
 * - W każdym kroku obliczaj:
 * a) Normę: N(t) = c*^T * O * c. Powinna być stała (bliska 1.0)[cite: 35].
 * b) Prawdopodobieństwa stanów: |p_i(t)|^2 = |c_i^T * O * c(t)|^2, 
 * gdzie c_i to i-ty stacjonarny wektor własny[cite: 35, 37].
 * - Zapisz t oraz |p_i(t)|^2 do pliku tekstowego, aby móc wygenerować wykresy[cite: 37].
 *
 * 6. REALIZACJA DODATKOWYCH SCENARIUSZY (Deliverables):
 * - Wykonaj symulację dla nominalnego omega (rezonans)[cite: 37].
 * - Powtórz symulację zmieniając omega na 0.5 * omega[cite: 38].
 * - Zaobserwuj "leakage" (wyciek) do wyższych stanów wzbudzonych przy różnych eF[cite: 39].
 */