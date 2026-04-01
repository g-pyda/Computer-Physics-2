from data_processing import write_to_file, save_grid_to_file
from math import pi, sin

# -----------------------
# GRID CREATION PART
# -----------------------

def create_grid(N):
    grid = [[0 for _ in range(N)] for _ in range(N)]
    delta_x = 3 / (N - 1)
    for i in range(N):
        for j in range(N):
            if i == 0 or i == N - 1 or j == 0 or j == N - 1:
                grid[i][j] = 0
            else:
                grid[i][j] = poisson_equation(i * delta_x, j * delta_x)
    return grid, delta_x

# -----------------------
# RELAXATION PART
# -----------------------

energy_exact = -1*(3*pi/2)**2


def poisson_equation(x, y):
    return 2*pi**2 * sin(pi*x) * sin(pi*y)

def laplace_5_stencil(u, i, j, N, delta_x):
    if i == 0 or i == N - 1 or j == 0 or j == N - 1:
        return 0
    return -1 * (
        (u[i+1][j] - 2*u[i][j] + u[i-1][j]) / delta_x**2
        + (u[i][j+1] - 2*u[i][j] + u[i][j-1]) / delta_x**2)

    
def residuum(u, i, j, N, delta_x):
    res = (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4*u[i][j]) / delta_x**2 \
        + laplace_5_stencil(u, i, j, N, delta_x)
    return res  


def energy_functional(u, N, delta_x):
    energy = 0
    for i in range(1, N-1):
        for j in range(1, N-1):
            energy += 0.5*((u[i+1][j] - u[i-1][j]) / (2*delta_x))**2 \
                + 0.5*((u[i][j+1] - u[i][j-1]) / (2*delta_x))**2 \
                - laplace_5_stencil(u, i, j, N, delta_x) * u[i][j]
    return energy * delta_x**2


def gauss_seidel(u, N, delta_x, max_iterations=10, tolerance=1e-6):
    for iteration in range(max_iterations):
        max_diff = 0
        for i in range(1, N-1):
            for j in range(1, N-1):
                u[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) \
                    + (delta_x**2) * laplace_5_stencil(u, i, j, N, delta_x)

                energy = energy_functional(u, N, delta_x)
                diff = abs(energy - energy_exact)
                max_diff = max(max_diff, diff)
        # saving the residual (1/10th of the iterations)
        if iteration % (max_iterations // 10) == 0:
            for i in range(1, N-1):
                for j in range(1, N-1):
                    res = residuum(u, i, j, N, delta_x)
                    write_to_file(f"./data/{N}/residual_{iteration}.txt", f"{res:.6e}\n")
        
        # saving the energy functional
        write_to_file(f"./data/{N}/energy_functional.txt", f"{iteration} {energy:.6f} {energy_exact:.6f}\n")
        if max_diff < tolerance:
            print(f"Converged after {iteration+1} iterations.")
            break
    else:
        print("Max iterations reached without convergence.")

    save_grid_to_file(f"./data/{N}/solution_grid.txt", u, N)

