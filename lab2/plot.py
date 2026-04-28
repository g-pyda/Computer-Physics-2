import numpy as np
import matplotlib.pyplot as plt
import os

# ========== SHAPE FUNCTIONS ========== #

def g_func(i, xi1, xi2):
    f1 = lambda x: (1.0 - x) / 2.0
    f2 = lambda x: (1.0 + x) / 2.0
    if i == 1: return f1(xi1) * f1(xi2)
    if i == 2: return f2(xi1) * f1(xi2)
    if i == 3: return f1(xi1) * f2(xi2)
    if i == 4: return f2(xi1) * f2(xi2)
    return 0

def h_func(i, xi1, xi2):
    q1 = lambda x: x * (x - 1.0) / 2.0
    q2 = lambda x: (1.0 - x) * (1.0 + x)
    q3 = lambda x: x * (x + 1.0) / 2.0
    funcs = [
        lambda x, y: q1(x)*q1(y), lambda x, y: q3(x)*q1(y),
        lambda x, y: q1(x)*q3(y), lambda x, y: q3(x)*q3(y),
        lambda x, y: q2(x)*q1(y), lambda x, y: q3(x)*q2(y),
        lambda x, y: q1(x)*q2(y), lambda x, y: q2(x)*q3(y),
        lambda x, y: q2(x)*q2(y)
    ]
    return funcs[i-1](xi1, xi2)

# ========== PHYSICS ========== #

def exact_potential(x, y):
    xc, yc = -1.0, 0.0
    r = np.sqrt((x - xc)**2 + (y - yc)**2)
    return -1.0 / (2.0 * np.pi) * np.log(r)

# ========== VISUALIZATION ENGINE ========== #

def plot_potential_comparison(folder, N, biparabolic=False, diff=True):
    path = f"./data/{folder}/{N}/"
    
    # loading data from C++ exports
    nodes = np.loadtxt(path + "nodes.txt")[:, 1:] 
    nlg = np.loadtxt(path + "nlg.txt", dtype=int)[:, 1:] - 1 
    psi_fem_raw = np.loadtxt(path + "solution.txt")
    psi_fem = psi_fem_raw[:, 1] if psi_fem_raw.ndim > 1 else psi_fem_raw

    # selecting shape functions based on basis type
    shape_fn = h_func if biparabolic else g_func
    num_local = 9 if biparabolic else 4

    # scanning elements with step 0.2
    xi_range = np.arange(-1.0, 1.1, 0.2)
    res = len(xi_range)
    
    # creating grids for potential mapping
    grid_size = N * (res - 1) + 1
    fem_grid = np.zeros((grid_size, grid_size))
    exact_grid = np.zeros((grid_size, grid_size))

    # iterating over elements to fill the grid
    for ej in range(N):
        for ei in range(N):
            k = ej * N + ei
            node_indices = nlg[k]
            node_potentials = psi_fem[node_indices]

            for jj, xi2 in enumerate(xi_range):
                for ii, xi1 in enumerate(xi_range):
                    # mapping from reference to real space (Eq. 4, 5)
                    x_real = sum(nodes[node_indices[i]][0] * g_func(i+1, xi1, xi2) for i in range(4))
                    y_real = sum(nodes[node_indices[i]][1] * g_func(i+1, xi1, xi2) for i in range(4))
                    
                    # calculating fem and exact potential values (Eq. 24)
                    val_fem = sum(node_potentials[i] * shape_fn(i+1, xi1, xi2) for i in range(num_local))
                    val_exact = exact_potential(x_real, y_real)

                    # filling grid matrices
                    gi = ei * (res - 1) + ii
                    gj = ej * (res - 1) + jj
                    fem_grid[gj, gi] = val_fem
                    exact_grid[gj, gi] = val_exact

    # calculating error grid
    diff_grid = fem_grid - exact_grid

    # plotting setup for two subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    extent = [0, 5, 0, 5]
    
    # first plot: always FEM result
    im1 = axes[0].imshow(fem_grid, origin='lower', extent=extent, cmap='viridis', interpolation='nearest')
    axes[0].set_title(f"FEM Potential (N={N}, {folder})")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.colorbar(im1, ax=axes[0])

    # second plot: based on diff flag
    if diff:
        im2 = axes[1].imshow(diff_grid, origin='lower', extent=extent, cmap='coolwarm', interpolation='nearest')
        axes[1].set_title("Difference (FEM - Exact)")
        file_suffix = "diff"
    else:
        im2 = axes[1].imshow(exact_grid, origin='lower', extent=extent, cmap='viridis', interpolation='nearest')
        axes[1].set_title("Exact Solution")
        file_suffix = "exact"
        
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    fig.colorbar(im2, ax=axes[1])

    # saving and closing figure
    output_path = f"plots/{folder}_N{N}_{file_suffix}.png"
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)
    
    print(f"Saved: {output_path}")
    return np.max(np.abs(diff_grid))

# ========== D.I.5 / D.II.2: ERROR VS N ========== #

def plot_error_convergence(N_values, biparabolic=False):
    errors = []
    folder = "biparabolic" if biparabolic else "bilinear"
    
    for N in N_values:
        try:
            max_err = plot_potential_comparison(folder, N, biparabolic)
            errors.append(max_err)
        except:
            print(f"Data for N={N} ({folder}) not found.")

    plt.figure()
    plt.plot(N_values[:len(errors)], errors, 'o-', label=folder)
    plt.xlabel("N (elements per side)")
    plt.ylabel("Max Absolute Error")
    plt.yscale('log')
    plt.xscale('log', base=2)
    plt.legend()
    plt.title("Convergence Analysis")
    plt.grid(True, which="both", ls="-")
    plt.savefig(f"plots/convergence_{folder}.png")

# ========== EXECUTION ========== #

if __name__ == "__main__":
    # For convergence plot (D.I.5)
    plot_potential_comparison("bilinear", 4, diff=False)
    plot_error_convergence([4, 8, 16, 32], biparabolic=False)
    plot_error_convergence([4, 8, 16, 32], biparabolic=True)
    