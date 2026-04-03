# ------------------------------
# DATA VISUALIZATION
# ------------------------------

import matplotlib.pyplot as plt
import numpy as np
import math

def visualize_grids(study, var_name, cases, output_filename):
    num_cases = len(cases)
    if num_cases == 0:
        print("Brak danych do wyświetlenia.")
        return

    cols = math.ceil(math.sqrt(num_cases))
    rows = math.ceil(num_cases / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), squeeze=False)
    axes = axes.flatten()

    for i, case in enumerate(cases):
        try:
            if study == "relaxation":
                data = np.loadtxt(f"./data/{study}/{case}/solution_grid.txt")
            elif study == "overrelaxation":
                data = np.loadtxt(f"./data/{study}/{case:.6f}/solution_grid.txt")
            im = axes[i].imshow(data, cmap='viridis', origin='lower')
            
            axes[i].set_title(f'Case: {var_name}={case}')
            axes[i].set_xlabel('X-axis')
            axes[i].set_ylabel('Y-axis')
            
            fig.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)
        except Exception as e:
            axes[i].set_title(f'Error: {case}')
            print(f"Nie udało się wczytać danych dla {case}: {e}")

    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()


def plot_energy_functional(max_iterations, study, cases, output_filename):
    import matplotlib.pyplot as plt
    import numpy as np

    exact_generated = False

    for case in cases:
        if study == "relaxation":
            data = np.loadtxt(f"./data/{study}/{case}/energy_functional.txt")
        elif study == "overrelaxation":
            data = np.loadtxt(f"./data/{study}/{case:.6f}/energy_functional.txt")
        iterations = data[:, 0]
        energy_values = data[:, 1]
        energy_exact_values = data[:, 2]
        energy_differences = np.abs(energy_values - energy_exact_values)

        plt.plot(iterations, energy_differences, label=f"EF_{case}")
        # if len(iterations) == max_iterations and not exact_generated:
        #     plt.plot(iterations, energy_differences, label='Exact Energy difference', linestyle='--')
        #     exact_generated = True

    plt.xlabel('Iteration')
    plt.yscale('log')
    plt.ylabel('Energy')
    plt.title('Energy Functional Convergence')
    plt.legend()
    plt.savefig(output_filename)
    plt.close()
        

def visualize_residuals(study, var_name, cases, output_filename):
    num_cases = len(cases)
    if num_cases == 0:
        print("Brak danych do wyświetlenia.")
        return

    cols = math.ceil(math.sqrt(num_cases))
    rows = math.ceil(num_cases / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), squeeze=False)
    axes = axes.flatten()

    for i, case in enumerate(cases):
        try:
            data = np.loadtxt(f"./data/{study}/residual_{case}.txt")
            im = axes[i].imshow(data, cmap='viridis', origin='lower')
            axes[i].set_title(f'Case: {var_name}={case}')
            axes[i].set_xlabel('X-axis')
            axes[i].set_ylabel('Y-axis')
            
            fig.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)
        except Exception as e:
            axes[i].set_title(f'Error: {case}')
            print(f"Nie udało się wczytać danych dla {case}: {e}")

    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()