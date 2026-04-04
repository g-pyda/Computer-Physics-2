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
            elif study == "grid_refinement":
                data = np.loadtxt(f"./data/{study}/{case}/solution_grid.txt")
            im = axes[i].imshow(data, cmap='viridis', origin='lower')
            
            axes[i].set_title(f'Case: {var_name}={case:.1f}')
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


def plot_energy_functional(study, cases, output_filename1, output_filename2, itlog1=False, ylim2=(-200, 200), itlog2=False):
    import matplotlib.pyplot as plt
    import numpy as np

    if study == "grid_refinement":
        last_iter = 0
        for case in cases:
            data = np.loadtxt(f"./data/{study}/{case}/energy_functional.txt")

            if len(data.shape) == 1:
                iterations = np.array([data[0]]) + last_iter
                energy_values = np.array([data[1]])
                energy_exact_values = np.array([data[2]])
            else:
                iterations = data[:, 0] + last_iter
                #print(iterations, iterations[len(iterations)-1], last_iter)
                energy_values = data[:, 1]
                energy_exact_values = data[:, 2]
                last_iter = iterations[len(iterations)-1] + 1
            energy_differences = np.abs(energy_values - energy_exact_values)
            
            plt.plot(iterations, energy_differences, label=f"EF_{case:.1f}")

    else:
        for case in cases:
            if study == "relaxation":
                data = np.loadtxt(f"./data/{study}/{case}/energy_functional.txt")
            elif study.startswith("overrelaxation"):
                case_str = f"{case:.6f}" if isinstance(case, float) else str(case)
                data = np.loadtxt(f"./data/{study}/{case_str}/energy_functional.txt")

            if len(data.shape) == 1:
                iterations = np.array([data[0]])
                energy_values = np.array([data[1]])
                energy_exact_values = np.array([data[2]])
            else:
                iterations = data[:, 0]
                energy_values = data[:, 1]
                energy_exact_values = data[:, 2]
            energy_differences = np.abs(energy_values - energy_exact_values)

            plt.plot(iterations, energy_differences, label=f"EF_{case:.1f}")

    plt.xlabel('Iteration')
    plt.yscale('log')
    if itlog1:
        plt.xscale('log')
    plt.ylabel('Energy')
    plt.title("Energy Functional Convergence - | $S_{exact} - S_{numerical}$ |")
    plt.legend()
    plt.savefig(output_filename1)
    plt.close()


    max_iter = 0
    exact = 0
    if study == "grid_refinement":
        for case in cases:
            data = np.loadtxt(f"./data/{study}/{case}/energy_functional.txt")
    
            if len(data.shape) == 1:
                iterations = np.array([data[0]])
                energy_values = np.array([data[1]])
                energy_exact_values = np.array([data[2]])
            else:
                iterations = data[:, 0] + max_iter
                energy_values = data[:, 1]
                energy_exact_values = data[:, 2]
    
            if len(iterations) > max_iter:
                max_iter = len(iterations)
                exact = energy_exact_values[0]
            max_iter = int(iterations[-1] + 1)
    
            plt.plot(iterations, energy_values, label=f"EF_{case:.1f}")

    else:
        for case in cases:
            if study == "relaxation" or study == "grid_refinement":
                data = np.loadtxt(f"./data/{study}/{case}/energy_functional.txt")
            elif study == "overrelaxation":
                case_str = f"{case:.6f}" if isinstance(case, float) else str(case)
                data = np.loadtxt(f"./data/{study}/{case_str}/energy_functional.txt")

            if len(data.shape) == 1:
                iterations = np.array([data[0]])
                energy_values = np.array([data[1]])
                energy_exact_values = np.array([data[2]])
            else:
                iterations = data[:, 0]
                energy_values = data[:, 1]
                energy_exact_values = data[:, 2]

            if len(iterations) > max_iter:
                max_iter = len(iterations)
                exact = energy_exact_values[0]

            plt.plot(iterations, energy_values, label=f"EF_{case:.1f}")
    
    plt.plot(range(max_iter), [exact] * max_iter, label=f"EF_exact", linestyle='--')
    plt.xlabel('Iteration')
    plt.ylabel('Energy')
    if itlog2:
        plt.xscale('log')
    plt.ylim(ylim2)
    plt.title("Energy Functional Convergence - $S_{numerical}$")
    plt.legend()
    plt.savefig(output_filename2)
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
            data = np.loadtxt(f"./data/{study}/{case}.txt")
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