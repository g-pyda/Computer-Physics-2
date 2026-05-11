import numpy as np
import matplotlib.pyplot as plt
import os

# Set directory path where .out files are located
data_path = "./data/excitation/"

def plot_convergence(k_levels):
    plt.figure(figsize=(6, 10))
    for k in k_levels:
        # Load energy data: col 0 is iter, col 1 is E_calc, col 2 is E_exact
        file_path = os.path.join(data_path, f"{k}_energy.out")
        if k == 0: file_path = os.path.join(data_path, "grnd_energy.out")
        if os.path.exists(file_path):
            data = np.loadtxt(file_path)
            plt.plot(data[:, 0], data[:, 1], label=f"State K={k} (cal)")
            plt.axhline(y=data[0, 2], color='r', linestyle='--', alpha=0.3)
    
    plt.xlabel("Iteration Number")
    plt.ylabel("Energy Expectation Value")
    plt.title("Energy Convergence for Excited States")
    plt.legend()
    plt.grid(True)
    plt.savefig("./plots/excitation_conv.png")

def plot_convergence_cropped(k_levels):
    plt.figure(figsize=(6, 10))
    for k in k_levels:
        # Load energy data: col 0 is iter, col 1 is E_calc, col 2 is E_exact
        file_path = os.path.join(data_path, f"{k}_energy.out")
        if k == 0: file_path = os.path.join(data_path, "grnd_energy.out")
        if os.path.exists(file_path):
            data = np.loadtxt(file_path)
            plt.plot(data[:, 0], data[:, 1], label=f"State K={k} (cal)")
            plt.axhline(y=data[0, 2], color='r', linestyle='--', alpha=0.3)
    
    plt.xlabel("Iteration Number")
    plt.ylabel("Energy Expectation Value")
    plt.ylim((1.5, 6.0))
    plt.title("Energy Convergence for Excited States")
    plt.legend()
    plt.grid(True)
    plt.savefig("./plots/excitation_conv_cropped.png")

def plot_wavefunctions(k_levels):
    for k in k_levels:
        file_path = os.path.join(data_path, f"{k}_wf.out")
        if os.path.exists(file_path):
            # Load 2D grid data
            wf = np.loadtxt(file_path)
            
            plt.figure(figsize=(6, 5))
            # Plotting the heatmap of the wavefunction
            im = plt.imshow(wf, cmap='RdBu', interpolation='bilinear')
            plt.colorbar(im, label=r'$\Psi_K(x,y)$')
            plt.title(f"Visualization of Wavefunction State K={k}")
            plt.xlabel("Grid X")
            plt.ylabel("Grid Y")
            plt.savefig(f"./plots/excitation_{k}.png")


if __name__ == "__main__":
    # Specify the excited states to visualize
    levels = [0, 1, 2, 3, 4] 
    
    # Generate convergence plots
    plot_convergence(levels)
    plot_convergence_cropped(levels)
    
    # Generate wavefunction visualizations
    plot_wavefunctions(levels)