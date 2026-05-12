import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import os


def visualize_and_save(N, basis_type="bilinear", L=5.0):
    data_path = f"./data/{basis_type}/{N}"
    plot_path = f"./plots/{basis_type}/{N}"
    
    os.makedirs(plot_path, exist_ok=True)
    
    try:
        nodes = np.loadtxt(os.path.join(data_path, "nodes.txt"))
        x, y = nodes[:, 1], nodes[:, 2]
        
        eigenvalues = np.loadtxt(os.path.join(data_path, "eigenvalues.txt"))
        energies = eigenvalues[:, 1]
        
        eigenvectors = np.loadtxt(os.path.join(data_path, "eigenvectors.txt"))
        
        physical_indices = np.where(energies > 0)[0]
        physical_energies = energies[physical_indices]
        physical_vectors = eigenvectors[:, physical_indices]
        
        sort_idx = np.argsort(physical_energies)
        physical_energies = physical_energies[sort_idx]
        physical_vectors = physical_vectors[:, sort_idx]

        def get_exact_energy(nx, ny, side_l):
            return ((nx**2 + ny**2) * np.pi**2) / (2 * side_l**2)

        exact_energies = sorted([get_exact_energy(nx, ny, L) for nx in range(1, 6) for ny in range(1, 6)])[:10]

        plt.figure(figsize=(10, 6))
        plt.plot(physical_energies[:10], 'ro', label="FEM")
        plt.plot(exact_energies, 'k+', markersize=10, label="EXACT")
        plt.title(f"The spectrum (N={N}, {basis_type})")
        plt.xlabel("State number")
        plt.ylabel("Energy E")
        plt.legend()
        plt.grid(True)
        
        plt.savefig(os.path.join(plot_path, "spectrum_comparison.png"))
        plt.close()

        tri = Triangulation(x, y)
        num_to_plot = min(10, len(physical_energies))

        for i in range(num_to_plot):
            plt.figure(figsize=(8, 7))
            psi = physical_vectors[:, i]
            energy = physical_energies[i]
            
            cnt = plt.tricontourf(tri, psi, levels=30, cmap='RdBu_r')
            plt.colorbar(cnt)
            plt.title(f"Wave function - State {i+1}\nE = {energy:.6f} (N={N}, {basis_type})")
            plt.gca().set_aspect('equal')
            
            plt.savefig(os.path.join(plot_path, f"wf_state_{i+1}.png"))
            plt.close()

        print(f"Plots saved at: {plot_path}")

    except Exception as e:
        print(f"Error occured: {basis_type} N={N}: {e}")

def get_exact_energy(nx, ny, L=5.0):
    return ((nx**2 + ny**2) * np.pi**2) / (2 * L**2)

def calculate_convergence():
    L = 5.0
    N_values = [4, 8, 16, 32] 
    basis_types = ["bilinear", "biparabolic"]
    
    exact_e1 = get_exact_energy(1, 1, L)
    
    plt.figure(figsize=(10, 7))
    
    for basis in basis_types:
        errors = []
        valid_N = []
        
        for N in N_values:
            path = f"./data/{basis}/{N}/eigenvalues.txt"
            if not os.path.exists(path):
                continue
                
            try:
                data = np.loadtxt(path)
                energies = data[:, 1]
                physical_energies = sorted(energies[energies > 0])
                
                if physical_energies:
                    calculated_e1 = physical_energies[0]
                    rel_error = abs(calculated_e1 - exact_e1) / exact_e1
                    errors.append(rel_error)
                    valid_N.append(N)
            except Exception as e:
                print(f"Error occurred while loading N={N} ({basis}): {e}")

        if errors:
            plt.loglog(valid_N, errors, 'o-', label=f"Basis: {basis}")


    if valid_N:
        x_ref = np.array(valid_N)
        plt.loglog(x_ref, 0.1/x_ref**2, '--', color='gray', alpha=0.5, label="O(1/N²)")
        plt.loglog(x_ref, 0.1/x_ref**4, ':', color='gray', alpha=0.5, label="O(1/N⁴)")

    plt.title("Convergence of the ground state (Convergence with N)")
    plt.xlabel("Side element number (N)")
    plt.ylabel("Relative error |(E_num - E_ex) / E_ex|")
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    
    os.makedirs("./plots", exist_ok=True)
    plt.savefig("./plots/convergence_analysis.png")
    print("Convergence plot saved at: ./plots/convergence_analysis.png")

if __name__ == "__main__":
    test_N = [4, 8, 16, 32]
    for n in test_N:
        visualize_and_save(n, "bilinear")
        visualize_and_save(n, "biparabolic")
    calculate_convergence()