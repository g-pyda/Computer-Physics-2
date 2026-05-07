import pandas as pd
import matplotlib.pyplot as plt
import os

# 1. Load the data from the C++ generated file
# Expected columns: time | norm | energy | avg_x | avg_y

script_dir = os.path.dirname(os.path.abspath(__file__))
try:
    cases = ['proper_init', 'inproper_init']
    for case in cases:
        os.makedirs(os.path.join(script_dir, 'plots', case))
        file_path = os.path.join(script_dir, 'data', case, 'deliverables.txt')
        columns = ['t', 'norm', 'energy', 'avg_x', 'avg_y']
        data = pd.read_csv(file_path, sep=r'\s+', header=None, names=columns)

        # Global plot styling
        plt.rcParams.update({'font.size': 10, 'grid.alpha': 0.3})

        # --- Plot 1: Norm vs Time ---
        # This checks the stability of the Crank-Nicolson scheme (should stay near 1.0)
        plt.figure(figsize=(8, 4))
        plt.plot(data['t'], data['norm'], color='blue', label='Norm')
        plt.xlabel(r'Time [t]')
        plt.ylabel(r'Norm $|\Psi|^2$')
        plt.title(r'Wavefunction Norm Stability over Time')
        plt.grid(True)
        plt.savefig(f"./plots/{case}/norm_vs_t.png", dpi=300)
        plt.close()

        # --- Plot 2: Mean Energy vs Time ---
        # In a conservative system, this should remain constant
        plt.figure(figsize=(8, 4))
        plt.plot(data['t'], data['energy'], color='red', label='Mean Energy')
        plt.xlabel(r'Time [t]')
        plt.ylabel(r'$\langle E \rangle$')
        plt.title(r'Mean Energy of the System')
        plt.grid(True)
        plt.savefig(f"./plots/{case}/energy_vs_t.png", dpi=300)
        plt.close()

        # --- Plot 3: Average X and Y vs Time ---
        # Shows the evolution of the real/imaginary averages (or coordinates)
        plt.figure(figsize=(8, 4))
        plt.plot(data['t'], data['avg_x'], label='Avg X (Real)', alpha=0.8)
        plt.plot(data['t'], data['avg_y'], label='Avg Y (Imag)', alpha=0.8)
        plt.xlabel(r'Time [t]')
        plt.ylabel(r'Expected Values')
        plt.title(r'Evolution of Average X and Y')
        plt.legend()
        plt.grid(True)
        plt.savefig(f"./plots/{case}/avg_xy_vs_t.png", dpi=300)
        plt.close()

        # --- Plot 4: Phase Trajectory (Avg Y vs Avg X) ---
        # Shows the correlation/trajectory of the system's state
        plt.figure(figsize=(6, 6))
        plt.plot(data['avg_x'], data['avg_y'], color='purple', lw=1)
        plt.scatter(data['avg_x'].iloc[0], data['avg_y'].iloc[0], color='green', label='Start')
        plt.scatter(data['avg_x'].iloc[-1], data['avg_y'].iloc[-1], color='orange', label='End')
        plt.xlabel(r'$\langle x \rangle$')
        plt.ylabel(r'$\langle y \rangle$')
        plt.title(r'Trajectory: $\langle y \rangle$ vs $\langle x \rangle$')
        plt.legend()
        plt.axis('equal') # Maintain 1:1 aspect ratio for the trajectory
        plt.grid(True)
        plt.savefig(f"./plots/{case}/trajectory_y_x.png", dpi=300)
        plt.close()

        print("Success: Plots generated (norm_vs_t.png, energy_vs_t.png, avg_xy_vs_t.png, trajectory_y_x.png)")

except Exception as e:
    print(f"An unexpected error occurred: {e}")