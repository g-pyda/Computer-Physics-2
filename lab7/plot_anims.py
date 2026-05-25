import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def find_cycle_indices(Z_frames, N):
    # Track the center node's displacement to find exact phases of the first cycle
    center_idx = N // 2
    center_z = [frame[center_idx, center_idx] for frame in Z_frames]
    
    # Fallback default indices if a full cycle isn't completed in the data
    default_indices = [0, len(Z_frames)//4, len(Z_frames)//2, 3*len(Z_frames)//4, len(Z_frames)-1]
    
    if len(center_z) < 10:
        return default_indices

    # Find the first zero-crossing (going down) -> "Flat"
    idx_flat1 = None
    for i in range(1, len(center_z)):
        if center_z[i-1] > 0 and center_z[i] <= 0:
            idx_flat1 = i
            break
            
    if idx_flat1 is None: return default_indices
        
    # Peak Up is the maximum before the first zero-crossing
    idx_peak_up = np.argmax(center_z[:idx_flat1])
    
    # Find the second zero-crossing (going up) -> "Flat cycle end"
    idx_flat2 = None
    for i in range(idx_flat1 + 1, len(center_z)):
        if center_z[i-1] < 0 and center_z[i] >= 0:
            idx_flat2 = i
            break
            
    if idx_flat2 is None: return default_indices
        
    # Peak Down is the minimum between the two zero-crossings
    idx_peak_down = idx_flat1 + np.argmin(center_z[idx_flat1:idx_flat2])
    
    # Return indices for: Start, Peak Up, Flat 1, Peak Down, Flat 2
    return [0, idx_peak_up, idx_flat1, idx_peak_down, idx_flat2]

def generate_cycle_snapshots(X, Y, Z_frames, steps, exp_name, plots_dir, z_lim):
    N = X.shape[0]
    indices = find_cycle_indices(Z_frames, N)
    titles = ["T = 0 (Start)", "Peak Up", "Flat", "Peak Down", "Flat (Cycle End)"]
    
    # Increased width slightly to accommodate colorbars nicely
    fig = plt.figure(figsize=(28, 5))
    
    for i, idx in enumerate(indices):
        ax = fig.add_subplot(1, 5, i+1, projection='3d')
        ax.set_zlim(-z_lim, z_lim)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Displacement (Psi)")
        ax.set_title(f"{titles[i]}\nStep: {steps[idx]}")
        
        # Plot the surface for the specific phase
        surf = ax.plot_surface(X, Y, Z_frames[idx], cmap='viridis', edgecolor='none')
        
        # Lock the color scale to the global limits so colors are consistent across all 5 plots
        surf.set_clim(-z_lim, z_lim)
        
        # Add colorbar to each individual subplot
        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.1, label="Psi")
        
    plt.tight_layout()
    output_file = os.path.join(plots_dir, f"{exp_name}_snapshots.png")
    plt.savefig(output_file, dpi=300)
    plt.close(fig)
    print(f"  Saved cycle snapshots to: {output_file}")

def main():
    data_dir = "data"
    plots_dir = "plots"
    os.makedirs(plots_dir, exist_ok=True)
    
    experiments = [d for d in glob.glob(os.path.join(data_dir, "*")) if os.path.isdir(d)]
    
    if not experiments:
        print(f"No data found in the '{data_dir}' directory.")
        return

    for exp_path in experiments:
        exp_name = os.path.basename(exp_path)
        
        # Skip the convergence error directory if it exists
        if "convergence" in exp_name.lower():
            continue
            
        print(f"Processing experiment: {exp_name}...")
        
        files = glob.glob(os.path.join(exp_path, "step_*.dat"))
        if not files:
            print(f"  No data files found in {exp_path}. Skipping.")
            continue
            
        # Extract step numbers and sort files properly
        get_step = lambda f: int(os.path.basename(f).replace("step_", "").replace(".dat", ""))
        files.sort(key=get_step)
        steps = [get_step(f) for f in files]
        
        # Initialize grid
        data0 = np.loadtxt(files[0])
        N = int(np.sqrt(len(data0)))
        X = data0[:, 0].reshape((N, N))
        Y = data0[:, 1].reshape((N, N))
        
        print("  Loading frames into memory...")
        Z_frames = []
        for f in files:
            data = np.loadtxt(f)
            Z = data[:, 2].reshape((N, N))
            Z_frames.append(Z)
            
        z_min = np.min(Z_frames)
        z_max = np.max(Z_frames)
        z_abs_max = max(abs(z_min), abs(z_max))
        z_lim = z_abs_max * 1.1 if z_abs_max != 0 else 1.0

        # --- 1. GENERATE STATIC 5-PHASE SNAPSHOT ---
        print("  Generating 5-phase cycle snapshot...")
        generate_cycle_snapshots(X, Y, Z_frames, steps, exp_name, plots_dir, z_lim)

        # --- 2. GENERATE GIF ANIMATION ---
        print("  Generating animation...")
        fig_ani = plt.figure(figsize=(9, 6)) # slightly wider for the colorbar
        ax_ani = fig_ani.add_subplot(111, projection='3d')
        
        # Create an initial dummy surface to setup the colorbar once before the animation loop starts
        surf_init = ax_ani.plot_surface(X, Y, Z_frames[0], cmap='viridis', edgecolor='none')
        surf_init.set_clim(-z_lim, z_lim)
        fig_ani.colorbar(surf_init, ax=ax_ani, shrink=0.5, aspect=10, pad=0.1, label="Psi")
        
        def update(frame_idx):
            ax_ani.clear()
            ax_ani.set_zlim(-z_lim, z_lim)
            ax_ani.set_xlabel("X")
            ax_ani.set_ylabel("Y")
            ax_ani.set_zlabel("Displacement (Psi)")
            ax_ani.set_title(f"Experiment: {exp_name} | Step: {steps[frame_idx]}")
            
            surf = ax_ani.plot_surface(X, Y, Z_frames[frame_idx], cmap='viridis', edgecolor='none')
            # Enforce global color limits in every frame
            surf.set_clim(-z_lim, z_lim) 
            return surf,

        ani = animation.FuncAnimation(
            fig_ani, 
            update, 
            frames=len(Z_frames), 
            interval=50, 
            blit=False
        )
        
        output_gif = os.path.join(plots_dir, f"{exp_name}.gif")
        ani.save(output_gif, writer='pillow', fps=15)
        plt.close(fig_ani)
        
        print(f"  Successfully saved animation to: {output_gif}\n")

if __name__ == "__main__":
    main()