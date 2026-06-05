import os
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Define directories
    data_dir = 'data'
    plots_dir = 'plots'
    
    # Create the plots directory if it does not exist
    os.makedirs(plots_dir, exist_ok=True)
    
    # Data format in text files:
    # Column 0: Time
    # Column 1: Min_u
    # Column 2: Max_u
    # Column 3: Center_X
    # Column 4: Center_Y

    print("Starting data visualization...")

    # ---------------------------------------------------------
    # 1. Pure Diffusion Analysis (D=0.1, vx=vy=0)
    # ---------------------------------------------------------
    filename_diff = os.path.join(data_dir, 'results_pure_diffusion.txt')
    if os.path.exists(filename_diff):
        data = np.loadtxt(filename_diff)
        time = data[:, 0]
        min_u = data[:, 1]
        max_u = data[:, 2]
        
        plt.figure(figsize=(8, 6))
        plt.plot(time, max_u, label='Max value of u', color='red', linewidth=2)
        plt.plot(time, min_u, label='Min value of u', color='blue', linewidth=2)
        plt.title('Pure Diffusion: Min and Max values over Time')
        plt.xlabel('Time (t)')
        plt.ylabel('Value of u')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        
        out_path = os.path.join(plots_dir, 'pure_diffusion_min_max.png')
        plt.savefig(out_path, dpi=300)
        plt.close()
        print(f"Saved plot: {out_path}")
    else:
        print(f"Warning: Data file {filename_diff} not found.")

    # ---------------------------------------------------------
    # 2. Pure Advection Analysis (D=0, vx=vy=1)
    # ---------------------------------------------------------
    filename_adv = os.path.join(data_dir, 'results_pure_advection.txt')
    if os.path.exists(filename_adv):
        data = np.loadtxt(filename_adv)
        time = data[:, 0]
        max_u = data[:, 2]
        center_x = data[:, 3]
        center_y = data[:, 4]
        
        # Plot A: Peak Value over Time (Shape Constancy Check)
        # If the shape is perfectly constant, max_u should remain completely flat.
        plt.figure(figsize=(8, 6))
        plt.plot(time, max_u, label='Max value of u', color='green', linewidth=2)
        plt.title('Pure Advection: Peak Value over Time\n(Checking if shape remains constant)')
        plt.xlabel('Time (t)')
        plt.ylabel('Max value of u')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        
        out_path_adv1 = os.path.join(plots_dir, 'pure_advection_shape_check.png')
        plt.savefig(out_path_adv1, dpi=300)
        plt.close()
        print(f"Saved plot: {out_path_adv1}")
        
        # Plot B: Center Position vs Time
        plt.figure(figsize=(8, 6))
        plt.plot(time, center_x, label='Center X position', linestyle='-', color='purple')
        plt.plot(time, center_y, label='Center Y position', linestyle='--', color='orange')
        plt.title('Pure Advection: Center of Mass Trajectory')
        plt.xlabel('Time (t)')
        plt.ylabel('Coordinate Position')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        
        out_path_adv2 = os.path.join(plots_dir, 'pure_advection_trajectory.png')
        plt.savefig(out_path_adv2, dpi=300)
        plt.close()
        print(f"Saved plot: {out_path_adv2}")
    else:
        print(f"Warning: Data file {filename_adv} not found.")

    # ---------------------------------------------------------
    # 3. Advection-Diffusion Analysis (D=0.1, vx=vy=1)
    # ---------------------------------------------------------
    filename_adv_diff = os.path.join(data_dir, 'results_advection_diffusion.txt')
    if os.path.exists(filename_adv_diff):
        data = np.loadtxt(filename_adv_diff)
        time = data[:, 0]
        center_x = data[:, 3]
        center_y = data[:, 4]
        
        # Plot: Center Position vs Time (Velocity Check)
        # The slope of these lines represents the calculated velocity of the packet.
        plt.figure(figsize=(8, 6))
        plt.plot(time, center_x, label='Center X', color='teal', linewidth=2)
        plt.plot(time, center_y, label='Center Y', color='magenta', linewidth=2, linestyle=':')
        plt.title('Advection-Diffusion: Center Position vs Time\n(Velocity field verification)')
        plt.xlabel('Time (t)')
        plt.ylabel('Position')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        
        out_path_advdiff = os.path.join(plots_dir, 'advection_diffusion_velocity.png')
        plt.savefig(out_path_advdiff, dpi=300)
        plt.close()
        print(f"Saved plot: {out_path_advdiff}")
    else:
        print(f"Warning: Data file {filename_adv_diff} not found.")

    print("Visualization complete.")

if __name__ == "__main__":
    main()