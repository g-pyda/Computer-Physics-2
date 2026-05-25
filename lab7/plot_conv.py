import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def main():
    # Define directories
    exp_dir = "data/convergence_omega_pi_L"
    plots_dir = "plots"
    os.makedirs(plots_dir, exist_ok=True)
    
    # Get all error data files
    files = glob.glob(os.path.join(exp_dir, "error_dt_*.dat"))
    
    if not files:
        print(f"No error files found in {exp_dir}.")
        return

    plt.figure(figsize=(10, 6))
    
    # Process and plot each file
    for f in sorted(files, reverse=True): # Sort descending to have 0.05 first in legend
        dt_value = f.split('_dt_')[-1].replace('.dat', '')
        
        # Read data using pandas
        data = pd.read_csv(f, sep=' ')
        plt.plot(data['time'], data['error'], label=f'dt = {dt_value}')

    # Formatting the plot
    plt.title("Numerical Error Over Time for Different Time Steps (dt)")
    plt.xlabel("Simulation Time (s)")
    plt.ylabel("Root Mean Square Error (RMSE)")
    plt.yscale("log") # Logarithmic scale highlights differences nicely
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    
    # Save the plot
    output_file = os.path.join(plots_dir, "convergence_plot.png")
    plt.savefig(output_file, dpi=300)
    print(f"Successfully saved convergence plot to: {output_file}")

if __name__ == "__main__":
    main()