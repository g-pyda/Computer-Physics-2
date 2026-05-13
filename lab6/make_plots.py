import numpy as np
import matplotlib.pyplot as plt
import os

# Create the 'plots' directory if it doesn't exist
os.makedirs('plots', exist_ok=True)

def load_data(subpath):
    # Assumes the script is run in the parent directory of 'data'
    filepath = os.path.join('data', subpath, 'projections.txt')
    if not os.path.exists(filepath):
        print(f"Warning: File {filepath} not found")
        return None
    return np.loadtxt(filepath)

def plot_transitions(subpath, title, output_filename):
    data = load_data(subpath)
    if data is None:
        return
    
    t = data[:, 0]
    p0 = data[:, 1]
    p1 = data[:, 2]
    p2 = data[:, 3]
    leakage = data[:, 4]

    plt.figure(figsize=(10, 6))
    plt.plot(t, p0, label='Ground State $|p_0|^2$', linewidth=2)
    plt.plot(t, p1, label='1st Excited (deg 1) $|p_1|^2$', linewidth=2)
    plt.plot(t, p2, label='1st Excited (deg 2) $|p_2|^2$', linewidth=2, linestyle='--')
    plt.plot(t, leakage, label='Leakage', linewidth=2)
    
    plt.title(title, fontsize=14)
    plt.xlabel('Time [a.u.]', fontsize=12)
    plt.ylabel('Probability', fontsize=12)
    plt.ylim([-0.05, 1.05]) # Probability should be in [0, 1] range
    plt.legend()
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    
    # Save plot in the 'plots' directory
    save_path = os.path.join('plots', output_filename)
    plt.savefig(save_path, dpi=300)
    print(f"Saved plot: {save_path}")
    plt.close()

def plot_leakage_comparison():
    L = 5.0
    # Base factors used to calculate eF amplitudes
    factors = [0.1, 0.5, 1.0, 2.0, 5.0]
    
    plt.figure(figsize=(10, 6))
    
    for factor in factors:
        eF_val = factor / L
        
        # Formatting to match the exact C++ output: ef/0.02, ef/0.20, etc.
        folder_path = f"ef/{eF_val:.2f}"
        
        data = load_data(folder_path)
        if data is None:
            continue
            
        t = data[:, 0]
        leakage = data[:, 4]
        
        label = f'eF = {factor}/L ({eF_val:.2f})'
        if factor == 1.0:
            label += ' (Nominal)'
            
        plt.plot(t, leakage, label=label, linewidth=1.5)
        
    plt.title('Wave Function Leakage vs. Field Amplitude', fontsize=14)
    plt.xlabel('Time [a.u.]', fontsize=12)
    plt.ylabel('Leakage Probability', fontsize=12)
    plt.yscale('log')
    plt.legend()
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    
    # Save plot in the 'plots' directory
    save_path = os.path.join('plots', 'leakage_comparison.png')
    plt.savefig(save_path, dpi=300)
    print(f"Saved plot: {save_path}")
    plt.close()

if __name__ == "__main__":
    L = 5.0
    nominal_eF = 1.0 / L  # 0.20
    
    # 1. First case: nominal omega (resonance) -> eF = 0.20
    plot_transitions(f'ef/{nominal_eF:.2f}', 
                     r'Transition Probabilities (Resonant Frequency $\omega$)', 
                     'transitions_resonant.png')
    
    # 2. Second case: half omega 
    # (Assuming you saved it as 'data/half_omega/projections.txt'. If you saved it 
    # under the 'ef/' tree, change this path to match, e.g., 'ef/0.20_half')
    plot_transitions('half_omega', 
                     r'Transition Probabilities (Half Resonance $\omega/2$)', 
                     'transitions_half_omega.png')
    
    # 3. Third case: leakage analysis depending on field amplitude
    plot_leakage_comparison()