import os
import matplotlib.pyplot as plt

base_dir = "./data"
output_dir = "./plots"
os.makedirs(output_dir, exist_ok=True)

folders = [
    "grnd_alpha_0.028125",
    "grnd_alpha_0.030000",
    "grnd_alpha_0.031000",
    "grnd_alpha_0.031250",
]

plt.figure(figsize=(10, 6))

# store all exact values by iteration
exact_dict = {}

for folder in folders:
    file_path = os.path.join(base_dir, folder, "grnd_energy.out")
    
    iters = []
    calc_vals = []

    with open(file_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            it, calc, exact = map(float, line.split())
            
            iters.append(it)
            calc_vals.append(calc)
            
            # store exact (avoid duplicates overriding inconsistently)
            if it not in exact_dict:
                exact_dict[it] = exact

    label_base = folder.replace("grnd_alpha_", "α=")
    plt.plot(iters, calc_vals, label=f"{label_base} (calc)")

# build full exact curve
exact_iters = sorted(exact_dict.keys())
exact_vals = [exact_dict[it] for it in exact_iters]

plt.plot(exact_iters, exact_vals, linestyle="--", color="black", label="exact")

plt.xlabel("Iteration")
plt.ylabel("Energy")
plt.title("Ground State Energy Convergence")
plt.legend()
plt.grid(True)

output_path = os.path.join(output_dir, "energy_plot.png")
plt.savefig(output_path, dpi=300)
plt.close()

print(f"Plot saved to {output_path}")