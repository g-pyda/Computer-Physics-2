from data_processing import visualize_grids, plot_energy_functional, visualize_residuals

def main():
    # -------------------------------
    # RELAXATION
    # experiment over the grid size for n in [2, 3, 4, 5, 6]
    # -------------------------------
    n_values = [2, 3, 4, 5, 6]
    max_iterations = 300

    plot_energy_functional(max_iterations, "relaxation", [2**n + 1 for n in n_values], f"./plots/relaxation/energy_functionals.png")
    visualize_grids("relaxation", "N", [2**n + 1 for n in n_values], f"./plots/relaxation/solution_grids.png")
    visualize_residuals("relaxation/65", "iter", [i for i in range(0, max_iterations, max_iterations // 10)], f"./plots/relaxation/residuals.png")

    # -------------------------------
    # OVERRELAXATION
    # experiment for n = 5 over the weight w in [0.5, 1.0, 1.5, 1.9]
    # -------------------------------
    weights = [1.0, 1.5, 1.9]
    max_iterations = 1000

    plot_energy_functional(max_iterations, "overrelaxation", weights, f"./plots/overrelaxation/energy_functionals.png")
    visualize_grids("overrelaxation", "w", weights, f"./plots/overrelaxation/solution_grids.png")
    visualize_residuals("overrelaxation/1.500000", "iter", [i for i in range(0, max_iterations, max_iterations // 10)], f"./plots/overrelaxation/residuals.png")

if __name__ == "__main__":
    main()