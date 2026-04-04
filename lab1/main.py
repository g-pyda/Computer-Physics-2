from data_processing import visualize_grids, plot_energy_functional, visualize_residuals

def main():
    # -------------------------------
    # RELAXATION
    # experiment over the grid size for n in [2, 3, 4, 5, 6, 7]
    # -------------------------------
    n_values = [2, 3, 4, 5, 6, 7]
    max_iterations = 1000
    last_residuals = [
        "5/residual_17",
        "9/residual_32",
        "17/residual_72",
        "33/residual_122",
        "65/residual_576",
        "129/residual_900",
    ]

    plot_energy_functional("relaxation",
                           [2**n + 1 for n in n_values],
                           f"./plots/relaxation/energy_functionals_differences.png",
                           f"./plots/relaxation/energy_functionals.png")
    visualize_grids("relaxation",
                    "N",
                    [2**n + 1 for n in n_values],
                    f"./plots/relaxation/solution_grids.png")
    visualize_residuals("relaxation/129",
                        "iter",
                        [f"residual_{i}" for i in range(0, max_iterations, max_iterations // 10)],
                        f"./plots/relaxation/residuals_129.png")
    visualize_residuals("relaxation",
                        "N - last iter",
                        last_residuals,
                        "./plots/relaxation/residuals.png")
    
    print("Relaxation visualization ended")

    # -------------------------------
    # OVERRELAXATION - weight
    # experiment for n = 5 over the weight w in [0.5, 1.0, 1.5, 1.9]
    # -------------------------------
    weights = [w * 0.1 for w in range(5, 20)]

    plot_energy_functional("overrelaxation-w",
                           weights,
                           f"./plots/overrelaxation/w_energy_functionals_differences.png",
                           f"./plots/overrelaxation/w_energy_functionals.png",
                           ylim2=(-40, 0))

    # -------------------------------
    # OVERRELAXATION - size
    # experiment for n = 1.5 over the size n in [2, 3, 4, 5, 6, 7]
    # -------------------------------
    n_values = [2, 3, 4, 5, 6, 7]

    plot_energy_functional("overrelaxation-n",
                           [2**n + 1 for n in n_values],
                           f"./plots/overrelaxation/n_energy_functionals_differences.png",
                           f"./plots/overrelaxation/n_energy_functionals.png",
                           ylim2=(-40, 0))
    
    print("Overrelaxation visualization ended")
    
    # -------------------------------
    # GRID REFINEMENT
    # experiment for n = 2, 3, 4, 5, 6 over the refinement level in [0, 1, 2, 3, 4, 5]
    # -------------------------------

    n_values = [3, 4, 5, 6, 7, 8]
    last_residuals = [
        "9/residual_900",
        "17/residual_34",
        "33/residual_47",
        "65/residual_116",
        "129/residual_480",
        "257/residual_900",
    ]
    plot_energy_functional("grid_refinement",
                           [2**n + 1 for n in n_values],
                           f"./plots/grid_refinement/energy_functionals_differences.png",
                           f"./plots/grid_refinement/energy_functionals_close.png",
                           ylim2=(-60, 0))
    plot_energy_functional("grid_refinement",
                           [2**n + 1 for n in n_values],
                           f"./plots/grid_refinement/energy_functionals_differences.png",
                           f"./plots/grid_refinement/energy_functionals.png")
    visualize_grids("grid_refinement",
                    "N",
                    [2**n + 1 for n in n_values],
                    "./plots/grid_refinement/solution_grids.png")
    visualize_residuals("grid_refinement",
                        "last residual",
                        last_residuals,
                        "./plots/grid_refinement/residuals.png")
    
    print("Grid refinement visualization ended")

if __name__ == "__main__":
    main()