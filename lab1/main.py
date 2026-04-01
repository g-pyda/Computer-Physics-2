import time
from engine import gauss_seidel, create_grid

def main():
    # -------------------------------
    # RELAXATION
    # experiment over the grid size for n in [3, 4, 5, 6]
    # -------------------------------

    for n in range(3, 4):
        # power of 2 for the simulation: N = 2**n
        N = 2 ** n + 1

        # creating a grid
        grid, delta_x = create_grid(N)

        # relaxation through gauss-seidel method
        start_time = time.time()

        gauss_seidel(grid, N, delta_x, max_iterations=10, tolerance=1e-6)

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"n = {n} | Elapsed time: {elapsed_time:.6f} seconds")


if __name__ == "__main__":
    main()