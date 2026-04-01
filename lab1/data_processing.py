import os

def write_to_file(filename, content):
    dir = os.path.dirname(filename)
    if not os.path.exists(dir):
        os.makedirs(dir)

    with open(filename, 'w') as file:
        file.write(content)

def save_grid_to_file(filename, grid, N):
    dir = os.path.dirname(filename)
    if not os.path.exists(dir):
        os.makedirs(dir)

    with open(filename, 'w') as file:
        for row in range(N):
            for col in range(N):
                file.write(f"{grid[row][col]:.6f} ")
            file.write("\n")
        