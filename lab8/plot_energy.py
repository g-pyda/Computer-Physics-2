import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def read_average_energy(file_path):
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                return float(line.split()[0])
            except ValueError:
                continue
    raise ValueError(f"Could not parse average energy from {file_path}")


def parse_omega_from_dirname(name):
    match = re.match(r"^omega_([0-9]+(?:\.[0-9]+)?)_pi$", name)
    if not match:
        raise ValueError(f"Invalid omega directory name format: {name}")
    factor = float(match.group(1))
    return factor * math.pi


def parse_damping_from_dirname(name):
    match = re.match(r"^damping_([0-9]+(?:\.[0-9]+)?)$", name)
    if not match:
        raise ValueError(f"Invalid damping directory name format: {name}")
    return float(match.group(1))


def parse_hybrid_from_dirname(name):
    match = re.match(r"^omega_([0-9]+(?:\.[0-9]+)?)_pi_damping_([0-9]+(?:\.[0-9]+)?)$", name)
    if not match:
        raise ValueError(f"Invalid hybrid directory name format: {name}")
    omega_factor = float(match.group(1))
    damping = float(match.group(2))
    return omega_factor * math.pi, damping


def collect_experiment_points(data_dir, parser):
    data_dir = Path(data_dir)
    if not data_dir.exists() or not data_dir.is_dir():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    points = []
    for entry in sorted(data_dir.iterdir()):
        if not entry.is_dir():
            continue
        average_file = entry / "average_energy.txt"
        if not average_file.exists():
            continue
        x_value = parser(entry.name)
        mean_energy = read_average_energy(average_file)
        points.append((x_value, mean_energy, entry.name))

    if not points:
        raise RuntimeError(f"No average_energy.txt files found under {data_dir}")
    return points


def plot_omega_experiment(data_dir, out_file, show_plot=False):
    points = collect_experiment_points(data_dir, parse_omega_from_dirname)
    points.sort(key=lambda item: item[0])
    omegas = [item[0] for item in points]
    mean_energies = [item[1] for item in points]

    plt.figure(figsize=(10, 6))
    plt.plot(omegas, mean_energies, marker="o", linestyle="-", color="tab:blue")
    plt.xlabel(r"Driving frequency $\omega$ [rad/s]")
    plt.ylabel(r"Mean energy $\langle E \rangle$")
    plt.title("Mean energy vs driving frequency")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    print(f"Saved omega experiment plot to: {out_path}")
    if show_plot:
        plt.show()
    plt.close()


def plot_damping_experiment(data_dir, out_file, show_plot=False):
    points = collect_experiment_points(data_dir, parse_damping_from_dirname)
    points.sort(key=lambda item: item[0])
    damping_values = [item[0] for item in points]
    mean_energies = [item[1] for item in points]

    plt.figure(figsize=(10, 3))
    plt.plot(damping_values, mean_energies, marker="o", linestyle="-", color="tab:orange")
    plt.xlabel(r"Damping coefficient $d$")
    plt.ylabel(r"Mean energy $\langle E \rangle$")
    plt.title("Mean energy vs damping coefficient")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    print(f"Saved damping experiment plot to: {out_path}")
    if show_plot:
        plt.show()
    plt.close()


def collect_hybrid_grid(data_dir):
    data_dir = Path(data_dir)
    if not data_dir.exists() or not data_dir.is_dir():
        raise FileNotFoundError(f"Hybrid data directory not found: {data_dir}")

    points = []
    omega_values = set()
    damping_values = set()
    for entry in sorted(data_dir.iterdir()):
        if not entry.is_dir():
            continue
        average_file = entry / "average_energy.txt"
        if not average_file.exists():
            continue
        omega, damping = parse_hybrid_from_dirname(entry.name)
        mean_energy = read_average_energy(average_file)
        points.append((omega, damping, mean_energy))
        omega_values.add(omega)
        damping_values.add(damping)

    if not points:
        raise RuntimeError(f"No average_energy.txt files found under {data_dir}")

    omega_values = sorted(omega_values)
    damping_values = sorted(damping_values)
    omega_index = {omega: idx for idx, omega in enumerate(omega_values)}
    damping_index = {d: idx for idx, d in enumerate(damping_values)}

    grid = [[float("nan") for _ in omega_values] for _ in damping_values]
    for omega, damping, mean_energy in points:
        grid[damping_index[damping]][omega_index[omega]] = mean_energy

    return omega_values, damping_values, grid


def read_energy_time(file_path):
    times = []
    energies = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                times.append(float(parts[0]))
                energies.append(float(parts[1]))
            except ValueError:
                continue

    if not times:
        raise ValueError(f"Could not parse energy time series from {file_path}")
    return np.array(times), np.array(energies)


def plot_energy_evolution_omega_pi(data_dir, out_file, show_plot=False):
    omega_dir = Path(data_dir) / "omega" / "omega_1.000000_pi"
    energy_file = omega_dir / "energy_time.txt"
    if not energy_file.exists():
        raise FileNotFoundError(f"Energy evolution file not found: {energy_file}")

    times, energies = read_energy_time(energy_file)

    plt.figure(figsize=(10, 6))
    plt.plot(times, energies, marker="", linestyle="-", color="tab:purple")
    plt.xlabel(r"Time $t$")
    plt.ylabel(r"Energy $E(t)$")
    plt.title(r"Energy evolution for $\omega = \pi$")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()

    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    print(f"Saved omega = pi energy evolution plot to: {out_path}")
    if show_plot:
        plt.show()
    plt.close()


def plot_hybrid_heatmap(data_dir, out_file, show_plot=False):
    omega_values, damping_values, grid = collect_hybrid_grid(data_dir)
    Z = np.array(grid, dtype=float)
    positive_mask = np.isfinite(Z) & (Z > 0)
    if not np.any(positive_mask):
        raise RuntimeError("Hybrid heatmap requires at least one positive mean energy value for logarithmic scaling.")
    vmin = float(np.nanmin(Z[positive_mask]))
    vmax = float(np.nanmax(Z[positive_mask]))

    plt.figure(figsize=(12, 8))
    mesh = plt.imshow(
        Z,
        origin="lower",
        aspect="auto",
        extent=[omega_values[0], omega_values[-1], damping_values[0], damping_values[-1]],
        cmap=cm.viridis,
        norm=LogNorm(vmin=vmin, vmax=vmax),
    )
    plt.colorbar(mesh, label=r"Mean energy $\langle E \rangle$")
    plt.xlabel(r"Driving frequency $\omega$ [rad/s]")
    plt.ylabel(r"Damping coefficient $d$")
    plt.title("Hybrid mean energy map: frequency vs damping (log scale)")
    plt.tight_layout()

    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    print(f"Saved hybrid heatmap to: {out_path}")
    if show_plot:
        plt.show()
    plt.close()


def plot_hybrid_surface(data_dir, out_file, show_plot=False):
    omega_values, damping_values, grid = collect_hybrid_grid(data_dir)
    X, Y = np.meshgrid(omega_values, damping_values)
    Z = np.array(grid, dtype=float)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        X,
        Y,
        Z,
        cmap=cm.viridis,
        edgecolor="k",
        linewidth=0.2,
        antialiased=True,
    )
    fig.colorbar(surf, shrink=0.6, aspect=12, label=r"Mean energy $\langle E \rangle$")
    ax.set_xlabel(r"Driving frequency $\omega$ [rad/s]")
    ax.set_ylabel(r"Damping coefficient $d$")
    ax.set_zlabel(r"Mean energy $\langle E \rangle$")
    ax.set_zscale("log")
    ax.set_title("Hybrid mean energy surface (log scale)")
    plt.tight_layout()

    out_path = Path(out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    print(f"Saved hybrid 3D surface to: {out_path}")
    if show_plot:
        plt.show()
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Plot mean energy experiments from lab8 simulation output.")
    parser.add_argument(
        "--data-root",
        default="data",
        help="Root lab8 data directory containing omega/, damping/, and hybrid/ subdirectories (default: data).",
    )
    parser.add_argument(
        "--output-root",
        default="plots",
        help="Directory where output image files will be written (default: plots).",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plots after generation.",
    )
    parser.add_argument(
        "--no-omega",
        action="store_true",
        help="Skip omega experiment plot.",
    )
    parser.add_argument(
        "--no-damping",
        action="store_true",
        help="Skip damping experiment plot.",
    )
    parser.add_argument(
        "--no-hybrid",
        action="store_true",
        help="Skip hybrid experiment plots.",
    )
    parser.add_argument(
        "--no-omega-pi-evolution",
        action="store_true",
        help="Skip the omega = pi energy evolution plot.",
    )
    args = parser.parse_args()

    root = Path(args.data_root)
    out_root = Path(args.output_root)

    if not args.no_omega:
        plot_omega_experiment(root / "omega", out_root / "mean_energy_vs_omega.png", show_plot=args.show)
    if not args.no_damping:
        plot_damping_experiment(root / "damping", out_root / "mean_energy_vs_damping.png", show_plot=args.show)
    if not args.no_omega_pi_evolution:
        plot_energy_evolution_omega_pi(root, out_root / "energy_evolution_omega_pi.png", show_plot=args.show)
    if not args.no_hybrid:
        plot_hybrid_heatmap(root / "hybrid", out_root / "hybrid_mean_energy_heatmap.png", show_plot=args.show)
        plot_hybrid_surface(root / "hybrid", out_root / "hybrid_mean_energy_surface.png", show_plot=args.show)


if __name__ == "__main__":
    main()
