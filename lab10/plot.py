import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path


def get_snapshot_files(data_dir="data"):
    """
    Return snapshot files sorted by simulation time.
    """
    return sorted(
        Path(data_dir).glob("*.txt"),
        key=lambda p: float(p.stem)
    )


def load_snapshot(filename):
    """
    Load a snapshot from a text file.
    """
    return np.loadtxt(filename)


def save_selected_snapshots(
        t_values,
        data_dir="data",
        output_dir="snapshots",
        cmap="inferno"):
    """
    Save selected simulation times as PNG images.

    Parameters
    ----------
    t_values : list[float]
        Simulation times to save.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    for t in t_values:

        filename = Path(data_dir) / f"{t}00000.txt"

        if not filename.exists():
            print(f"Skipping t={t}: file {filename} not found")
            continue

        data = load_snapshot(filename)

        plt.figure(figsize=(6, 6))

        plt.imshow(
            data,
            origin="lower",
            cmap=cmap
        )

        plt.colorbar(label="u")
        plt.title(f"t = {t}")
        plt.xlabel("x")
        plt.ylabel("y")

        plt.tight_layout()

        output_file = output_dir / f"snapshot_t_{t}.png"
        plt.savefig(output_file, dpi=150)
        plt.close()

        print(f"Saved {output_file}")


def create_animation(
        data_dir="data",
        output_file="simulation.gif",
        step=4,
        fps=10,
        cmap="inferno"):
    """
    Create a GIF using every 'step'-th snapshot.
    """

    files = get_snapshot_files(data_dir)

    files = files[::step]

    if not files:
        raise RuntimeError("No snapshot files found.")

    frames = [load_snapshot(f) for f in files]

    vmin = min(frame.min() for frame in frames)
    vmax = max(frame.max() for frame in frames)

    fig, ax = plt.subplots(figsize=(6, 6))

    image = ax.imshow(
        frames[0],
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax
    )

    plt.colorbar(image, ax=ax, label="u")

    title = ax.set_title(
        f"t = {float(files[0].stem):.3f}"
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    def update(frame_index):

        image.set_array(frames[frame_index])

        t = float(files[frame_index].stem)

        title.set_text(f"t = {t:.3f}")

        return image, title

    animation = FuncAnimation(
        fig,
        update,
        frames=len(frames),
        interval=1000 / fps,
        blit=False
    )

    animation.save(
        output_file,
        writer="pillow",
        fps=fps
    )

    plt.close()

    print(f"Saved animation: {output_file}")

if __name__ == "__main__":
    # Save selected snapshots
    save_selected_snapshots(
        t_values=[1.0, 50.0, 100.0, 150.0, 200.0, 250.0]
    )

    # Create GIF using every 4th file
    create_animation(
        step=20,
        fps=12,
        output_file="./plots/simulation.gif"
    )