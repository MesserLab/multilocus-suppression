import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# make plot theme more like R
custom_params = {"font.family": "DejaVu Sans", "axes.facecolor": "white", "axes.edgecolor": "white",
    "axes.linewidth": 1.2, "axes.grid": True, "grid.color": "#EAEAEA", "grid.linestyle": "-",
    "grid.linewidth": 0.5, "xtick.color": "black", "ytick.color": "black", "xtick.major.size": 0,
    "ytick.major.size": 0}
sns.set_style("whitegrid", rc=custom_params)
plt.rcParams.update(custom_params)

# read the data
file = "/SSD/reg259/gene_drive_sims/multilocus-suppression/data/s-vs-n-r1-high.csv"
save_dir = "/SSD/reg259/gene_drive_sims/multilocus-suppression/figures-r1-high"
df = pd.read_csv(file)


def make_grid_edges(n_vals, s_vals):
    # cell edges for pcolormesh, in real data coordinates (so contour lines
    # computed from a formula line up with the heatmap cells)
    n_step = np.diff(n_vals).min()
    s_step = np.diff(s_vals).min()
    n_edges = np.concatenate([n_vals - n_step / 2, [n_vals[-1] + n_step / 2]])
    s_edges = np.concatenate([s_vals - s_step / 2, [s_vals[-1] + s_step / 2]])
    return n_edges, s_edges


def make_max_load_heatmap(save_dir):
    pivot = df.pivot(index="target_fitness_cost", columns="num_targets", values="avg_max_load")
    n_vals = pivot.columns.values.astype(float)
    s_vals = pivot.index.values.astype(float)
    n_edges, s_edges = make_grid_edges(n_vals, s_vals)

    fig, ax = plt.subplots(figsize=(6, 4.875))

    mesh = ax.pcolormesh(n_edges, s_edges, pivot.values, cmap="viridis", vmin=0, vmax=1, shading="flat")
    cbar = fig.colorbar(mesh, ax=ax, shrink=0.7, pad=0.08)
    cbar.ax.set_title("Avg max\ngenetic load", fontsize=10, pad=10)

    # theoretical maximum load, evaluated on a fine grid for smooth contour lines
    n_fine = np.linspace(n_vals.min(), n_vals.max(), 400)
    s_fine = np.linspace(s_vals.min(), s_vals.max(), 400)
    n, s = np.meshgrid(n_fine, s_fine)

    # theoretical max load = 1 - ((1 - drive_cost) * ((1 - s) ** n))
    drive_cost = df["drive_fitness_cost"].unique()[0]
    theoretical_max_load = 1 - ((1 - drive_cost) * ((1 - s) ** n))

    contour_levels = [0.5, 0.75, 0.95, 0.999]
    cs = ax.contour(n, s, theoretical_max_load, levels=contour_levels, colors="black",
        linewidths=1, linestyles="dashed")

    # place each label at the midpoint (by arc length) of its contour line
    label_positions = []
    for level_segs in cs.allsegs:
        if not level_segs:
            continue
        longest = max(level_segs, key=len)
        seg_lengths = np.hypot(*np.diff(longest, axis=0).T)
        cum_len = np.concatenate([[0], np.cumsum(seg_lengths)])
        midpoint = longest[np.searchsorted(cum_len, cum_len[-1] / 2)]
        label_positions.append(tuple(midpoint))

    level_fmt = {level: str(level) for level in contour_levels}
    labels = ax.clabel(cs, manual=label_positions, inline=True, fmt=level_fmt, fontsize=9)
    for label in labels:
        label.set_color("black")

    ax.set_xticks(n_vals)
    ax.set_yticks(s_vals[::2])
    ax.set_xlabel("Number of target sites (n)")
    ax.set_ylabel("Disrupted target site fitness cost (s)")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(f"{save_dir}/s-vs-n-max-load-heatmap.png", dpi=300)


def make_recovery_time_heatmap(save_dir):
    pivot = df.pivot(index="target_fitness_cost", columns="num_targets",
        values="avg_time_to_genetic_load_below_marker")
    n_vals = pivot.columns.values.astype(float)
    s_vals = pivot.index.values.astype(float)
    n_edges, s_edges = make_grid_edges(n_vals, s_vals)

    fig, ax = plt.subplots(figsize=(6, 4.875))

    # -1 marks parameter combos that theoretically never recover (genetic load
    # never falls back below the marker) -- mask these out to light gray
    plot_values = np.where(pivot.values == -1, np.nan, pivot.values)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#DCDCDC")

    mesh = ax.pcolormesh(n_edges, s_edges, plot_values, cmap=cmap, shading="flat")
    cbar = fig.colorbar(mesh, ax=ax, shrink=0.7, pad=0.08)
    cbar.ax.set_title("Avg time\nto recover", fontsize=10, pad=10)

    ax.set_xticks(n_vals)
    ax.set_yticks(s_vals[::2])
    ax.set_xlabel("Number of target sites (n)")
    ax.set_ylabel("Disrupted target site fitness cost (s)")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(f"{save_dir}/s-vs-n-time-to-recover-heatmap.png", dpi=300)


make_max_load_heatmap(save_dir)
make_recovery_time_heatmap(save_dir)
