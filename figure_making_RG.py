import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

# make plot theme more like R
custom_params = {"font.family": "DejaVu Sans", "axes.facecolor": "white", "axes.edgecolor": "white",
    "axes.linewidth": 1.2, "axes.grid": True, "grid.color": "#EAEAEA", "grid.linestyle": "-",
    "grid.linewidth": 0.5, "xtick.color": "black", "ytick.color": "black", "xtick.major.size": 0,
    "ytick.major.size": 0}
sns.set_style("whitegrid", rc=custom_params)
plt.rcParams.update(custom_params)


def add_panel_label(fig, ax, label, dx_in=-0.5, dy_in=0.08):
    # anchored to each axes' own bounding box with a constant *absolute* (inch)
    # offset, rather than an axes-fraction offset -- axes-fraction x/y scales with
    # that axes' own width/height, so panels of different widths (e.g. heatmaps
    # with a colorbar carving out space vs. plain trajectory plots) end up with
    # visibly different absolute offsets even when given the "same" -0.1 fraction
    pos = ax.get_position()
    x = pos.x0 + dx_in / fig.get_figwidth()
    y = pos.y1 + dy_in / fig.get_figheight()
    fig.text(x, y, label, fontsize=18, fontweight="bold", va="bottom", ha="left")


def make_grid_edges(n_vals, s_vals):
    # cell edges for pcolormesh, in real data coordinates (so contour lines
    # computed from a formula line up with the heatmap cells)
    n_step = np.diff(n_vals).min()
    s_step = np.diff(s_vals).min()
    n_edges = np.concatenate([n_vals - n_step / 2, [n_vals[-1] + n_step / 2]])
    s_edges = np.concatenate([s_vals - s_step / 2, [s_vals[-1] + s_step / 2]])
    return n_edges, s_edges


def plot_max_load_heatmap(ax, df, cbar_pad=0.08, cbar_fraction=0.15):
    pivot = df.pivot(index="target_fitness_cost", columns="num_targets", values="avg_max_load")
    n_vals = pivot.columns.values.astype(float)
    s_vals = pivot.index.values.astype(float)
    n_edges, s_edges = make_grid_edges(n_vals, s_vals)

    mesh = ax.pcolormesh(n_edges, s_edges, pivot.values, cmap="viridis", vmin=0, vmax=1, shading="flat")
    cbar = ax.figure.colorbar(mesh, ax=ax, shrink=0.7, pad=cbar_pad, fraction=cbar_fraction)
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


def make_max_load_heatmap(df, save_dir):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    plot_max_load_heatmap(ax, df)
    fig.tight_layout()
    fig.savefig(f"{save_dir}/s-vs-n-max-load-heatmap.png", dpi=450)


def plot_recovery_time_heatmap(ax, df, cbar_pad=0.08, cbar_fraction=0.15):
    pivot = df.pivot(index="target_fitness_cost", columns="num_targets",
        values="avg_time_to_genetic_load_below_marker")
    n_vals = pivot.columns.values.astype(float)
    s_vals = pivot.index.values.astype(float)
    n_edges, s_edges = make_grid_edges(n_vals, s_vals)

    # -1 marks parameter combos that theoretically never recover (genetic load
    # never falls back below the marker) -- mask these out to light gray
    plot_values = np.where(pivot.values == -1, np.nan, pivot.values)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#DCDCDC")

    mesh = ax.pcolormesh(n_edges, s_edges, plot_values, cmap=cmap, shading="flat")
    cbar = ax.figure.colorbar(mesh, ax=ax, shrink=0.7, pad=cbar_pad, fraction=cbar_fraction)
    cbar.ax.set_title("Avg time\nto recover", fontsize=10, pad=10)

    ax.set_xticks(n_vals)
    ax.set_yticks(s_vals[::2])
    ax.set_xlabel("Number of target sites (n)")
    ax.set_ylabel("Disrupted target site fitness cost (s)")
    ax.grid(False)


def make_recovery_time_heatmap(df, save_dir):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    plot_recovery_time_heatmap(ax, df)
    fig.tight_layout()
    fig.savefig(f"{save_dir}/s-vs-n-time-to-recover-heatmap.png", dpi=450)


def load_avg_trajectory(data_dir, n, s):
    # each replicate stops once it reaches an absorbing state (drive lost/fixed),
    # so its final genetic load holds steady for any generation beyond that --
    # forward-fill each replicate out to the longest-running replicate before
    # averaging, rather than letting the average silently drop replicates
    folder = os.path.join(data_dir, f"s{s}_n{n}")
    rep_files = sorted(glob.glob(os.path.join(folder, "*.csv")))

    reps = []
    for rep_file in rep_files:
        rep_df = pd.read_csv(rep_file, skipinitialspace=True)
        reps.append(rep_df.set_index("gen_num")["genetic_load"])

    max_gen = int(max(rep.index.max() for rep in reps))
    full_index = np.arange(1, max_gen + 1)
    reps = [rep.reindex(full_index).ffill() for rep in reps]

    return pd.concat(reps, axis=1).mean(axis=1)


def plot_trajectory_lines(ax, df, data_dir, param_combos, colors):
    # param_combos: list of (n, s) tuples, one per line, in the same order as colors
    drive_cost = df["drive_fitness_cost"].unique()[0]

    for (n, s), color in zip(param_combos, colors):
        avg_load = load_avg_trajectory(data_dir, n, s)
        ax.plot(avg_load.index, avg_load.values, color=color, linewidth=1.5, zorder=3)

        theoretical_max_load = 1 - ((1 - drive_cost) * ((1 - s) ** n))
        ax.axhline(theoretical_max_load, color=color, linestyle="dashed", linewidth=1.2, zorder=2)

    ax.set_xlabel("Generation")
    ax.set_ylabel("Average genetic load")
    ax.set_xlim(0, 100)


def add_side_trajectory_legends(ax, colors, value_labels, fixed_label=None, inside=False):
    # anchored at the axes' own top-right/top-left corner (rather than a fixed
    # x fraction past the right edge) so a wider legend -- e.g. long "n = 10,
    # s = 0.25" value labels -- grows back into the axes instead of overflowing
    # past the figure boundary and getting clipped, which a right-edge anchor
    # would let happen once labels get longer than the short single-parameter case
    loc = "upper right" if inside else "upper left"
    anchor_x = 0.98 if inside else 1.02

    # three separate legends stacked to the right of the axes: varying-parameter
    # colors, the fixed-parameter value, and the realized/potential line-style key
    # -- when nothing is held fixed (e.g. both n and s vary together), fixed_label
    # is omitted and the value legend drops down to take its place
    value_handles = [Line2D([0], [0], color=color, linewidth=1.5) for color in colors]
    value_legend = ax.legend(value_handles, value_labels, loc=loc, bbox_to_anchor=(anchor_x, 1.0),
        frameon=True, framealpha=1.0, handlelength=1.5)
    ax.add_artist(value_legend)

    style_y = 0.7
    if fixed_label is not None:
        fixed_legend = ax.legend([Line2D([0], [0], color="none")], [fixed_label], loc=loc,
            bbox_to_anchor=(anchor_x, 0.8), frameon=True, framealpha=1.0, handlelength=0, handletextpad=0)
        ax.add_artist(fixed_legend)
    else:
        style_y = 0.8

    style_handles = [Line2D([0], [0], color="black", linestyle="solid"),
        Line2D([0], [0], color="black", linestyle="dashed")]
    style_labels = ["realized", "potential"]
    ax.legend(style_handles, style_labels, loc=loc, bbox_to_anchor=(anchor_x, style_y), frameon=True, framealpha=1.0)


def make_vary_n_fixed_s_trajectories(df, data_dir, save_dir):
    n_vals = [1, 5, 10]
    s = 0.3
    colors = ["#2a78d6", "#1baf7a", "#eda100"]

    fig, ax = plt.subplots(figsize=(7.5, 5))
    plot_trajectory_lines(ax, df, data_dir, [(n, s) for n in n_vals], colors)
    add_side_trajectory_legends(ax, colors, [f"n = {n}" for n in n_vals], f"s = {s}")

    fig.tight_layout()
    fig.subplots_adjust(right=0.78)
    fig.savefig(f"{save_dir}/vary-n-fixed-s-trajectories.png", dpi=450)


def make_vary_s_fixed_n_trajectories(df, data_dir, save_dir):
    n = 10
    s_vals = [0.1, 0.2, 0.3]
    colors = ["#2a78d6", "#1baf7a", "#eda100"]

    fig, ax = plt.subplots(figsize=(7.5, 5))
    plot_trajectory_lines(ax, df, data_dir, [(n, s) for s in s_vals], colors)
    add_side_trajectory_legends(ax, colors, [f"s = {s}" for s in s_vals], f"n = {n}")

    fig.tight_layout()
    fig.subplots_adjust(right=0.78)
    fig.savefig(f"{save_dir}/vary-s-fixed-n-trajectories.png", dpi=450)


def make_vary_n_and_s_trajectories(df, data_dir, save_dir):
    param_combos = [(10, 0.25), (1, 0.95)]
    colors = ["#2a78d6", "#eda100"]

    fig, ax = plt.subplots(figsize=(7.5, 5))
    plot_trajectory_lines(ax, df, data_dir, param_combos, colors)
    value_labels = [f"n = {n}\ns = {s}" for n, s in param_combos]
    add_side_trajectory_legends(ax, colors, value_labels)

    fig.tight_layout()
    fig.subplots_adjust(right=0.78)
    fig.savefig(f"{save_dir}/vary-n-and-s-trajectories.png", dpi=450)


def plot_percent_of_potential_heatmap(ax, df, cbar_pad=0.08, cbar_fraction=0.15):
    # percent of the theoretical maximum genetic load actually realized
    # (avg_max_load / theoretical_potential_load), already precomputed as
    # fraction_of_max_load_attained -- pinning color scale to 0-1 across both
    # parameter sets keeps panels comparable even though their data ranges differ
    pivot = df.pivot(index="target_fitness_cost", columns="num_targets", values="fraction_of_max_load_attained")
    n_vals = pivot.columns.values.astype(float)
    s_vals = pivot.index.values.astype(float)
    n_edges, s_edges = make_grid_edges(n_vals, s_vals)

    mesh = ax.pcolormesh(n_edges, s_edges, pivot.values, cmap="viridis", vmin=0, vmax=1, shading="flat")
    cbar = ax.figure.colorbar(mesh, ax=ax, shrink=0.7, pad=cbar_pad, fraction=cbar_fraction)
    cbar.ax.set_title(r"$\frac{\text{realized}}{\text{potential}}$", fontsize=14, pad=16)

    ax.set_xticks(n_vals)
    ax.set_yticks(s_vals[::2])
    ax.set_xlabel("Number of target sites (n)")
    ax.set_ylabel("Disrupted target site fitness cost (s)")
    ax.grid(False)


def make_percent_of_potential_heatmap(df, save_dir):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    plot_percent_of_potential_heatmap(ax, df)
    fig.tight_layout()
    fig.savefig(f"{save_dir}/s-vs-n-percent-of-potential-heatmap.png", dpi=450)


def make_combined_figure(df, data_dir, save_dir):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    colors = ["#2a78d6", "#1baf7a", "#eda100"]

    # a tighter colorbar (less pad/fraction reserved for it) leaves more of each
    # heatmap's own axes for the pcolormesh itself, widening A/B; this doesn't
    # touch C/D since they have no colorbar, and doesn't touch tight_layout's
    # own spacing since that's set by the pad= below, unchanged
    plot_max_load_heatmap(axes[0, 0], df, cbar_fraction=0.1)
    plot_recovery_time_heatmap(axes[0, 1], df, cbar_fraction=0.1)

    plot_trajectory_lines(axes[1, 0], df, data_dir, [(1, 0.3), (5, 0.3), (10, 0.3)], colors)
    add_side_trajectory_legends(axes[1, 0], colors, ["n = 1", "n = 5", "n = 10"], "s = 0.3", inside=True)

    plot_trajectory_lines(axes[1, 1], df, data_dir, [(10, 0.1), (10, 0.2), (10, 0.3)], colors)
    add_side_trajectory_legends(axes[1, 1], colors, ["s = 0.1", "s = 0.2", "s = 0.3"], "n = 10", inside=True)

    # rect reserves a small top margin (uniformly, so relative row/column spacing
    # is untouched) -- without it, row 1's axes sit right at the top edge of the
    # page and the panel labels below get pushed past the figure boundary and clipped.
    # w_pad/h_pad add extra breathing room between the four subplots themselves.
    # right=0.995 (rather than past 1.0) leaves a sliver of margin so panel D's
    # "100" x-tick label, which sits flush against its axes' right edge, isn't cropped
    fig.tight_layout(pad=0.15, w_pad=3.0, h_pad=3.0, rect=[0, 0, 0.995, 0.965])

    # added after tight_layout, once each axes' final position is known, so all
    # four labels line up on the same absolute offset instead of an axes-fraction one
    for label, ax in zip(["A", "B", "C", "D"], axes.flat):
        add_panel_label(fig, ax, label)

    fig.savefig(f"{save_dir}/combined-figure.pdf")


def make_combined_n_and_s_percent_figure(df_default, data_dir_default, df_r1_high, data_dir_r1_high, save_dir):
    # spans both parameter sets side by side (default in column 0, high r1 rate
    # in column 1), so unlike make_combined_figure it can't be driven by a single
    # df/data_dir pair
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    colors = ["#2a78d6", "#eda100"]
    param_combos = [(10, 0.25), (1, 0.95)]
    value_labels = [f"n = {n}\ns = {s}" for n, s in param_combos]

    plot_trajectory_lines(axes[0, 0], df_default, data_dir_default, param_combos, colors)
    add_side_trajectory_legends(axes[0, 0], colors, value_labels, inside=True)

    plot_trajectory_lines(axes[0, 1], df_r1_high, data_dir_r1_high, param_combos, colors)
    add_side_trajectory_legends(axes[0, 1], colors, value_labels, inside=True)

    plot_percent_of_potential_heatmap(axes[1, 0], df_default, cbar_fraction=0.1)
    plot_percent_of_potential_heatmap(axes[1, 1], df_r1_high, cbar_fraction=0.1)

    # titled on both rows, not just the top -- a title on row 0 alone reads as
    # describing only that panel, not the whole column beneath it
    for ax in axes[:, 0]:
        ax.set_title("r1 rate = 0.01")
    for ax in axes[:, 1]:
        ax.set_title("r1 rate = 1/3")

    # top=0.985 (vs. the plain make_combined_figure's 0.965): tight_layout already
    # reserves room for the per-axes titles above row 0 as part of laying out the
    # axes themselves, so the extra top margin that make_combined_figure needs for
    # its panel labels (which have no title to make way for) is mostly redundant
    # here and was leaving a band of blank space above panels A/B
    fig.tight_layout(pad=0.15, w_pad=3.0, h_pad=3.0, rect=[0, 0, 0.995, 0.985])

    for label, ax in zip(["A", "B", "C", "D"], axes.flat):
        add_panel_label(fig, ax, label)

    fig.savefig(f"{save_dir}/combined-vary-n-and-s-percent-of-potential-figure.pdf")


# create figures for the default parameter set
df_default = pd.read_csv("data/s-vs-n-default.csv")
data_dir_default = "data/default"
save_dir_default = "figures/default"

make_max_load_heatmap(df_default, save_dir_default)
make_recovery_time_heatmap(df_default, save_dir_default)
make_vary_n_fixed_s_trajectories(df_default, data_dir_default, save_dir_default)
make_vary_s_fixed_n_trajectories(df_default, data_dir_default, save_dir_default)
make_vary_n_and_s_trajectories(df_default, data_dir_default, save_dir_default)
make_percent_of_potential_heatmap(df_default, save_dir_default)
make_combined_figure(df_default, data_dir_default, save_dir_default)


# repeat for high r1 rate
df_r1_high = pd.read_csv("data/s-vs-n-r1-high.csv")
data_dir_r1_high = "data/r1-high"
save_dir_r1_high = "figures/r1-high"

make_max_load_heatmap(df_r1_high, save_dir_r1_high)
make_recovery_time_heatmap(df_r1_high, save_dir_r1_high)
make_vary_n_fixed_s_trajectories(df_r1_high, data_dir_r1_high, save_dir_r1_high)
make_vary_s_fixed_n_trajectories(df_r1_high, data_dir_r1_high, save_dir_r1_high)
make_vary_n_and_s_trajectories(df_r1_high, data_dir_r1_high, save_dir_r1_high)
make_percent_of_potential_heatmap(df_r1_high, save_dir_r1_high)
make_combined_figure(df_r1_high, data_dir_r1_high, save_dir_r1_high)


# 4-panel figure spanning both parameter sets: vary-n-and-s trajectories (top)
# and percent-of-potential heatmaps (bottom), default vs high r1 rate
make_combined_n_and_s_percent_figure(df_default, data_dir_default, df_r1_high, data_dir_r1_high, "figures")