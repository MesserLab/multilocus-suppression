import glob
import os
import re
import string
import zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
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


MAX_LOAD_CBAR_LABEL = "Avg max\ngenetic load"
HALF_MAX_CBAR_LABEL = "Avg time to\nhalf max load"


def plot_max_load_heatmap(ax, df, cbar_pad=0.08, cbar_fraction=0.15, cbar=True):
    # cbar=False skips the per-axes colorbar -- for panels sharing one colorbar
    # that the caller draws from the returned mesh (the 0-1 range is already
    # fixed, so every panel maps load to color identically)
    pivot = df.pivot(index="target_fitness_cost", columns="num_targets", values="avg_max_load")
    n_vals = pivot.columns.values.astype(float)
    s_vals = pivot.index.values.astype(float)
    n_edges, s_edges = make_grid_edges(n_vals, s_vals)

    mesh = ax.pcolormesh(n_edges, s_edges, pivot.values, cmap="viridis", vmin=0, vmax=1, shading="flat")
    if cbar:
        cb = ax.figure.colorbar(mesh, ax=ax, shrink=0.7, pad=cbar_pad, fraction=cbar_fraction)
        cb.ax.set_title(MAX_LOAD_CBAR_LABEL, fontsize=10, pad=10)

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
    return mesh


def plot_half_max_time_heatmap(ax, df, cbar_pad=0.08, cbar_fraction=0.15, gamma=0.5, vmin=None, vmax=None,
        cbar=True):
    # vmin/vmax default to this df's own range; pass them explicitly (along with
    # cbar=False, as in plot_max_load_heatmap) for panels sharing one color scale
    pivot = df.pivot(index="target_fitness_cost", columns="num_targets",
        values="avg_time_to_below_half_max_genetic_load")
    n_vals = pivot.columns.values.astype(float)
    s_vals = pivot.index.values.astype(float)
    n_edges, s_edges = make_grid_edges(n_vals, s_vals)

    # -1 marks parameter combos that never fall below half max -- mask these out
    plot_values = np.where(pivot.values == -1, np.nan, pivot.values)

    # most cells fall well under 40 generations, but a handful run past 100 and
    # would otherwise compress everything else into near-black. rather than
    # changing the (linear, full-range) scale itself, resample magma so more of
    # its colors are spent at the low end: gamma < 1 packs colors toward small
    # values, gamma = 1 is plain magma
    cmap = ListedColormap(plt.get_cmap("magma")(np.linspace(0, 1, 256) ** gamma))
    cmap.set_bad("#DCDCDC")

    mesh = ax.pcolormesh(n_edges, s_edges, plot_values, cmap=cmap, vmin=vmin, vmax=vmax, shading="flat")
    if cbar:
        cb = ax.figure.colorbar(mesh, ax=ax, shrink=0.7, pad=cbar_pad, fraction=cbar_fraction)
        cb.ax.set_title(HALF_MAX_CBAR_LABEL, fontsize=10, pad=10)

    ax.set_xticks(n_vals)
    ax.set_yticks(s_vals[::2])
    ax.set_xlabel("Number of target sites (n)")
    ax.set_ylabel("Disrupted target site fitness cost (s)")
    ax.grid(False)
    return mesh


def find_replicate_csvs(data_dir, n, s):
    # replicate CSVs live inside zips whose folder names encode the parameters,
    # either as "-n10-s0.25-trajectories/" or as a parent/child pair like
    # "-n10-trajectories/s0.1/" -- pull n and s out of each member's folder path
    # (anchored on - or / so the "s-vs-n" prefix and "saturated" don't match)
    matches = []
    for zip_path in sorted(glob.glob(os.path.join(data_dir, "*.zip"))):
        with zipfile.ZipFile(zip_path) as zf:
            for member in zf.namelist():
                if not member.endswith(".csv") or "__MACOSX" in member:
                    continue
                folder = os.path.dirname(member)
                n_match = re.findall(r"(?:^|[-/])n(\d+)(?=[-/]|$)", folder)
                s_match = re.findall(r"(?:^|[-/])s(\d*\.?\d+)(?=[-/]|$)", folder)
                if n_match and s_match and int(n_match[-1]) == n and float(s_match[-1]) == s:
                    matches.append((zip_path, member))
    if not matches:
        raise FileNotFoundError(f"no replicate CSVs for n={n}, s={s} in {data_dir}")
    return matches


def load_avg_trajectory(data_dir, n, s):
    # each replicate stops once it reaches an absorbing state (drive lost/fixed),
    # so its final genetic load holds steady for any generation beyond that --
    # forward-fill each replicate out to the longest-running replicate before
    # averaging, rather than letting the average silently drop replicates
    reps = []
    for zip_path, member in find_replicate_csvs(data_dir, n, s):
        with zipfile.ZipFile(zip_path) as zf, zf.open(member) as f:
            rep_df = pd.read_csv(f, skipinitialspace=True)
        rep_df["gen_num"] = rep_df["gen_num"].astype(int)
        reps.append(rep_df.set_index("gen_num")["genetic_load"])

    min_gen = int(min(rep.index.min() for rep in reps))
    max_gen = int(max(rep.index.max() for rep in reps))
    full_index = np.arange(min_gen, max_gen + 1)
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


def make_combined_figure(df, data_dir, save_dir):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    colors = ["#2a78d6", "#1baf7a", "#eda100"]

    # a tighter colorbar (less pad/fraction reserved for it) leaves more of each
    # heatmap's own axes for the pcolormesh itself, widening A/B; this doesn't
    # touch C/D since they have no colorbar, and doesn't touch tight_layout's
    # own spacing since that's set by the pad= below, unchanged
    plot_max_load_heatmap(axes[0, 0], df, cbar_fraction=0.1)
    plot_half_max_time_heatmap(axes[0, 1], df, cbar_fraction=0.1)

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


def make_combined_n_and_s_figure(df_low, data_dir_low, df_high, data_dir_high, save_path):
    # spans two parameter sets side by side (r1 rate = 0.01 in column 0, high r1
    # rate in column 1), so unlike make_combined_figure it can't be driven by a
    # single df/data_dir pair
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = ["#2a78d6", "#eda100"]
    param_combos = [(10, 0.25), (1, 0.95)]
    value_labels = [f"n = {n}\ns = {s}" for n, s in param_combos]

    for ax, df, data_dir, title in zip(axes, [df_low, df_high], [data_dir_low, data_dir_high],
            ["r1 rate = 0.01", "r1 rate = 1/3"]):
        plot_trajectory_lines(ax, df, data_dir, param_combos, colors)
        add_side_trajectory_legends(ax, colors, value_labels, inside=True)
        ax.set_title(title)

    # left=0.02 leaves room for panel label A, which sits a fixed 0.5 in left of
    # its axes and otherwise gets clipped at the figure's left edge
    fig.tight_layout(pad=0.15, w_pad=3.0, rect=[0.02, 0, 0.995, 0.97])

    for label, ax in zip(["A", "B"], axes):
        add_panel_label(fig, ax, label)

    fig.savefig(save_path)


def make_recombination_heatmap_figure(recomb_dfs, save_path):
    # recomb_dfs: list of (row title, df) pairs, one row per recombination rate
    # (top to bottom). max load in the left column, time to half max in the right,
    # each column with one horizontal colorbar underneath
    fig, axes = plt.subplots(len(recomb_dfs), 2, figsize=(10, 5 * len(recomb_dfs) + 0.3), layout="constrained")
    # rect leaves a margin on the left and top for the panel labels, which sit a
    # fixed offset outside each axes and aren't accounted for by the layout itself
    fig.get_layout_engine().set(rect=(0.02, 0, 0.98, 0.99), h_pad=0.15, w_pad=0.15)

    # one colorbar per column shared by all rows, so a color means the same value
    # in every row -- time to half max needs a shared vmin/vmax for that, since
    # its range grows a lot as recombination drops (max load is already fixed 0-1)
    half_times = pd.concat([df["avg_time_to_below_half_max_genetic_load"] for _, df in recomb_dfs])
    half_times = half_times[half_times != -1]

    for row, (title, df) in enumerate(recomb_dfs):
        load_mesh = plot_max_load_heatmap(axes[row, 0], df, cbar=False)
        half_mesh = plot_half_max_time_heatmap(axes[row, 1], df, vmin=half_times.min(), vmax=half_times.max(),
            cbar=False)
        for ax in axes[row]:
            ax.set_title(title)
            # square plotting area regardless of figsize -- figsize above is just
            # tuned so the square panels fill the page without leftover whitespace
            ax.set_box_aspect(1)

    # pad is a fraction of the whole column's height (all rows), so the default
    # 0.05 leaves a wide gap under the bottom row -- keep it small
    for col, mesh, label in [(0, load_mesh, MAX_LOAD_CBAR_LABEL), (1, half_mesh, HALF_MAX_CBAR_LABEL)]:
        cb = fig.colorbar(mesh, ax=axes[:, col], location="bottom", shrink=0.8, fraction=0.025, aspect=40,
            pad=0.015)
        cb.set_label(label.replace("\n", " "))

    # constrained layout only settles axes positions at draw time, so draw once
    # before reading them for the panel labels
    fig.canvas.draw()
    for label, ax in zip(string.ascii_uppercase, axes.flat):
        add_panel_label(fig, ax, label)

    fig.savefig(save_path)


# summary csv for each parameter set -- the key doubles as both the folder under
# data/ holding its trajectory zips and the folder under figures/ to save into
datasets = {
    "default": "data/default/updated-s-vs-n-default.csv",
    "r1-high": "data/r1-high/updated-s-vs-n-r1-high.csv",
    "Cas9-default": "data/Cas9-default/updated-s-vs-n-Cas9-saturated.csv",
    "Cas9-r1-high": "data/Cas9-r1-high/updated-s-vs-n-Cas9-saturated-r1-high.csv",
}

dfs = {}
for name, csv_path in datasets.items():
    dfs[name] = pd.read_csv(csv_path)
    save_dir = f"figures/{name}"
    os.makedirs(save_dir, exist_ok=True)
    make_combined_figure(dfs[name], f"data/{name}", save_dir)


# 2-panel figures spanning both r1 rates: vary-n-and-s trajectories, low vs high r1 rate
make_combined_n_and_s_figure(dfs["default"], "data/default", dfs["r1-high"], "data/r1-high",
    "figures/combined-vary-n-and-s-trajectories.pdf")
make_combined_n_and_s_figure(dfs["Cas9-default"], "data/Cas9-default", dfs["Cas9-r1-high"], "data/Cas9-r1-high",
    "figures/combined-vary-n-and-s-trajectories-Cas9.pdf")


# heatmap figures across recombination rates, one row per rate from unlinked
# (top) down to the lowest rate (bottom)
recomb_rows = {
    "0.5": ("Recombination rate = 0.5", dfs["default"]),
    **{exp: (f"Recombination rate = $10^{{-{exp}}}$", pd.read_csv(f"data/recomb10_{exp}-s-vs-n.csv"))
        for exp in ["4", "5", "7"]},
}

# 3x2: unlinked, a middle rate, and the lowest rate -- compare 10^-4 vs 10^-5;
# each gets its own output file
for middle_rate in ["4", "5"]:
    make_recombination_heatmap_figure([recomb_rows[r] for r in ["0.5", middle_rate, "7"]],
        f"figures/recombination-heatmaps-1e-{middle_rate}.pdf")

# 4x2: all four rates
make_recombination_heatmap_figure([recomb_rows[r] for r in ["0.5", "4", "5", "7"]],
    "figures/recombination-heatmaps-all.pdf")