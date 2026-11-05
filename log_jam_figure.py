"""Log-jam haplotype figure: haplotype frequency vs generation, with haplotype
cartoons at time points between the arrivals of persistent haplotypes. Reads
the per-generation HAPLOTYPES blocks that SLiM prints when
PRINT_COMMON_HAPLOTYPES=T (see log_jam_sims.sh).

Usage: python log_jam_figure.py N   (data/log-jam-sims/log-jam-simN.txt ->
                                     figures/log-jam-sims/log-jam-simN.pdf)"""

import os
import re
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

plt.rcParams.update({
    "font.size": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})

# set up the argument parser
parser = argparse.ArgumentParser()
parser.add_argument("hap_run", type=int, help="The haplotype run number (integer)")
args = parser.parse_args()

HAP_FILE = f"data/log-jam-sims/log-jam-sim{args.hap_run}.txt"
SAVE_PATH = f"figures/log-jam-sims/log-jam-sim{args.hap_run}.pdf"

# True: haplotypes include the drive site (d_rddr and r_rddr are different
#       haplotypes, and the cartoons draw the drive site as well)
# False: frequencies are summed over the drive site, so only the target sites
#        are shown (d_rddr + r_rddr -> rddr)
SHOW_DRIVE_SITE = False

# Target-site states: r = resistance, d = disrupted, w = wild-type. SLiM prints
# disrupted sites as "b"; parse_haplotypes converts them to "d".
# Drive-site states: d = drive, r = resistance, w = wild-type.
# Haplotype colours are assigned in order of first appearance.
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9",
           "#882255", "#999933", "#332288", "#44AA99"]
TARGET_BLUE = "#2b8fd1"
RESIST_GREEN = "#4c9a3c"
DISRUPT_RED = "#e8302a"
DRIVE_RED = "#e8302a"
DRIVE_WT_BLACK = "#111111"

# a haplotype is persistent if its frequency stays >= MIN_FREQ for at least
# MIN_GENS generations in a row; it arises at the start of that stretch.
# Arrivals closer together than MIN_GENS count as one (the founders and early
# recombinants all appear while the drive spreads). Time points go midway
# between arrivals; their cartoons show the persistent haplotypes still there.
MIN_FREQ = 0.05
MIN_GENS = 200
# boxes are sized as if there were at least this many cartoons in the row
MIN_CARTOON_COLS = 4

# figure layout, in inches. The cartoon boxes get taller with the largest
# number of haplotypes shown in any cartoon, and the figure height follows.
FIG_WIDTH_IN = 12.0
PLOT_HEIGHT_IN = 3.9         # frequency plot
GAP_IN = 0.68                # between the cartoon row and the frequency plot
MARGIN_TOP_IN = 0.2
MARGIN_BOTTOM_IN = 0.7       # room for the x-axis label
BOX_TITLE_IN = 0.30          # box top -> centre of the "Generation N" title
BOX_FIRST_ROW_IN = 0.65      # box top -> first haplotype row
BOX_ROW_IN = 0.35            # spacing between haplotype rows
BOX_BOTTOM_IN = 0.31         # last haplotype row -> box bottom
KEY_LINE_IN = 0.27           # line spacing in the site key
# the site key is top-aligned with the cartoon boxes and right-aligned with
# the line-plot legend


def parse_haplotypes(path, show_drive=SHOW_DRIVE_SITE):
    """Return a DataFrame (index = generation, columns = haplotype) of
    haplotype frequencies. Haplotypes are "<drive>_<targets>" (e.g. "d_rddr")
    if show_drive, otherwise just "<targets>" summed over the drive site.
    Only the top NUM_HAPLOTYPES_PRINT haplotypes are printed each generation.
    A haplotype missing from a generation's list is at frequency 0 only if
    the list wasn't full; if it was full, the haplotype may just have ranked
    lower, so its frequency there is unknown and is interpolated from the
    generations around it (a haplotype never exits because of the cutoff)."""
    rows = {}
    n_printed = {}
    gen = None
    hap_re = re.compile(r"^([a-z])_([a-z]+):\s+(\d+)\s+([\d.eE+-]+)")
    with open(path) as fh:
        for line in fh:
            if line.startswith("GEN::"):
                gen = int(line.split()[1])
                rows[gen] = {}
                n_printed[gen] = 0
                continue
            m = hap_re.match(line)
            if m and gen is not None:
                n_printed[gen] += 1
                sites = m.group(2).replace("b", "d")
                key = f"{m.group(1)}_{sites}" if show_drive else sites
                rows[gen][key] = rows[gen].get(key, 0.0) + float(m.group(4))
    df = pd.DataFrame.from_dict(rows, orient="index").sort_index()
    full_range = range(df.index.min(), df.index.max() + 1)
    df = df.reindex(full_range)
    n_printed = pd.Series(n_printed).reindex(full_range).fillna(0)
    full = n_printed >= n_printed.max()    # the list hit the print limit
    # missing from a list that wasn't full: truly 0; from a full list: unknown
    df = df.where(df.notna() | full.values[:, None], 0.0)
    return df.interpolate(limit_area="inside").fillna(0.0)


def target_sites(hap):
    """Target-site part of a haplotype name ("d_rddr" -> "rddr")."""
    return hap.split("_")[-1]


def drive_site(hap):
    """Drive-site state of a haplotype name ("d_rddr" -> "d"), or None."""
    return hap.split("_")[0] if "_" in hap else None


def stretches(mask):
    """(start, stop) index pairs of the runs of True in a boolean array."""
    edges = np.flatnonzero(np.diff(np.r_[0, mask.astype(int), 0]))
    return list(zip(edges[::2], edges[1::2]))


def persistent_episodes(df, min_freq=MIN_FREQ, min_gens=MIN_GENS):
    """Split each haplotype's history into episodes (unbroken stretches of
    frequency > 0; it exits when it drops to 0). An episode is persistent if
    the frequency stays >= `min_freq` for `min_gens` generations in a row
    within it. A haplotype that exits and re-enters is judged afresh.

    Returns (keep, arrivals): keep is a boolean DataFrame like df, True inside
    persistent episodes; arrivals is a list of (gen, haplotype), one per
    persistent episode, at the start of its first stretch >= `min_freq`."""
    gens = df.index.values
    keep = pd.DataFrame(False, index=df.index, columns=df.columns)
    arrivals = []
    for hap in df.columns:
        freq = df[hap].values
        for a, b in stretches(freq > 0):
            high = [(s, e) for s, e in stretches(freq[a:b] >= min_freq)
                    if e - s >= min_gens]
            if high:
                keep.iloc[a:b, keep.columns.get_loc(hap)] = True
                arrivals.append((gens[a + high[0][0]], hap))
    return keep, sorted(arrivals)


def detect_events(df, min_freq=MIN_FREQ, min_gens=MIN_GENS):
    """Find the time points to depict.

    Arrivals of persistent episodes (see persistent_episodes) closer together
    than `min_gens` are merged. A time point goes midway between one arrival
    and the next, and midway between the last arrival and the end of the run
    (or at the end, if a haplotype is fixed there). Its cartoon shows every
    haplotype in a persistent episode at that generation.

    Returns (events, colors, keep): events is a list of (gen, haplotypes);
    colors maps each haplotype in a cartoon to its colour; keep is from
    persistent_episodes() (where to draw the lines).
    """
    keep, arrivals = persistent_episodes(df, min_freq, min_gens)
    if not arrivals:
        raise RuntimeError("no persistent haplotypes; try lowering MIN_FREQ "
                           "or MIN_GENS")
    first_arrival = {}
    for t, h in arrivals:
        first_arrival.setdefault(h, t)
    times = [t for t, _ in arrivals]
    merged = [times[0]]                    # last arrival of each merged group
    for t in times[1:]:
        if t - merged[-1] < min_gens:
            merged[-1] = t
        else:
            merged.append(t)
    bounds = merged + [df.index[-1]]

    def segregating(gen):
        return sorted(keep.columns[keep.loc[gen].values])

    events = [(int((a + b) // 2), segregating((a + b) // 2))
              for a, b in zip(bounds[:-1], bounds[1:])]
    # if the run ends with a haplotype fixed, the last time point is the
    # final generation (showing just the fixed haplotype)
    last = df.index[-1]
    if df.loc[last].max() >= 1.0 - 1e-9:
        events[-1] = (int(last), segregating(last))

    def sort_key(h):
        t = target_sites(h)
        return (t.index("r") if "r" in t else len(t), t, drive_site(h) or "")

    # founders (first time point) by position of their r, then the other
    # cartoon haplotypes in the order they arise
    order = sorted(events[0][1], key=sort_key)
    later = {h for _, haps in events for h in haps} - set(order)
    order += sorted(later, key=lambda h: (first_arrival[h], h))
    colors = {h: PALETTE[i % len(PALETTE)] for i, h in enumerate(order)}
    return events, colors, keep


def draw_site(ax, x0, x1, y, state, lw=3, ms=9):
    """One target site in the style of the model schematic: a blue segment,
    with a green dot (resistance) or red X (disrupted) on top."""
    ax.plot([x0, x1], [y, y], color=TARGET_BLUE, lw=lw, solid_capstyle="butt",
            zorder=2)
    xm = (x0 + x1) / 2
    if state == "r":
        ax.plot(xm, y, "o", color=RESIST_GREEN, ms=ms, zorder=3)
    elif state == "d":
        ax.plot(xm, y, "x", color=DISRUPT_RED, ms=ms, mew=2.5, zorder=3)


def draw_drive_site(ax, x0, x1, y, state, lw=3, ms=9):
    """Drive site in the style of the model schematic: red line for the drive,
    black line for wild-type, black line + green dot for resistance."""
    color = DRIVE_RED if state == "d" else DRIVE_WT_BLACK
    ax.plot([x0, x1], [y, y], color=color, lw=lw, solid_capstyle="butt",
            zorder=2)
    if state == "r":
        ax.plot((x0 + x1) / 2, y, "o", color=RESIST_GREEN, ms=ms, zorder=3)


def haplotype_layout(hap):
    """Segment width, gap and extra gap after the drive site; segments are
    shrunk when the drive site is drawn so the row still fits in the box."""
    if drive_site(hap) is None:
        return 0.11, 0.02, 0.0
    return 0.10, 0.02, 0.04


def haplotype_width(hap):
    site_w, gap, drive_gap = haplotype_layout(hap)
    n = len(target_sites(hap)) + (drive_site(hap) is not None)
    return n * site_w + (n - 1) * gap + drive_gap


def box_height_in(n_rows):
    """Height (inches) of a cartoon box holding n_rows haplotypes."""
    return BOX_FIRST_ROW_IN + (n_rows - 1) * BOX_ROW_IN + BOX_BOTTOM_IN


def draw_haplotype(ax, x, y, hap, color, tab_h=0.16):
    """Coloured identity tab (matches the line plot) + the sites (drive site
    first, if the haplotype has one). The target sites are linked, so they
    sit on one continuous blue line."""
    site_w, gap, drive_gap = haplotype_layout(hap)
    ax.add_patch(FancyBboxPatch((x - 0.1, y - tab_h / 2), 0.05, tab_h,
                                boxstyle="round,pad=0,rounding_size=0.015",
                                fc=color, ec="none", zorder=2))
    if drive_site(hap) is not None:
        draw_drive_site(ax, x, x + site_w, y, drive_site(hap))
        x += site_w + gap + drive_gap
    n = len(target_sites(hap))
    ax.plot([x, x + n * site_w + (n - 1) * gap], [y, y], color=TARGET_BLUE,
            lw=3, solid_capstyle="butt", zorder=2)
    for i, s in enumerate(target_sites(hap)):
        x0 = x + i * (site_w + gap)
        draw_site(ax, x0, x0 + site_w, y, s)


def draw_cartoon(ax, gen, haps, colors, height_in):
    """x runs 0-1 across the box; y is in inches from the box bottom."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, height_in)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.96,
                                boxstyle="round,pad=0,rounding_size=0.04",
                                fc="#f7f7f5", ec="#cfcfcf", lw=1.0,
                                transform=ax.transAxes))
    ax.text(0.5, height_in - BOX_TITLE_IN, f"Generation {gen}", fontsize=13,
            fontweight="bold", ha="center", va="center")
    x = 0.5 - haplotype_width(haps[0]) / 2 + 0.025
    for k, hap in enumerate([h for h in colors if h in haps]):
        y = height_in - BOX_FIRST_ROW_IN - k * BOX_ROW_IN
        draw_haplotype(ax, x, y, hap, colors[hap])


def make_haplotype_figure(hap_file=HAP_FILE, save_path=SAVE_PATH,
                          show_drive=SHOW_DRIVE_SITE):
    df = parse_haplotypes(hap_file, show_drive=show_drive)
    events, colors, keep = detect_events(df)
    shown = df.where(keep)                 # NaN outside persistent episodes
    for gen, haps in events:
        print(f"time point at gen {gen}: {len(haps)} haplotypes ({', '.join(haps)})")

    # site key entries: (draw function or None for a header, state, label)
    key = []
    if show_drive:
        drive_states = sorted({drive_site(h) for h in colors},
                              key="dwr".index)
        drive_labels = {"d": "drive (d)", "r": "resistance (r)",
                        "w": "wild-type (w)"}
        key += [(None, None, "Drive site")]
        key += [(draw_drive_site, s, drive_labels[s]) for s in drive_states]
    key += [(None, None, "Target site"),
            (draw_site, "r", "resistance (r)"),
            (draw_site, "d", "disrupted (d)")]
    key_h_in = len(key) * KEY_LINE_IN

    # vertical layout in inches; boxes fit the largest cartoon, and the gap
    # grows if the site key would otherwise hang down into the plot legend
    box_h_in = box_height_in(max(len(haps) for _, haps in events))
    gap_in = max(GAP_IN, key_h_in - box_h_in + 0.2)
    fig_h_in = (MARGIN_TOP_IN + box_h_in + gap_in + PLOT_HEIGHT_IN
                + MARGIN_BOTTOM_IN)
    fig = plt.figure(figsize=(FIG_WIDTH_IN, fig_h_in))
    left, right = 0.07, 0.9            # figure fractions (width is fixed)
    plot_y0 = MARGIN_BOTTOM_IN / fig_h_in
    box_y0 = (MARGIN_BOTTOM_IN + PLOT_HEIGHT_IN + gap_in) / fig_h_in
    box_top = box_y0 + box_h_in / fig_h_in
    ax = fig.add_axes([left, plot_y0, right - left, PLOT_HEIGHT_IN / fig_h_in])

    # cartoon boxes are sized for at least MIN_CARTOON_COLS per row; with
    # fewer steady states they keep that width and are centred
    n_cols = max(len(events), MIN_CARTOON_COLS)
    wspace = 0.12                      # gap between boxes / box width
    box_w = (right - left) / (n_cols + wspace * (n_cols - 1))
    pitch = box_w * (1 + wspace)
    x_start = left + (n_cols - len(events)) / 2 * pitch
    # each box is sized to its own cartoon and hangs from a shared top edge;
    # the tallest box sets the space reserved above the plot
    for i, (gen, haps) in enumerate(events):
        h_in = box_height_in(len(haps))
        cax = fig.add_axes([x_start + i * pitch, box_top - h_in / fig_h_in,
                            box_w, h_in / fig_h_in])
        draw_cartoon(cax, gen, haps, colors, h_in)

    # only haplotypes that appear in a cartoon are plotted, and only during
    # their persistent episodes (not before entering, after exiting, or when
    # they re-enter briefly)
    gens = df.index.values
    for hap, color in colors.items():
        ax.plot(gens, shown[hap].values, color=color, lw=2, label=hap,
                zorder=2)

    # generation labels; nudge apart neighbours closer than 6% of the x-range
    event_gens = [gen for gen, _ in events]
    close = 0.06 * gens[-1]
    for i, gen in enumerate(event_gens):
        ha = "center"
        if i + 1 < len(event_gens) and event_gens[i + 1] - gen < close:
            ha = "right"
        elif i > 0 and gen - event_gens[i - 1] < close:
            ha = "left"
        ax.axvline(gen, color="#9a9a9a", ls="--", lw=1, zorder=0,
                   clip_on=False)      # stays whole at the right edge
        ax.text(gen, 1.02, str(gen), transform=ax.get_xaxis_transform(),
                ha=ha, va="bottom", fontsize=12, fontweight="bold")

    ax.set_xlim(0, gens[-1])
    # headroom above the highest line; up to 1.03 so lines at fixation
    # aren't clipped while the top tick stays at 1.0
    ax.set_ylim(0, min(1.03, 1.08 * np.nanmax(shown[list(colors)].values)))
    ax.set_xlabel("Generation", fontsize=14)
    ax.set_ylabel("Haplotype frequency", fontsize=14)
    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.06, 1.0),
                    frameon=False, prop={"family": "monospace", "size": 11},
                    handlelength=1.5)

    # site key, drawn with the same marks as the cartoons; its top is level
    # with the box tops (offset by the box outline inset)
    key_w = 0.1                        # figure fraction; fixes the mark size
    inset = 0.02 * box_h_in / fig_h_in
    kax = fig.add_axes([0.905, box_top - inset - key_h_in / fig_h_in, key_w,
                        key_h_in / fig_h_in])
    kax.set_xlim(0, 1)
    kax.set_ylim(0, len(key))
    kax.axis("off")
    key_texts = []
    for k, (draw, state, label) in enumerate(key):
        y = len(key) - k - 0.5
        if draw is None:
            key_texts.append(kax.text(0.0, y, label, fontsize=11,
                                      fontweight="bold", va="center"))
        else:
            draw(kax, 0.0, 0.22, y, state)
            key_texts.append(kax.text(0.3, y, label, fontsize=10,
                                      va="center"))

    # shift the key so its right edge lines up with the legend's right edge
    fig.canvas.draw()
    to_fig = fig.transFigure.inverted()
    leg_right = max(to_fig.transform(t.get_window_extent())[1, 0]
                    for t in leg.get_texts())
    key_right = max(to_fig.transform(t.get_window_extent())[1, 0]
                    for t in key_texts)
    pos = kax.get_position()
    kax.set_position([pos.x0 + leg_right - key_right, pos.y0, pos.width,
                      pos.height])

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"saved {save_path}")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    make_haplotype_figure()
