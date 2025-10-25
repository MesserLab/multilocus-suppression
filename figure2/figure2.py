# -*- coding: utf-8 -*-
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# change to the directory of the data
path = ("/SSD/reg259/gene_drive_sims/20251021_figure2A-B/")
os.chdir(path)

# list out all the files in the directory
files = sorted([i for i in os.listdir(path) if "output" in i or "parameters" in i])
output, parameters = files[::2], files[1::2]

# collect all the data from the files
all_data = []
sim_id = 0

for out, par in zip(output, parameters):
    # collect parameters
    with open(par, "r") as text:
        par_dict = {}
        for line in text.readlines():
            key, value = line.replace("\n", "").split(": ")
            par_dict[key] = float(value)
    
    # load output data
    data = pd.read_csv(out, sep=", ", engine="python")
    data["r1"] = par_dict["resistance_rate"]
    data["sim_id"] = sim_id
    sim_id += 1
    all_data.append(data[["gen_num", "mean_fit", "wt_rate", "drive_rate", "resist_rate", "r1", "sim_id"]])

# combine all into one DataFrame
df = pd.concat(all_data, ignore_index=True)

# figure out which simulations end in resistance
pivot_df = df.pivot_table(index="sim_id", values="resist_rate", aggfunc="last")
pivot_df["outcome_resist"] = pivot_df["resist_rate"] > 0.5
df["outcome_resist"] = df["sim_id"].map(pivot_df["outcome_resist"])

# set up the plot
fig, axes = plt.subplots(2, 2, figsize=(15, 11), layout="constrained")
fig.set_constrained_layout_pads(hspace=0.05)
plt.rcParams["font.size"] = 13

####################################################
# color palette for R1 values
r1_values = sorted(df["r1"].unique())
cmap = plt.colormaps["rocket_r"]
uniform_colors = [cmap(0.05 + 0.9 * (i / (len(r1_values) - 1))) for i in range(len(r1_values))]
r1_color_map = {r:color for r,color in zip(r1_values, uniform_colors)}

# plot the first panel
for r1 in r1_values:
    # plot simulations that end in resistance
    subset = df[(df["r1"] == r1) & (df["outcome_resist"] == True)]
    num_resist = subset["sim_id"].nunique()

    # plot average line
    avg_data = subset.groupby("gen_num")["mean_fit"].mean().reset_index()
    axes[0,0].plot(avg_data["gen_num"], avg_data["mean_fit"], color=r1_color_map[r1], linewidth=2, label=f"{r1:.0e}")

    # plot simulations that end in fixation
    subset = df[(df["r1"] == r1) & (df["outcome_resist"] == False)]
    num_fix = subset["sim_id"].nunique()

    # plot average line
    avg_data = subset.groupby("gen_num")["mean_fit"].mean().reset_index()
    axes[0,0].plot(avg_data["gen_num"], avg_data["mean_fit"], color=r1_color_map[r1], linewidth=2, linestyle="--", label=f"{r1:.0e}")

    # add the fraction text
    if not avg_data.empty:
        last_gen = avg_data["gen_num"].max()
        last_val = avg_data["mean_fit"].iloc[-1]
        offset = {1e-06:-0.025, 1e-05:0.025}
        axes[0,0].text(last_gen + 0.5, last_val + offset[r1], f"{num_fix}/{num_fix + num_resist}", color=r1_color_map[r1], fontsize=11, va="center")

# axis formatting
axes[0,0].set_title("dsx drive, s=1, h=0", fontsize=15, loc="center")
axes[0,0].set_xlabel("Generations Since Drive Release", fontsize=13)
axes[0,0].set_ylabel("Mean Population Fitness", fontsize=13)
axes[0,0].set_xlim(0, 50)
axes[0,0].set_ylim(-0.05, 1.05)
axes[0,0].tick_params(axis="both", labelsize=13)
axes[0,0].xaxis.set_major_locator(ticker.MultipleLocator(5))
axes[0,0].yaxis.set_major_locator(ticker.MultipleLocator(0.1))

# add both legends
legend1 = [Line2D([0], [0], color=r1_color_map[r], lw=2, label=f"{r:.0e}") for r in r1_values]
legend2 = [Line2D([0], [0], color="gray", lw=2, linestyle="-", label="resist"),
           Line2D([0], [0], color="gray", lw=2, linestyle="--", label="fix")]

axes[0,0].add_artist(axes[0,0].legend(handles=legend1, title="r1 rate", loc="lower right"))
axes[0,0].legend(handles=legend2, loc="center right", bbox_to_anchor=(1, 0.47))

####################################################
# colormaps for each allele
gray_color_map = {r: plt.colormaps["Greys"](0.3 + 0.6 * (i / (len(r1_values) - 1))) for i, r in enumerate(r1_values)}
wt_color_map = {r: plt.colormaps["Blues"](0.3 + 0.6 * (i / (len(r1_values) - 1))) for i, r in enumerate(r1_values)}
drive_color_map = {r: plt.colormaps["Reds"](0.3 + 0.6 * (i / (len(r1_values) - 1))) for i, r in enumerate(r1_values)}
resist_color_map = {r: plt.colormaps["Greens"](0.3 + 0.6 * (i / (len(r1_values) - 1))) for i, r in enumerate(r1_values)}

# plot the second panel
for r1 in r1_values:
    # plot simulations that end in resistance
    subset = df[(df["r1"] == r1) & (df["outcome_resist"] == True)]
    num_resist = subset["sim_id"].nunique()

    # plot mean allele frequency trajectories
    avg_data = subset.groupby("gen_num")[["wt_rate", "drive_rate", "resist_rate"]].mean().reset_index()
    axes[0,1].plot(avg_data["gen_num"], avg_data["wt_rate"], color=wt_color_map[r1], linewidth=1.5)
    axes[0,1].plot(avg_data["gen_num"], avg_data["drive_rate"], color=drive_color_map[r1], linewidth=1.5)
    axes[0,1].plot(avg_data["gen_num"], avg_data["resist_rate"], color=resist_color_map[r1], linewidth=1.5)

    # plot simulations that end in fixation
    subset = df[(df["r1"] == r1) & (df["outcome_resist"] == False)]
    num_fix = subset["sim_id"].nunique()

    # plot mean allele frequency trajectories
    avg_data = subset.groupby("gen_num")[["wt_rate", "drive_rate", "resist_rate"]].mean().reset_index()
    axes[0,1].plot(avg_data["gen_num"], avg_data["wt_rate"], color=wt_color_map[r1], linewidth=1.5, linestyle="--")
    axes[0,1].plot(avg_data["gen_num"], avg_data["drive_rate"], color=drive_color_map[r1], linewidth=1.5, linestyle="--")
    axes[0,1].plot(avg_data["gen_num"], avg_data["resist_rate"], color=resist_color_map[r1], linewidth=1.5, linestyle="--")

    # add the fraction text
    if not avg_data.empty:
        last_gen = avg_data["gen_num"].max()
        last_val = avg_data["drive_rate"].iloc[-1]
        offset = {1e-06:-0.025, 1e-05:0.025}
        axes[0,1].text(last_gen + 0.5, last_val + offset[r1], f"{num_fix}/{num_fix + num_resist}", color=drive_color_map[r1], fontsize=11, va="center")

# axis formatting
axes[0,1].set_title("dsx drive, s=1, h=0", fontsize=15, loc="center")
axes[0,1].set_xlabel("Generations Since Drive Release", fontsize=13)
axes[0,1].set_ylabel("Allele Frequency", fontsize=13)
axes[0,1].set_xlim(0, 50)
axes[0,1].set_ylim(-0.05, 1.05)
axes[0,1].tick_params(axis="both", labelsize=13)
axes[0,1].xaxis.set_major_locator(ticker.MultipleLocator(5))
axes[0,1].yaxis.set_major_locator(ticker.MultipleLocator(0.1))

# allele type legend
allele_legend_handles = [Line2D([0], [0], color="tab:blue", lw=2, label="wild-type"),
                         Line2D([0], [0], color="tab:red", lw=2, label="drive"),
                         Line2D([0], [0], color="tab:green", lw=2, label="resistance")]
axes[0,1].add_artist(axes[0,1].legend(handles=allele_legend_handles, loc="lower right", bbox_to_anchor=(1, 0.12)))

# resist or fix legend
resist_fix = [Line2D([0], [0], color="gray", lw=2, linestyle="-", label="resist"),
              Line2D([0], [0], color="gray", lw=2, linestyle="--", label="fix")]
axes[0,1].add_artist(axes[0,1].legend(handles=resist_fix, loc="upper right", bbox_to_anchor=(1, 0.88)))

# r1 rate legend
r1_legend_handles = [Line2D([0], [0], color=gray_color_map[r], lw=2, label=f"{r:.0e}") for r in r1_values]
axes[0,1].legend(handles=r1_legend_handles, title="r1 rate", loc="center right", bbox_to_anchor=(1, 0.53))

####################################################
# change to the directory of the data
path = ("/SSD/reg259/gene_drive_sims/20251021_figure2C-D/")
os.chdir(path)

# list out all the files in the directory
files = sorted([i for i in os.listdir(path) if "output" in i or "parameters" in i])
output, parameters = files[::2], files[1::2]

# collect all the data from the files
all_data = []
sim_id = 0

for out, par in zip(output, parameters):
    # collect parameters
    with open(par, "r") as text:
        par_dict = {}
        for line in text.readlines():
            key, value = line.replace("\n", "").split(": ")
            par_dict[key] = float(value)
    
    # load output data
    data = pd.read_csv(out, sep=", ", engine="python")
    data["s"] = par_dict["drive_coeff"]
    data["sim_id"] = sim_id
    sim_id += 1
    all_data.append(data[["gen_num", "mean_fit", "wt_rate", "drive_rate", "resist_rate", "s", "sim_id"]])

# combine all into one DataFrame
df = pd.concat(all_data, ignore_index=True)

# figure out which simulations end in resistance
pivot_df = df.pivot_table(index="sim_id", values="resist_rate", aggfunc="last")
pivot_df["outcome_resist"] = pivot_df["resist_rate"] > 0.5
df["outcome_resist"] = df["sim_id"].map(pivot_df["outcome_resist"])

####################################################
# color palette for s values
s_values = sorted(df["s"].unique(), reverse=True)
cmap = plt.colormaps["rocket_r"]
uniform_colors = [cmap(0.05 + 0.9 * (i / (len(s_values) - 1))) for i in range(len(s_values))]
s_color_map = {s:color for s,color in zip(s_values, uniform_colors)}

# plot the third panel
for s_val in s_values:
    # plot simulations that end in resistance
    subset = df[(df["s"] == s_val) & (df["outcome_resist"] == True)]
    num_resist = subset["sim_id"].nunique()

    # plot average line
    avg_data = subset.groupby("gen_num")["mean_fit"].mean().reset_index()
    axes[1,0].plot(avg_data["gen_num"], avg_data["mean_fit"], color=s_color_map[s_val], linewidth=2, label=s_val)

    # plot simulations that end in fixation
    subset = df[(df["s"] == s_val) & (df["outcome_resist"] == False)]
    num_fix = subset["sim_id"].nunique()

    # plot average line
    avg_data = subset.groupby("gen_num")["mean_fit"].mean().reset_index()
    axes[1,0].plot(avg_data["gen_num"], avg_data["mean_fit"], color=s_color_map[s_val], linewidth=2, linestyle="--", label=s_val)

    # add the fraction text
    if not avg_data.empty:
        last_gen = avg_data["gen_num"].max()
        last_val = avg_data["mean_fit"].iloc[-1]
        axes[1,0].text(last_gen + 2, last_val + 0.025, f"{num_fix}/{num_fix + num_resist}", color=s_color_map[s_val], fontsize=11, va="center")

# axis formatting
axes[1,0].set_title("homing drive, h=0.5, r1=0.001", fontsize=15, loc="center")
axes[1,0].set_xlabel("Generations Since Drive Release", fontsize=13)
axes[1,0].set_ylabel("Mean Population Fitness", fontsize=13)
axes[1,0].set_xlim(0, 250)
axes[1,0].set_ylim(-0.05, 1.05)
axes[1,0].tick_params(axis="both", labelsize=13)
axes[1,0].xaxis.set_major_locator(ticker.MultipleLocator(25))
axes[1,0].yaxis.set_major_locator(ticker.MultipleLocator(0.1))

# add both legends
legend1 = [Line2D([0], [0], color=s_color_map[s], lw=2, label=s) for s in s_values]
legend2 = [Line2D([0], [0], color="gray", lw=2, linestyle="-", label="resist"),
           Line2D([0], [0], color="gray", lw=2, linestyle="--", label="fix")]

axes[1,0].add_artist(axes[1,0].legend(handles=legend1, title="s", loc="lower right"))
axes[1,0].legend(handles=legend2, loc="center right", bbox_to_anchor=(1, 0.71))

####################################################
# colormaps for each allele
gray_color_map = {r: plt.colormaps["Greys"](0.3 + 0.6 * (i / (len(s_values) - 1))) for i, r in enumerate(s_values)}
wt_color_map = {r: plt.colormaps["Blues"](0.3 + 0.6 * (i / (len(s_values) - 1))) for i, r in enumerate(s_values)}
drive_color_map = {r: plt.colormaps["Reds"](0.3 + 0.6 * (i / (len(s_values) - 1))) for i, r in enumerate(s_values)}
resist_color_map = {r: plt.colormaps["Greens"](0.3 + 0.6 * (i / (len(s_values) - 1))) for i, r in enumerate(s_values)}

# plot the fourth panel
for s_val in s_values:
    # plot simulations that end in resistance
    subset = df[(df["s"] == s_val) & (df["outcome_resist"] == True)]
    num_resist = subset["sim_id"].nunique()

    # plot mean allele frequency trajectories
    avg_data = subset.groupby("gen_num")[["wt_rate", "drive_rate", "resist_rate"]].mean().reset_index()
    axes[1,1].plot(avg_data["gen_num"], avg_data["wt_rate"], color=wt_color_map[s_val], linewidth=1.5)
    axes[1,1].plot(avg_data["gen_num"], avg_data["drive_rate"], color=drive_color_map[s_val], linewidth=1.5)
    axes[1,1].plot(avg_data["gen_num"], avg_data["resist_rate"], color=resist_color_map[s_val], linewidth=1.5)

    # plot simulations that end in fixation
    subset = df[(df["s"] == s_val) & (df["outcome_resist"] == False)]
    num_fix = subset["sim_id"].nunique()

    # plot mean allele frequency trajectories
    avg_data = subset.groupby("gen_num")[["wt_rate", "drive_rate", "resist_rate"]].mean().reset_index()
    axes[1,1].plot(avg_data["gen_num"], avg_data["wt_rate"], color=wt_color_map[s_val], linewidth=1.5, linestyle="--")
    axes[1,1].plot(avg_data["gen_num"], avg_data["drive_rate"], color=drive_color_map[s_val], linewidth=1.5, linestyle="--")
    axes[1,1].plot(avg_data["gen_num"], avg_data["resist_rate"], color=resist_color_map[s_val], linewidth=1.5, linestyle="--")

    # add the fraction text
    if not avg_data.empty:
        last_gen = avg_data["gen_num"].max()
        last_val = avg_data["drive_rate"].iloc[-1]
        axes[1,1].text(last_gen + 2, last_val + 0.025, f"{num_fix}/{num_fix + num_resist}", color=drive_color_map[s_val], fontsize=11, va="center")

# axis formatting
axes[1,1].set_title("homing drive, h=0.5, r1=0.001", fontsize=15, loc="center")
axes[1,1].set_xlabel("Generations Since Drive Release", fontsize=13)
axes[1,1].set_ylabel("Allele Frequency", fontsize=13)
axes[1,1].set_xlim(0, 250)
axes[1,1].set_ylim(-0.05, 1.05)
axes[1,1].tick_params(axis="both", labelsize=13)
axes[1,1].xaxis.set_major_locator(ticker.MultipleLocator(25))
axes[1,1].yaxis.set_major_locator(ticker.MultipleLocator(0.1))

# allele type legend
allele_legend_handles = [Line2D([0], [0], color="tab:blue", lw=2, label="wild-type"),
                         Line2D([0], [0], color="tab:red", lw=2, label="drive"),
                         Line2D([0], [0], color="tab:green", lw=2, label="resistance")]
axes[1,1].add_artist(axes[1,1].legend(handles=allele_legend_handles, loc="lower right", bbox_to_anchor=(1, 0), framealpha=1))

# resist or fix legend
resist_fix = [Line2D([0], [0], color="gray", lw=2, linestyle="-", label="resist"),
              Line2D([0], [0], color="gray", lw=2, linestyle="--", label="fix")]
axes[1,1].add_artist(axes[1,1].legend(handles=resist_fix, loc="upper right", bbox_to_anchor=(1, 1), framealpha=1))

# s value legend
s_legend_handles = [Line2D([0], [0], color=gray_color_map[s], lw=2, label=s) for s in s_values]
axes[1,1].legend(handles=s_legend_handles, title="s", loc="center right", bbox_to_anchor=(1, 0.53))

# save figure
plt.savefig("../figure2.png", dpi=450, bbox_inches="tight")
plt.close()