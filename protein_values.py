import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# --------------------------------------------------
# 1. Load data
# --------------------------------------------------

df = pd.read_excel("Proteome - maternal BE.xlsx")


# --------------------------------------------------
# 2. Proteins to display
# --------------------------------------------------

proteins = [
    "Cldn34b",
    "Or51h1",
    "Or8d4",
    "Rsph10b",
    "Sh3rf1"
]


# --------------------------------------------------
# 3. Identify sample columns
# --------------------------------------------------

control_cols = [
    col for col in df.columns
    if str(col).startswith("Control")
]

obe_cols = [
    col for col in df.columns
    if str(col).startswith("Only M.binge")
]


# --------------------------------------------------
# 4. Create figure
# --------------------------------------------------

fig, axes = plt.subplots(
    1,
    5,
    figsize=(15, 5),
    sharey=True
)


# --------------------------------------------------
# 5. Plot each protein
# --------------------------------------------------

for ax, protein in zip(axes, proteins):

    row = df[df["Genes"] == protein]

    control_values = (
        row[control_cols]
        .values
        .flatten()
        .astype(float)
    )

    obe_values = (
        row[obe_cols]
        .values
        .flatten()
        .astype(float)
    )


    # Small horizontal jitter so all animals are visible
    control_jitter = np.linspace(
        -0.07, 0.07, len(control_values)
    )

    obe_jitter = np.linspace(
        -0.07, 0.07, len(obe_values)
    )

    x_control = control_jitter
    x_obe = 1 + obe_jitter


    # --------------------------------------------------
    # Individual values
    # --------------------------------------------------

    ax.scatter(
        x_control,
        control_values,
        s=65,
        color="tab:blue",
        edgecolors="black",
        linewidth=0.8,
        zorder=3
    )

    ax.scatter(
        x_obe,
        obe_values,
        s=65,
        color="tab:orange",
        edgecolors="black",
        linewidth=0.8,
        zorder=3
    )


    # --------------------------------------------------
    # Group means
    # --------------------------------------------------

    control_mean = np.mean(control_values)
    obe_mean = np.mean(obe_values)

    ax.hlines(
        control_mean,
        -0.18,
        0.18,
        color="tab:blue",
        linewidth=2.2,
        zorder=4
    )

    ax.hlines(
        obe_mean,
        0.82,
        1.18,
        color="tab:orange",
        linewidth=2.2,
        zorder=4
    )


    # --------------------------------------------------
    # Formatting
    # --------------------------------------------------

    ax.set_xticks([0, 1])

    ax.set_xticklabels(
        ["O-C", "O-BE"],
        fontsize=10
    )

    ax.set_title(
        protein,
        fontsize=12,
        fontstyle="italic",
        pad=10
    )

    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(13.5, 23.5)

    ax.tick_params(
        axis="y",
        labelsize=10
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# --------------------------------------------------
# 6. Shared Y-axis
# --------------------------------------------------

axes[0].set_ylabel(
    "Log₂ protein intensity",
    fontsize=12
)


# --------------------------------------------------
# 7. Figure title
# --------------------------------------------------

fig.suptitle(
    "Individual Protein Abundance in the ACC",
    fontsize=14,
    y=0.98
)


# --------------------------------------------------
# 8. Legend
# --------------------------------------------------

control_handle = plt.Line2D(
    [0],
    [0],
    marker="o",
    linestyle="None",
    markerfacecolor="tab:blue",
    markeredgecolor="black",
    markersize=7,
    label="O-C"
)

obe_handle = plt.Line2D(
    [0],
    [0],
    marker="o",
    linestyle="None",
    markerfacecolor="tab:orange",
    markeredgecolor="black",
    markersize=7,
    label="O-BE"
)

fig.legend(
    handles=[
        control_handle,
        obe_handle
    ],
    loc="upper center",
    bbox_to_anchor=(0.5, 0.91),
    ncol=2,
    frameon=False,
    fontsize=10
)


# --------------------------------------------------
# 9. Layout
# --------------------------------------------------

plt.tight_layout(
    rect=[0, 0, 1, 0.82]
)


# --------------------------------------------------
# 10. Save figure
# --------------------------------------------------

plt.savefig(
    "ACC_individual_protein_values_final.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()