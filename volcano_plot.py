import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# --------------------------------------------------
# 1. Load data
# --------------------------------------------------

df = pd.read_excel("Proteome - maternal BE.xlsx")

# Calculate -log10(p-value)
df["-log10_p"] = -np.log10(df["p value"])


# --------------------------------------------------
# 2. Define thresholds
# --------------------------------------------------

p_threshold = 0.001
fc_threshold = 1.5


# --------------------------------------------------
# 3. Select proteins passing both thresholds
# --------------------------------------------------

selected = df[
    (df["p value"] < p_threshold) &
    (df["Fold change"].abs() >= fc_threshold)
].copy()

# Separate proteins according to direction of change
up = selected[selected["Fold change"] >= fc_threshold]
down = selected[selected["Fold change"] <= -fc_threshold]


# Print selected proteins
print("\nProteins passing both thresholds:")
print(selected[["Genes", "p value", "Fold change"]])

print("\nNumber of selected proteins:", len(selected))
print("Higher in O-BE:", len(up))
print("Lower in O-BE:", len(down))


# --------------------------------------------------
# 4. Create Volcano Plot
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 7))


# All proteins
ax.scatter(
    df["Fold change"],
    df["-log10_p"],
    s=40,
    alpha=0.40,
    label="Other proteins"
)


# Proteins higher in O-BE
ax.scatter(
    up["Fold change"],
    up["-log10_p"],
    s=100,
    edgecolors="black",
    linewidth=1,
    label="Higher in O-BE",
    zorder=3
)


# Proteins lower in O-BE
ax.scatter(
    down["Fold change"],
    down["-log10_p"],
    s=100,
    edgecolors="black",
    linewidth=1,
    label="Lower in O-BE",
    zorder=3
)


# --------------------------------------------------
# 5. Threshold lines
# --------------------------------------------------

# p = 0.001
ax.axhline(
    y=-np.log10(p_threshold),
    linestyle="--",
    linewidth=1.2
)

# log2FC = +1.5
ax.axvline(
    x=fc_threshold,
    linestyle="--",
    linewidth=1.2
)

# log2FC = -1.5
ax.axvline(
    x=-fc_threshold,
    linestyle="--",
    linewidth=1.2
)


# --------------------------------------------------
# 6. Protein labels
# --------------------------------------------------

# Manual positions to prevent overlap
label_offsets = {
    "Cldn34b": (12, 8),
    "Or51h1": (12, 8),
    "Or8d4": (-12, 8),
    "Rsph10b": (-12, 8),
    "AABR07040789.1": (14, 14),
    "Sh3rf1": (12, 2)
}


for _, row in selected.iterrows():

    gene = str(row["Genes"])
    x = row["Fold change"]
    y = row["-log10_p"]

    dx, dy = label_offsets.get(gene, (8, 8))

    ax.annotate(
        gene,
        xy=(x, y),
        xytext=(dx, dy),
        textcoords="offset points",
        fontsize=9.5,
        ha="left" if x > 0 else "right",
        va="center",
        zorder=4
    )


# --------------------------------------------------
# 7. Axis labels and title
# --------------------------------------------------

ax.set_xlabel(
    "log₂ fold change (O-BE − O-C)",
    fontsize=12
)

ax.set_ylabel(
    "−log₁₀(p-value)",
    fontsize=12
)

ax.set_title(
    "Differential Protein Abundance in the ACC: O-BE vs O-C",
    fontsize=14,
    pad=14
)


# --------------------------------------------------
# 8. Plot formatting
# --------------------------------------------------

# Extra horizontal space for labels
ax.set_xlim(
    df["Fold change"].min() - 1.5,
    df["Fold change"].max() + 1.5
)

# Extra vertical space so upper labels do not touch the top
ax.set_ylim(
    1.5,
    df["-log10_p"].max() + 1.2
)

ax.tick_params(
    axis="both",
    labelsize=10
)

ax.legend(
    frameon=False,
    fontsize=10,
    loc="upper left"
)

# Cleaner scientific figure
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()


# --------------------------------------------------
# 9. Save high-resolution figure
# --------------------------------------------------

plt.savefig(
    "ACC_volcano_plot_final.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()