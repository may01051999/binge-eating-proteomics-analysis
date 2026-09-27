import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


# --------------------------------------------------
# 1. Load data
# --------------------------------------------------

df = pd.read_excel("Proteome - maternal BE.xlsx")


# --------------------------------------------------
# 2. Identify sample columns
# --------------------------------------------------

control_cols = [
    col for col in df.columns
    if str(col).startswith("Control")
]

obe_cols = [
    col for col in df.columns
    if str(col).startswith("Only M.binge")
]

sample_cols = control_cols + obe_cols


# --------------------------------------------------
# 3. Create protein × sample matrix
# --------------------------------------------------

protein_matrix = df[sample_cols].apply(
    pd.to_numeric,
    errors="coerce"
)

protein_matrix = protein_matrix.dropna()


# --------------------------------------------------
# 4. Transpose matrix
# Rows = samples
# Columns = proteins
# --------------------------------------------------

X = protein_matrix.T


# --------------------------------------------------
# 5. Check data
# --------------------------------------------------

print("Number of proteins used for PCA:", X.shape[1])
print("Number of samples used for PCA:", X.shape[0])


# --------------------------------------------------
# 6. PCA
# Data are already Log2-transformed
# --------------------------------------------------

pca = PCA(n_components=2)

pca_scores = pca.fit_transform(X)

explained_variance = pca.explained_variance_ratio_ * 100


print(
    f"PC1 explained variance: "
    f"{explained_variance[0]:.1f}%"
)

print(
    f"PC2 explained variance: "
    f"{explained_variance[1]:.1f}%"
)

print(
    f"Total explained variance (PC1 + PC2): "
    f"{explained_variance[0] + explained_variance[1]:.1f}%"
)


# --------------------------------------------------
# 7. Group information
# --------------------------------------------------

groups = [
    "O-C", "O-C", "O-C", "O-C", "O-C",
    "O-BE", "O-BE", "O-BE", "O-BE", "O-BE"
]

sample_numbers = [
    "1", "2", "3", "4", "5",
    "1", "2", "3", "4", "5"
]


# --------------------------------------------------
# 8. Identify group indices
# --------------------------------------------------

control_idx = [
    i for i, group in enumerate(groups)
    if group == "O-C"
]

obe_idx = [
    i for i, group in enumerate(groups)
    if group == "O-BE"
]


# --------------------------------------------------
# 9. Create PCA plot
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(8.5, 6.5))


# O-C samples
ax.scatter(
    pca_scores[control_idx, 0],
    pca_scores[control_idx, 1],
    s=220,
    color="tab:blue",
    edgecolors="black",
    linewidth=1,
    label="O-C",
    zorder=3
)


# O-BE samples
ax.scatter(
    pca_scores[obe_idx, 0],
    pca_scores[obe_idx, 1],
    s=220,
    color="tab:orange",
    edgecolors="black",
    linewidth=1,
    label="O-BE",
    zorder=3
)


# --------------------------------------------------
# 10. Add sample number inside each point
# --------------------------------------------------

for i, number in enumerate(sample_numbers):

    ax.text(
        pca_scores[i, 0],
        pca_scores[i, 1],
        number,
        ha="center",
        va="center",
        fontsize=9,
        fontweight="bold",
        color="white",
        zorder=4
    )


# --------------------------------------------------
# 11. Axis labels
# --------------------------------------------------

ax.set_xlabel(
    f"PC1 ({explained_variance[0]:.1f}% variance)",
    fontsize=12
)

ax.set_ylabel(
    f"PC2 ({explained_variance[1]:.1f}% variance)",
    fontsize=12
)


# --------------------------------------------------
# 12. Reference lines
# --------------------------------------------------

ax.axhline(
    y=0,
    linewidth=0.8,
    linestyle="--",
    alpha=0.35
)

ax.axvline(
    x=0,
    linewidth=0.8,
    linestyle="--",
    alpha=0.35
)


# --------------------------------------------------
# 13. Title
# --------------------------------------------------

ax.set_title(
    "PCA of ACC Proteomic Profiles",
    fontsize=15,
    pad=15
)


# --------------------------------------------------
# 14. Legend
# --------------------------------------------------

ax.legend(
    loc="lower left",
    fontsize=11,
    frameon=True,
    edgecolor="black",
    facecolor="white",
    framealpha=1
)


# --------------------------------------------------
# 15. Clean appearance
# --------------------------------------------------

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.tick_params(
    axis="both",
    labelsize=10
)

ax.margins(
    x=0.06,
    y=0.10
)


# --------------------------------------------------
# 16. Layout and save
# --------------------------------------------------

plt.tight_layout()

plt.savefig(
    "ACC_PCA_final.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()