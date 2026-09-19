"""
plot_mineral_maps.py
=====================
Turns the output of mineral_id.py into two spatial figures:

  1. abundance_maps_grid.png
     A 14-panel grid, one per endmember, each abundance map titled with its
     #1 mineral match and SAD score -- lets you visually cross-check strong
     abundance regions against the mineral ID (Next Steps item 4 in the
     handoff doc).

  2. dominant_mineral_map.png
     A single classified map of the whole scene: for every pixel, find
     which endmember has the highest abundance there, then color that
     pixel by that endmember's #1 mineral match. Comes with a legend.
     NOTE: this is a simplification -- it discards the actual mixing
     fractions and shows only the single dominant material per pixel, so
     it's useful for a quick visual overview, not a substitute for the
     abundance maps themselves.

Requires:
  - endmember_mineral_matches.csv  (produced by mineral_id.py)
  - abundance_maps.npy             (your FCLS output, shape (250, 190, 14) --
                                     save it with:
                                     np.save('abundance_maps.npy', abundance_maps))

Run:
    python plot_mineral_maps.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# ============================== CONFIG ======================================

MATCHES_CSV = "endmember_mineral_matches_p12.csv"
ABUNDANCE_MAPS_NPY = "abundance_maps_p12.npy"

GRID_OUTPUT_PNG = "abundance_maps_grid_p12.png"
DOMINANT_MAP_OUTPUT_PNG = "dominant_mineral_map_p12.png"

# ============================================================================


def load_top_matches(csv_path):
    """Return dict: endmember index (1-based) -> (mineral name, sad_degrees)
    using only the #1 rank match per endmember."""
    df = pd.read_csv(csv_path)
    top1 = df[df["rank"] == 1].set_index("endmember")
    return {
        idx: (row["mineral"], row["sad_degrees"])
        for idx, row in top1.iterrows()
    }


def plot_abundance_grid(abundance_maps, top_matches, n_cols=4):
    p = abundance_maps.shape[2]
    n_rows = int(np.ceil(p / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.4 * n_rows))
    axes = np.array(axes).reshape(-1)

    for i in range(p):
        ax = axes[i]
        im = ax.imshow(abundance_maps[:, :, i], cmap="viridis", vmin=0, vmax=1)
        mineral, sad_deg = top_matches.get(i + 1, ("unmatched", np.nan))
        title = f"Endmember {i+1}\n{mineral}"
        if not np.isnan(sad_deg):
            title += f" (SAD={sad_deg:.1f}\u00b0)"
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    for j in range(p, len(axes)):
        axes[j].axis("off")

    fig.suptitle("Endmember abundance maps with USGS mineral matches", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(GRID_OUTPUT_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved {GRID_OUTPUT_PNG}")


def plot_dominant_mineral_map(abundance_maps, top_matches):
    p = abundance_maps.shape[2]
    dominant_idx = np.argmax(abundance_maps, axis=2)  # (rows, cols), values 0..p-1

    # Build a distinct color per endmember and matching legend labels.
    cmap = plt.get_cmap("tab20", p)
    colors = [cmap(i) for i in range(p)]
    listed_cmap = ListedColormap(colors)

    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(dominant_idx, cmap=listed_cmap, vmin=-0.5, vmax=p - 0.5)
    ax.set_title("Dominant mineral per pixel (highest-abundance endmember)")
    ax.set_xticks([])
    ax.set_yticks([])

    legend_handles = []
    for i in range(p):
        mineral, sad_deg = top_matches.get(i + 1, ("unmatched", np.nan))
        label = f"E{i+1}: {mineral}"
        if not np.isnan(sad_deg):
            label += f" ({sad_deg:.1f}\u00b0)"
        legend_handles.append(Patch(facecolor=colors[i], label=label))

    ax.legend(
        handles=legend_handles,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8,
        title="Endmember: mineral match (SAD)",
    )
    fig.tight_layout()
    fig.savefig(DOMINANT_MAP_OUTPUT_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {DOMINANT_MAP_OUTPUT_PNG}")


def main():
    print("Loading top-1 mineral matches...")
    top_matches = load_top_matches(MATCHES_CSV)

    print("Loading abundance maps...")
    abundance_maps = np.load(ABUNDANCE_MAPS_NPY)
    print(f"  abundance_maps shape: {abundance_maps.shape}")

    plot_abundance_grid(abundance_maps, top_matches)
    plot_dominant_mineral_map(abundance_maps, top_matches)


if __name__ == "__main__":
    main()