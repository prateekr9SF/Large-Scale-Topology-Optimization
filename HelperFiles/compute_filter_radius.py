#!/usr/bin/env python3
"""
Simplified script:
- Reads a CSV with histogram bins (bin start, bin end, count).
- Computes the plain mean of bin centers (unweighted).
- Suggests filter radii as multiples of this mean.
- Plots the histogram.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Sequence

def compute_and_plot_filter_radius_from_csv(
    file_path: str,
    multiplier_values: Sequence[float] = (1.5, 2.0, 2.5, 3.0),
    out_png: str = "EdgeLengthDist.png",
    dpi: int = 300,
) -> None:
    # Load CSV
    df = pd.read_csv(file_path, header=0, index_col=[0, 1])

    bin_edges, counts = [], []
    for (l, r), c in df.iloc[:, 0].items():
        try:
            l = float(l)
            r = float(r)
            bin_center = 0.5 * (l + r)
        except Exception:
            continue
        bin_edges.append(bin_center)
        counts.append(float(c))

    bin_edges = np.array(bin_edges)
    counts = np.array(counts)

    if len(bin_edges) == 0:
        raise ValueError("No valid bins found in CSV.")

    # Plain (unweighted) mean
    mean_edge_length = bin_edges.mean()

    # Filter radii suggestions
    filter_radii = {f"{m}×": m * mean_edge_length for m in multiplier_values}

    # Plot
    plt.figure(figsize=(8, 5))
    plt.bar(bin_edges, counts, width=(bin_edges[1] - bin_edges[0] if len(bin_edges) > 1 else 0.01),
            edgecolor='k', align='center', alpha=0.8)
    plt.axvline(mean_edge_length, color='red', linestyle='--',
                label=f'Mean = {mean_edge_length:.3f} m')

    for label, value in filter_radii.items():
        plt.axvline(value, linestyle=':', label=f'{label}Mean = {value:.3f} m')

    plt.xlabel("Edge Length (m)")
    plt.ylabel("Count")
    plt.title("Edge Length Distribution and Filter Radius Suggestions")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)

    # Print results
    print(f"Mean edge length: {mean_edge_length:.5f} m")
    for label, value in filter_radii.items():
        print(f"Filter radius ({label}): {value:.5f} m")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute filter radius from edge length histogram CSV")
    parser.add_argument("csv_path", type=str, help="Path to the edge length histogram CSV file")
    parser.add_argument("--out", type=str, default="EdgeLengthDist.png", help="Output PNG path")
    parser.add_argument("--dpi", type=int, default=300, help="PNG DPI")
    parser.add_argument("--multipliers", type=float, nargs="+", default=[1.5, 2.0, 2.5, 3.0],
                        help="Multipliers applied to mean edge length")
    args = parser.parse_args()

    compute_and_plot_filter_radius_from_csv(
        args.csv_path,
        multiplier_values=args.multipliers,
        out_png=args.out,
        dpi=args.dpi,
    )
