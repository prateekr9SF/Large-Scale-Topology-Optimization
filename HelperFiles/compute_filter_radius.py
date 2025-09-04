#!/usr/bin/env python3
"""
Simplified script:
- Reads a CSV with histogram bins (bin start, bin end, count).
- Computes the plain mean of bin centers (unweighted).
- Plots the histogram with the mean marked.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use("Agg")

plt.style.use("seaborn-v0_8-deep")

def compute_and_plot_filter_radius_from_csv(
    file_path: str,
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

    # Plot
    fig1 = plt.figure()
    ax1 = fig1.gca()

    plt.bar(bin_edges, counts, width=(bin_edges[1] - bin_edges[0] if len(bin_edges) > 1 else 0.01),
            edgecolor='k', align='center', alpha=0.8)
    plt.axvline(mean_edge_length, color='black', linestyle='--',
                label=f'Mean = {mean_edge_length:.3f} m')

    plt.xlabel("Edge Length (m)", fontsize=22, fontname="Times New Roman")
    plt.ylabel("Count", fontsize=22, fontname="Times New Roman")
    plt.legend(frameon=False, loc ='upper right', prop={'size': 18, 'family': 'Times New Roman'}, ncol=2)

    ax1.grid(which='major', color='black', linestyle=':', linewidth='0.01')
    ax1.minorticks_on()
    ax1.grid(which='minor', color='black', linestyle=':', linewidth='0.01')

    ax1.tick_params(bottom=True, top=True, left=True, right=True)
    ax1.tick_params(which='major', length=10, width=1.2, direction='in')
    ax1.tick_params(which='minor', length=5, width=1.2, direction='in')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)

    plt.xticks(fontname="Times New Roman", fontsize=20)
    plt.yticks(fontname="Times New Roman", fontsize=20)
    ax1.set_xlim(0.0, 0.5)


    plt.tight_layout()
    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 1.5, Size[1] * 1.5, forward=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)

    # Print results
    print(f"Mean edge length: {mean_edge_length:.5f} m")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute mean edge length from histogram CSV and plot")
    parser.add_argument("csv_path", type=str, help="Path to the edge length histogram CSV file")
    parser.add_argument("--out", type=str, default="EdgeLengthDist.png", help="Output PNG path")
    parser.add_argument("--dpi", type=int, default=300, help="PNG DPI")
    args = parser.parse_args()

    compute_and_plot_filter_radius_from_csv(
        args.csv_path,
        out_png=args.out,
        dpi=args.dpi,
    )
