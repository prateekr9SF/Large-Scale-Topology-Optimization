
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def compute_and_plot_filter_radius_from_csv(file_path, multiplier_values=[1.5, 2.0]):
    df = pd.read_csv(file_path, header=0, index_col=[0, 1])

    bin_edges = []
    counts = []

    for idx, count in df.iloc[:, 0].items():
        try:
            left = float(idx[0])
            right = float(idx[1])
            bin_center = 0.5 * (left + right)
            bin_edges.append(bin_center)
            counts.append(count)
        except ValueError:
            continue

    bin_edges = np.array(bin_edges)
    counts = np.array(counts)

    mean_edge_length = np.average(bin_edges, weights=counts)
    filter_radii = {f"{m}×": m * mean_edge_length for m in multiplier_values}

    plt.figure(figsize=(8, 5))
    plt.bar(bin_edges, counts, width=0.015, edgecolor='k', align='center', alpha=0.8)
    plt.axvline(x=mean_edge_length, color='red', linestyle='--', label=f'Mean = {mean_edge_length:.3f} mm')

    for label, value in filter_radii.items():
        plt.axvline(x=value, linestyle=':', label=f'{label}Mean = {value:.3f} mm')

    plt.xlabel("Edge Length (mm)")
    plt.ylabel("Count")
    plt.title("Edge Length Distribution and Filter Radius Suggestions")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"Mean edge length: {mean_edge_length:.5f} mm")
    for label, value in filter_radii.items():
        print(f"Filter radius ({label}): {value:.5f} mm")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute filter radius from edge length histogram CSV")
    parser.add_argument("csv_path", type=str, help="Path to the edge length histogram CSV file")
    args = parser.parse_args()

    compute_and_plot_filter_radius_from_csv(args.csv_path)
