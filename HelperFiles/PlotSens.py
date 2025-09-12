#!/usr/bin/env python3
"""
Plot filtered sensitivities from CSV files in stacked subplots with a shared x-axis.

Usage:
  python plot_sensitivities.py \
      --compliance compliance_sens.csv \
      --volume volume_sens.csv \
      --cg cg_sens.csv \
      --out sensitivities.png \
      --dpi 300
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from typing import List

plt.style.use("seaborn-v0_8-deep")

def find_col(df: pd.DataFrame, must_contain_any: List[str]) -> str:
    cols = {c: c.strip() for c in df.columns}
    lowered = {c: c.strip().lower() for c in df.columns}
    patterns = [p.lower().split() for p in must_contain_any]

    for original, low in lowered.items():
        for tokens in patterns:
            if all(tok in low for tok in tokens):
                return cols[original]
    raise KeyError(
        f"Could not find any column containing one of {must_contain_any}.\n"
        f"Available columns: {list(df.columns)}"
    )

def main():
    ap = argparse.ArgumentParser(description="Plot sensitivities from CSVs in stacked subplots.")
    ap.add_argument("--compliance", required=True, help="Path to compliance_*.csv")
    ap.add_argument("--volume", required=True, help="Path to volume_*.csv")
    ap.add_argument("--cg", required=True, help="Path to CG sensitivities CSV")
    ap.add_argument("--out", default="sensitivities.png", help="Output image filename")
    ap.add_argument("--dpi", type=int, default=200, help="Figure DPI")
    ap.add_argument("--xmax", type=int, default=None, help="Optional max element index to plot")
    ap.add_argument("--xmin", type=int, default=0, help="Optional min element index to plot")
    args = ap.parse_args()

    comp = pd.read_csv(args.compliance)
    vol  = pd.read_csv(args.volume)
    cg   = pd.read_csv(args.cg)

    comp.columns = comp.columns.str.strip()
    vol.columns  = vol.columns.str.strip()
    cg.columns   = cg.columns.str.strip()

    comp_col = find_col(comp, ["compliance gradient filtered"])
    vol_col  = find_col(vol, ["volume gradient filtered"])
    cgx_col  = find_col(cg, ["cgx gradient"])
    cgy_col  = find_col(cg, ["cgy gradient"])
    cgz_col  = find_col(cg, ["cgz gradient"])

    s_comp = comp[comp_col].astype(float)
    s_vol  = vol[vol_col].astype(float)
    s_cgx  = cg[cgx_col].astype(float)
    s_cgy  = cg[cgy_col].astype(float)
    s_cgz  = cg[cgz_col].astype(float)

    n = min(len(s_comp), len(s_vol), len(s_cgx), len(s_cgy), len(s_cgz))
    xmin = max(0, args.xmin)
    xmax = min(n, args.xmax) if args.xmax is not None else n

    s_comp = s_comp.iloc[xmin:xmax].reset_index(drop=True)
    s_vol  = s_vol.iloc[xmin:xmax].reset_index(drop=True)
    s_cgx  = s_cgx.iloc[xmin:xmax].reset_index(drop=True)
    s_cgy  = s_cgy.iloc[xmin:xmax].reset_index(drop=True)
    s_cgz  = s_cgz.iloc[xmin:xmax].reset_index(drop=True)

    fig, axs = plt.subplots(5, 1, figsize=(12, 12), sharex=True)

    axs[0].plot(s_comp, linewidth=0.5)
    axs[0].set_ylabel(r"$\partial \mathcal{C}/ \partial \rho$")
    #axs[0].set_title("Compliance, Volume, and CG Sensitivities (Filtered)")

    axs[1].plot(s_vol, linewidth=0.5)
    axs[1].set_ylabel(r"$\partial \mathcal{V}/ \partial \rho$")

    axs[2].plot(s_cgx, linewidth=0.5)
    axs[2].set_ylabel("CGx")

    axs[3].plot(s_cgy, linewidth=0.5)
    axs[3].set_ylabel("CGy")

    axs[4].plot(s_cgz, linewidth=0.5)
    axs[4].set_ylabel("CGz")
    axs[4].set_xlabel("Element index")

    # Apply scientific notation to all y-axes
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    for ax in axs:
        #ax.yaxis.set_major_formatter(formatter)
        #ax.grid(True, linestyle="--", alpha=0.6)
        ax.grid(which='major', color='black', linestyle=':', linewidth='0.01')
        ax.minorticks_on()
        ax.grid(which='minor', color='black', linestyle=':', linewidth='0.01')


    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 1.5, Size[1] * 1.5, forward=True)
    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Saved figure to Sensitivity.png")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
