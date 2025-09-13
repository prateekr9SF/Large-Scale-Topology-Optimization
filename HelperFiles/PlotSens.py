#!/usr/bin/env python3
"""
Plot filtered sensitivities from CSV files in stacked subplots with a shared x-axis.

Usage:
  python plot_sensitivities.py \
      --compliance compliance_sens.csv \
      --volume volume_sens.csv \
      [--cg cg_sens.csv] \
      [--xmin 0 --xmax 50000] \
      --out sensitivities.png \
      --dpi 300
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from typing import List
import re

import matplotlib as mpl
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['agg.path.chunksize'] = 20000  

plt.style.use("seaborn-v0_8-deep")


import re

def offset_to_mathtext(offset: str) -> str:
    """
    Convert offset strings like '1e-2', '1e+03', '×10^-4', '×10⁻⁶' to LaTeX '$\\times 10^{...}$'.
    Returns '' if no usable exponent is found.
    """
    s = offset.strip()
    if not s:
        return ""

    # Normalize unicode minus
    s = s.replace('−', '-')

    # Pattern 1: '1e-04' or '1e+03'
    m = re.fullmatch(r"1e([+-]?\d+)", s)
    if m:
        exp = int(m.group(1))
        return rf"$\times 10^{{{exp}}}$"

    # Pattern 2: '×10^-4' or 'x10^-4' or '10^-4' / '10^-04'
    m = re.search(r"10\^?([+-]?\d+)", s)
    if m:
        exp = m.group(1)
        return rf"$\times 10^{{{exp}}}$"

    # Pattern 3: just an exponent with unicode superscripts, try to recover digits
    sup_map = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻", "0123456789+-")
    t = s.translate(sup_map)
    m = re.search(r"10\^?([+-]?\d+)", t)
    if m:
        return rf"$\times 10^{{{m.group(1)}}}$"

    return ""  # nothing usable


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
    ap.add_argument("--cg", required=False, help="Path to CG sensitivities CSV (optional)")
    ap.add_argument("--xmin", type=int, default=0, help="Optional min element index")
    ap.add_argument("--xmax", type=int, default=None, help="Optional max element index")
    ap.add_argument("--out", default="sensitivities.png", help="Output filename")
    ap.add_argument("--dpi", type=int, default=300, help="Figure DPI")
    args = ap.parse_args()

    # Load compliance and volume data
    comp = pd.read_csv(args.compliance)
    vol  = pd.read_csv(args.volume)

    comp.columns = comp.columns.str.strip()
    vol.columns  = vol.columns.str.strip()

    comp_col = find_col(comp, ["compliance gradient filtered"])
    vol_col  = find_col(vol, ["volume gradient filtered"])

    s_comp = comp[comp_col].astype(float)
    s_vol  = vol[vol_col].astype(float)

    # Load CG data if provided
    cg_present = args.cg is not None
    if cg_present:
        cg = pd.read_csv(args.cg)
        cg.columns = cg.columns.str.strip()
        cgx_col  = find_col(cg, ["cgx gradient"])
        cgy_col  = find_col(cg, ["cgy gradient"])
        cgz_col  = find_col(cg, ["cgz gradient"])
        s_cgx = cg[cgx_col].astype(float)
        s_cgy = cg[cgy_col].astype(float)
        s_cgz = cg[cgz_col].astype(float)

    # Determine slice range
    n = min(len(s_comp), len(s_vol), len(s_cgx)) if cg_present else min(len(s_comp), len(s_vol))
    xmin = max(0, args.xmin)
    xmax = min(n, args.xmax) if args.xmax is not None else n

    s_comp = s_comp.iloc[xmin:xmax].reset_index(drop=True)
    s_vol  = s_vol.iloc[xmin:xmax].reset_index(drop=True)

    if cg_present:
        s_cgx = s_cgx.iloc[xmin:xmax].reset_index(drop=True)
        s_cgy = s_cgy.iloc[xmin:xmax].reset_index(drop=True)
        s_cgz = s_cgz.iloc[xmin:xmax].reset_index(drop=True)

    # Setup subplots
    nrows = 5 if cg_present else 2
    fig, axs = plt.subplots(nrows, 1, figsize=(12, 3*nrows), sharex=True)

    if nrows == 2:
        axs = axs.ravel()
    else:
        axs = axs.ravel()

    axs[0].plot(s_comp, linewidth=0.5, color = "black", alpha = 0.8)
    axs[0].set_ylabel(r"$\partial \mathcal{C}/ \partial \rho$", fontsize=22, fontname="Times New Roman")

    axs[1].plot(s_vol, linewidth=0.5, color = "C0", alpha = 0.8)
    axs[1].set_ylabel(r"$\partial \mathcal{V}/ \partial \rho$", fontsize=22, fontname="Times New Roman")


    if cg_present:
        axs[2].plot(s_cgx, linewidth=0.5, color = "C1", alpha = 0.8)
        axs[2].set_ylabel(r"$\partial CG_x/ \partial \rho$", fontsize=22, fontname="Times New Roman")


        axs[3].plot(s_cgy, linewidth=0.5, color = "C2", alpha = 0.8)
        axs[3].set_ylabel(r"$\partial CG_y/ \partial \rho$", fontsize=22, fontname="Times New Roman")


        axs[4].plot(s_cgz, linewidth=0.5, color = "C3", alpha = 0.8)
        axs[4].set_ylabel(r"$\partial CG_z/ \partial \rho$", fontsize=22, fontname="Times New Roman")

        axs[4].set_xlabel("Element index", fontsize=22, fontname="Times New Roman")

    else:
        axs[1].set_xlabel("Element index", fontsize=22, fontname="Times New Roman")


    # Grid and scientific notation
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))
    formatter.set_useOffset(False) 


    for ax in axs:
        ax.yaxis.set_major_formatter(formatter)

    # draw once so offset texts are computed
    plt.tight_layout()
    plt.draw()



# Move offset into ylabel (LaTeX) and hide the corner text
    for ax in axs:
        off = ax.yaxis.get_offset_text().get_text()
        latex = offset_to_mathtext(off)
        if latex:
            ax.set_ylabel(ax.get_ylabel() + " " + latex)
            ax.yaxis.get_offset_text().set_visible(False)

    for ax in axs:
        ax.yaxis.set_major_formatter(formatter)
        ax.grid(which='major', color='black', linestyle=':', linewidth='0.01')
        ax.minorticks_on()
        ax.grid(which='minor', color='black', linestyle=':', linewidth='0.01')

        ax.tick_params(bottom=True, top=True, left=True, right=True)
        ax.tick_params(which='major', length=10, width=1.2, direction='in')
        ax.tick_params(which='minor', length=5, width=1.2, direction='in')

        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.5)

        #plt.xticks(fontname="Times New Roman", fontsize=20)
        #plt.yticks(fontname="Times New Roman", fontsize=20)        

    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 1.5, Size[1] * 1.5, forward=True)
    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Saved figure to {args.out}")

if __name__ == "__main__":
    main()
