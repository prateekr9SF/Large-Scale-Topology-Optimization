#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ipopt_log_to_csv.py
Parse an IPOPT log, extract the iteration table, and add a per-iteration
'free_mu_mode' boolean by scanning the log text between iteration rows.

USAGE:
  python ipopt_log_to_csv.py --log /path/to/ipopt.log --out table.csv
"""

import argparse
import os
import re
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

# ---------- utils ----------
_FLOAT_RX = re.compile(r"^[+-]?((\d+(\.\d*)?)|(\.\d+))([eE][+-]?\d+)?$")
def _is_num(tok: str) -> bool:
    t = tok.strip().lower()
    if t in ("nan", "inf", "-inf", "+inf"):
        return True
    return bool(_FLOAT_RX.match(t))

def _lower(s: str) -> str:
    return s.strip().lower()

REQUIRED = ["iter", "objective", "inf_pr", "inf_du"]
OPTIONAL_KNOWN = ["lg(mu)", "||d||", "lg(rg)", "alpha_du", "alpha_pr", "ls", "kkt", "kkt-error", "kappa_d"]

# ---------- parse headers/blocks ----------
def find_all_headers(lines: Sequence[str]) -> List[Tuple[int, List[str]]]:
    headers: List[Tuple[int, List[str]]] = []
    for i, line in enumerate(lines):
        s = _lower(line)
        if "iter" in s and "objective" in s and "inf_pr" in s and "inf_du" in s:
            cols = re.split(r"\s+", line.strip())
            headers.append((i, cols))
    return headers

def parse_block_from_header(lines: Sequence[str], header_idx: int, header_cols: List[str]) -> List[Tuple[Dict[str, Any], int]]:
    col_idx = {c.lower(): i for i, c in enumerate(header_cols)}
    if not all(k in col_idx for k in REQUIRED):
        return []

    rows: List[Tuple[Dict[str, Any], int]] = []
    j = header_idx + 1
    n = len(lines)
    while j < n:
        line = lines[j].rstrip("\n")
        if not line.strip():
            j += 1
            continue
        s = _lower(line)
        if "iter" in s and "objective" in s and "inf_pr" in s and "inf_du" in s:
            break  # next header

        toks = re.split(r"\s+", line.strip())
        try:
            itv = int(toks[col_idx["iter"]])
        except Exception:
            j += 1
            continue

        ok = True
        for req in REQUIRED:
            if req == "iter":
                continue
            pos = col_idx[req]
            if pos >= len(toks) or not _is_num(toks[pos]):
                ok = False
                break
        if not ok:
            j += 1
            continue

        rec: Dict[str, Any] = {"iter": itv}
        for name, pos in col_idx.items():
            if pos < len(toks):
                val = toks[pos]
                rec[name] = int(val) if name == "iter" else (float(val) if _is_num(val) else val)
        rows.append((rec, j))
        j += 1

    return rows

def parse_ipopt_table_with_positions(lines: Sequence[str]) -> List[Tuple[Dict[str, Any], int]]:
    """Return list of (row_dict, line_index). Keeps the last row per iter if repeated."""
    all_rows: List[Tuple[Dict[str, Any], int]] = []
    for hdr_idx, hdr_cols in find_all_headers(lines):
        all_rows.extend(parse_block_from_header(lines, hdr_idx, hdr_cols))
    if not all_rows:
        return []
    last: Dict[int, Tuple[Dict[str, Any], int]] = {}
    for rec, ln in all_rows:
        last[int(rec["iter"])] = (rec, ln)
    items = sorted(last.items(), key=lambda kv: kv[0])  # by iter
    return [pair for _, pair in items]  # (rec, ln)

# ---------- free μ mode detection ----------
_RX_FREE = [
    re.compile(r"\bfree\b.{0,80}\b(m[uµ])\b.{0,80}\bmode\b", re.IGNORECASE | re.DOTALL),
    re.compile(r"\bfree\s*(m[uµ])\s*mode\b", re.IGNORECASE),
]

def _window_has_free_mode(text: str) -> bool:
    for rx in _RX_FREE:
        if rx.search(text):
            return True
    return False

def detect_free_mu_flags(lines: Sequence[str], iter_positions: List[Tuple[Dict[str, Any], int]]) -> Dict[int, bool]:
    """
    Accepts iter_positions as [(row_dict, line_index), ...].
    Scans between each row and the next for any "free mu mode" message.
    Returns {iter_value: bool}.
    """
    # Normalize to [(iter_value, pos)]
    norm = []
    for rec, pos in iter_positions:
        it = int(rec["iter"])
        norm.append((it, pos))

    its = [it for it, _ in norm]
    poss = [pos for _, pos in norm]
    # Build windows
    ranges = {}
    for k, it in enumerate(its):
        start = poss[k]
        end = poss[k+1] if k+1 < len(poss) else len(lines)
        ranges[it] = (start, end)

    flags = {}
    for it, (s, e) in ranges.items():
        win_text = " ".join(lines[s:e])
        flags[it] = _window_has_free_mode(win_text)
    return flags

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Extract IPOPT iteration table to CSV with free_mu_mode.")
    ap.add_argument("--log", required=True, help="Path to IPOPT log file")
    ap.add_argument("--out", required=True, help="Output CSV path")
    args = ap.parse_args()

    with open(args.log, "r", errors="ignore") as f:
        lines = f.readlines()

    rows_with_pos = parse_ipopt_table_with_positions(lines)
    if not rows_with_pos:
        raise SystemExit("No IPOPT iteration table parsed. Check the log format.")

    rows = [rec for rec, _ in rows_with_pos]

    # Build dataframe
    all_keys = set()
    for r in rows:
        all_keys.update(r.keys())
    cols = [k for k in REQUIRED if k in all_keys]
    cols += [k for k in OPTIONAL_KNOWN if k in all_keys and k not in cols]
    cols += [k for k in sorted(all_keys) if k not in cols]

    df = pd.DataFrame(rows)
    if "iter" in df.columns:
        df["iter"] = pd.to_numeric(df["iter"], errors="coerce").astype("Int64")

    # free μ per iteration
    free_flags = detect_free_mu_flags(lines, rows_with_pos)
    df["free_mu_mode"] = df["iter"].astype(int).map(free_flags)

    # order columns
    ordered = [c for c in cols if c in df.columns]
    ordered.append("free_mu_mode")
    df = df[ordered]

    # write csv
    outdir = os.path.dirname(os.path.abspath(args.out))
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"Wrote {args.out} with {len(df)} rows and columns: {list(df.columns)}")

if __name__ == "__main__":
    main()
