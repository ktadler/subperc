#!/usr/bin/env python3
"""Plot simulated p_c bands (n^{1/3} .. n^{2/3}) against the analytic curve.

Each table<d>reg.tab has columns D, p13_*, p23_*.
  - shaded band = mean over runs of the two cutoffs
  - light dashed = the two cutoff curves
  - solid       = analytic E=1 numerics

    python3 plot_compare.py [outdir]
"""
import colorsys
import glob
import os
import re
import sys
from math import floor, log

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import pandas as pd


def E(p, D, d):
    k = floor(1 / D) - 1
    qk, qk1 = (k + 2) * D - 1, 1 - (k + 1) * D
    s = qk * p**k + qk1 * p**(k + 1)
    return (1 - D) * (d - 3) * p + 2 * p * p * (d - 2) / (1 - p) ** 2 * (
        1 - 2 * D - p + p * D + s
    )


def p_crit(D, d):
    if D == 0:
        return 1 / (d - 1)
    lo = hi = None
    prev = None
    for i in range(1, 999):
        p = i / 1000
        f = E(p, D, d) - 1
        if prev is not None and prev < 0 <= f:
            lo, hi = p - 0.001, p
            break
        prev = f
    if lo is None:
        return 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if E(mid, D, d) - 1 > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def p_crit_random(D, d):
    """Bond percolation + uniform random vertex retention: (1-D)p(d-1)=1."""
    if D >= 1:
        return 1.0
    return min(1.0, 1.0 / ((1.0 - D) * (d - 1)))


def lighten(color, amount=0.45):
    r, g, b = to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    return colorsys.hls_to_rgb(h, min(1.0, l + (1 - l) * amount), s)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "out"
    fig, ax = plt.subplots(figsize=(9, 6.5))
    palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    paths = []
    for path in glob.glob(os.path.join(outdir, "table*reg.tab")):
        df = pd.read_csv(path, sep="\t", nrows=1)
        if any(c.startswith(("p13_", "p12_", "p23_")) for c in df.columns):
            paths.append(path)
    paths.sort(key=lambda p: int(re.search(r"table(\d+)reg", p).group(1)))

    for i, path in enumerate(paths):
        d = int(re.search(r"table(\d+)reg", path).group(1))
        df = pd.read_csv(path, sep="\t")
        c13 = [c for c in df.columns if c.startswith("p13_")]
        c23 = [c for c in df.columns if c.startswith("p23_")]
        lo = df[c13].mean(axis=1)
        hi = df[c23].mean(axis=1)
        Ds = df["D"]
        color = palette[i % len(palette)]
        light = lighten(color, 0.4)

        ax.fill_between(Ds, lo, hi, color=color, alpha=0.18, linewidth=0, zorder=1)
        ax.plot(Ds, lo, linestyle="--", color=light, lw=1.0, zorder=2,
                label=fr"$n^{{1/3}}$, $n^{{2/3}}$ cutoffs  $d={d}$")
        ax.plot(Ds, hi, linestyle="--", color=light, lw=1.0, zorder=2)
        ax.plot(
            Ds, [min(p_crit(float(D), d), 1.0) for D in Ds],
            linestyle="-", color=color, lw=1.8, zorder=3,
            label=fr"structured numerics $d={d}$",
        )
        ax.plot(
            Ds, [p_crit_random(float(D), d) for D in Ds],
            linestyle=":", color=color, lw=1.5, zorder=3, alpha=0.9,
            label=fr"random deletion $d={d}$",
        )

    ax.set_xlabel(r"$D$", fontsize=12)
    ax.set_ylabel(r"$p$", fontsize=12)
    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=7, frameon=False, loc="upper left", ncol=1)
    fig.tight_layout()
    fn = os.path.join(outdir, "pvsD_compare.png")
    fig.savefig(fn, dpi=160, bbox_inches="tight")
    print("wrote", fn)

    # Marker-style figure matching overleaf pvsD.png: analytic + random, sim as markers.
    fig2, ax2 = plt.subplots(figsize=(9, 6.2))
    for i, path in enumerate(paths):
        d = int(re.search(r"table(\d+)reg", path).group(1))
        df = pd.read_csv(path, sep="\t")
        c23 = [c for c in df.columns if c.startswith("p23_")]
        Ds = df["D"]
        p23 = df[c23].mean(axis=1)
        color = palette[i % len(palette)]
        ax2.plot(
            Ds, [min(p_crit(float(D), d), 1.0) for D in Ds],
            linestyle="-", color=color, lw=1.8,
            label=fr"structured $d={d}$",
        )
        ax2.plot(
            Ds, [p_crit_random(float(D), d) for D in Ds],
            linestyle="--", color=color, lw=1.3, alpha=0.85,
            label=fr"random $d={d}$",
        )
        ax2.scatter(
            Ds[::2], p23[::2], s=12, color=color, zorder=4, alpha=0.85,
        )

    ax2.set_xlabel(r"$D$", fontsize=12)
    ax2.set_ylabel(r"$p$", fontsize=12)
    ax2.set_xlim(0, 0.8)
    ax2.set_ylim(0, 1.02)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(fontsize=7, frameon=False, loc="upper left", ncol=2)
    fig2.tight_layout()
    fn2 = os.path.join(outdir, "pvsD.png")
    fig2.savefig(fn2, dpi=160, bbox_inches="tight")
    print("wrote", fn2)

    # One paper-style figure per cutoff that is present in the tables.
    for prefix, latex, fname in (
        ("p13", r"n^{1/3}", "pvsD_n13.png"),
        ("p12", r"n^{1/2}", "pvsD_n12.png"),
        ("p23", r"n^{2/3}", "pvsD_n23.png"),
    ):
        paper_figure(paths, palette, prefix, latex, os.path.join(outdir, fname))
    size_figure(outdir, palette)
    size_vs_random(outdir, palette)


def paper_figure(paths, palette, prefix, latex, outfile):
    """Markers = mean simulated p_c at one cutoff; solid = E=1; dashed = random."""
    if not paths:
        return
    sample = pd.read_csv(paths[0], sep="\t", nrows=1)
    if not any(c.startswith(prefix + "_") for c in sample.columns):
        return
    fig, ax = plt.subplots(figsize=(9, 6.2))
    for i, path in enumerate(paths):
        d = int(re.search(r"table(\d+)reg", path).group(1))
        df = pd.read_csv(path, sep="\t")
        cols = [c for c in df.columns if c.startswith(prefix + "_")]
        if not cols:
            continue
        Ds = df["D"]
        pc = df[cols].mean(axis=1)
        color = palette[i % len(palette)]
        ax.plot(
            Ds, [min(p_crit(float(D), d), 1.0) for D in Ds],
            linestyle="-", color=color, lw=1.8,
            label=fr"structured $d={d}$",
        )
        ax.plot(
            Ds, [p_crit_random(float(D), d) for D in Ds],
            linestyle="--", color=color, lw=1.3, alpha=0.85,
            label=fr"random $d={d}$",
        )
        ax.scatter(Ds[::2], pc[::2], s=12, color=color, zorder=4, alpha=0.85)
    ax.set_xlabel(r"$D$", fontsize=12)
    ax.set_ylabel(r"$p$", fontsize=12)
    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=7, frameon=False, loc="upper left", ncol=2)
    fig.tight_layout()
    fig.savefig(outfile, dpi=160, bbox_inches="tight")
    print("wrote", outfile)


def giant_random(D, p, d):
    """Fraction of vertices in the giant after deleting a fraction D at random
    and occupying each remaining edge with probability p.
    """
    phi = 1.0 - D
    pi = p * phi
    if pi <= 0.0 or pi * (d - 1) <= 1.0:
        return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        rhs = (1.0 - pi) + pi * mid ** (d - 1)
        if rhs > mid:
            lo = mid
        else:
            hi = mid
    eta = 0.5 * (lo + hi)
    return phi * (1.0 - eta ** d)


def size_figure(outdir, palette):
    """Largest component at E=1-gap and E=1+gap, the two proved regimes."""
    paths = sorted(
        glob.glob(os.path.join(outdir, "size*_gap*.tab")),
        key=lambda p: int(re.search(r"size(\d+)_gap", os.path.basename(p)).group(1)),
    )
    if not paths:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    for i, path in enumerate(paths):
        d = int(re.search(r"size(\d+)_gap", os.path.basename(path)).group(1))
        with open(path) as fh:
            header = fh.readline()
        n_match = re.search(r"n=(\d+)", header)
        if n_match is None:
            continue
        n = int(n_match.group(1))
        df = pd.read_csv(path, sep="\t", comment="#")
        lo_cols = [c for c in df.columns if c.startswith("S_lo_")]
        hi_cols = [c for c in df.columns if c.startswith("S_hi_")]
        color = palette[i % len(palette)]
        Ds = df["D"]
        slo = df[lo_cols].mean(axis=1) / log(n)
        shi = df[hi_cols].mean(axis=1) / n
        axes[0].plot(Ds, slo, color=color, lw=1.6, label=fr"$d={d}$")
        axes[1].plot(Ds, shi, color=color, lw=1.6, label=fr"$d={d}$")
    axes[0].set_xlabel(r"$D$")
    axes[1].set_xlabel(r"$D$")
    axes[0].set_ylabel(r"$S/\log n$")
    axes[1].set_ylabel(r"$S/n$")
    for ax in axes:
        ax.set_xlim(0, 0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fn = os.path.join(outdir, "size_gap.png")
    fig.savefig(fn, dpi=160, bbox_inches="tight")
    print("wrote", fn)


def size_vs_random(outdir, palette):
    """Right-hand panel against random deletion of the same fraction D,
    at the same occupation probability p_hi (the p with structured E=1+gap).
    """
    paths = sorted(
        glob.glob(os.path.join(outdir, "size*_gap*.tab")),
        key=lambda p: int(re.search(r"size(\d+)_gap", os.path.basename(p)).group(1)),
    )
    if not paths:
        return
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    for i, path in enumerate(paths):
        d = int(re.search(r"size(\d+)_gap", os.path.basename(path)).group(1))
        with open(path) as fh:
            header = fh.readline()
        n_match = re.search(r"n=(\d+)", header)
        if n_match is None:
            continue
        n = int(n_match.group(1))
        df = pd.read_csv(path, sep="\t", comment="#")
        hi_cols = [c for c in df.columns if c.startswith("S_hi_")]
        color = palette[i % len(palette)]
        Ds = df["D"].to_numpy()
        shi = df[hi_cols].mean(axis=1).to_numpy() / n
        prad = [giant_random(float(D), float(p), d) if p == p else float("nan")
                for D, p in zip(Ds, df["p_hi"])]
        ax.plot(Ds, shi, color=color, lw=1.6, label=fr"spaced $d={d}$")
        ax.plot(Ds, prad, color=color, lw=1.2, ls="--", label=fr"random $d={d}$")
    ax.set_xlabel(r"$D$")
    ax.set_ylabel(r"$S/n$")
    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    fn = os.path.join(outdir, "size_gap_random.png")
    fig.savefig(fn, dpi=160, bbox_inches="tight")
    print("wrote", fn)


if __name__ == "__main__":
    main()
