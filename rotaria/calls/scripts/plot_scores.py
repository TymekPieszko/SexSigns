#!/usr/bin/env python3
"""Score histograms from a bcftools-query table, invariant and variant sites in
separate panels (invariant blue, variant red), streamed so memory stays flat
regardless of table size.

Usage: plot_scores.py <table.tsv> <out.png> [--calc_thresholds <file>]

table.tsv columns: TYPE  DP  MQ  QUAL  then per sample: GT,DP,GQ,SP
  TYPE from bcftools %TYPE: 'REF' = invariant, anything else = variant

Genotypes whose GT is missing are skipped: `bcftools filter -S .` blanks only
the GT subfield, leaving DP/GQ/SP in place, so a masked genotype would
otherwise still contribute its (failing) annotations to these histograms and
make a working filter look inert.

With --calc_thresholds, the 95th percentile of INFO/DP is written to <file> as
a bare number (so the filter rule can read it with cat) and drawn on the
DP_total panels. Without it the same plots are produced and no thresholds file
is written, so the script serves both the pre- and post-filter rules.
"""
import sys, os
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PCT = 0.95   # percentile of INFO/DP used as the site depth ceiling

args = sys.argv[1:]
calc = None
if "--calc_thresholds" in args:
    i = args.index("--calc_thresholds")
    calc = args[i + 1]
    del args[i:i + 2]
table, png = args[0], args[1]
title = os.path.basename(png).rsplit(".", 1)[0]

# metric -> bin width (values are binned as int(value / width))
METRICS = [("DP_total", 10), ("DP_sample", 1), ("QUAL", 5),
           ("MQ", 1), ("GQ_sample", 1), ("SP_sample", 1)]
WIDTH = dict(METRICS)

inv = {m: Counter() for m, _ in METRICS}   # invariant sites
var = {m: Counter() for m, _ in METRICS}   # variant sites
site_dp = Counter()                        # INFO/DP, for the depth ceiling

def add(cnt, metric, val):
    cnt[metric][int(val / WIDTH[metric])] += 1

with open(table) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        target = inv if p[0] == "REF" else var
        if p[1] != ".":
            v = float(p[1]); add(target, "DP_total", v); site_dp[int(v)] += 1
        if p[2] != ".":
            add(target, "MQ", float(p[2]))
        if p[3] != ".":
            add(target, "QUAL", float(p[3]))
        # per-sample GT,DP,GQ,SP -> per-site means over CALLED genotypes only
        sd = sg = ss = 0.0
        nd = ng = ns = 0
        for quad in p[4:]:
            gt, d, g, s = quad.split(",")
            if gt in ("./.", ".", ".|."): continue     # masked: skip entirely
            if d != ".": sd += float(d); nd += 1
            if g != ".": sg += float(g); ng += 1
            if s != ".": ss += float(s); ns += 1
        if nd: add(target, "DP_sample", sd / nd)
        if ng: add(target, "GQ_sample", sg / ng)
        if ns:
            m = ss / ns
            if 1 <= m <= 50:            # course convention for SP
                add(target, "SP_sample", m)

def quantile(h, q):
    n = sum(h.values())
    if not n:
        return 0
    cum = 0
    for k in sorted(h):
        cum += h[k]
        if cum >= q * n:
            return k
    return max(h)

dp_max = quantile(site_dp, PCT)

def draw(ax, cnt, metric, width, colour, label):
    if not cnt:
        ax.set_title(f"{metric} - {label} (no sites)", fontsize=10)
        return
    total = sum(cnt.values())
    cum, hi = 0, max(cnt)
    for b in sorted(cnt):               # cap x at this class's 99.5th percentile
        cum += cnt[b]
        if cum >= 0.995 * total:
            hi = b
            break
    bins = [b for b in sorted(cnt) if b <= hi]
    ax.bar([b * width for b in bins], [cnt[b] for b in bins],
           width=width, align="edge", color=colour)
    ax.set_title(f"{metric} - {label}", fontsize=10)
    ax.set_xlabel(metric, fontsize=8)
    ax.set_ylabel("sites", fontsize=8)
    ax.tick_params(labelsize=8)
    ax.text(0.97, 0.93, f"n={total:,}", transform=ax.transAxes,
            ha="right", va="top", fontsize=7)
    if calc and metric == "DP_total":
        ax.axvline(dp_max, color="black", ls="--", lw=1)
        ax.text(dp_max, ax.get_ylim()[1] * 0.75, f" p{PCT:.0%}={dp_max}", fontsize=7)

fig, axes = plt.subplots(3, 4, figsize=(16, 9))
flat = axes.flat
for metric, width in METRICS:
    draw(next(flat), inv[metric], metric, width, "steelblue", "invariant")
    draw(next(flat), var[metric], metric, width, "firebrick", "variant")
fig.suptitle(f"{title}: invariant (blue) and variant (red) sites")
fig.tight_layout()
os.makedirs(os.path.dirname(png) or ".", exist_ok=True)
fig.savefig(png, dpi=150)

if calc:
    with open(calc, "w") as out:
        out.write(f"{dp_max}\n")
    print(f"INFO/DP p{PCT:.0%} = {dp_max}")

n_site = sum(site_dp.values())
print(f"sites: {n_site:,}; "
      f"invariant={sum(inv['DP_total'].values()):,} variant={sum(var['DP_total'].values()):,}")

