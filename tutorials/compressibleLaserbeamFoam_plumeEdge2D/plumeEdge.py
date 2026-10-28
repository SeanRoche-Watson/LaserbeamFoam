#!/usr/bin/env python3
"""Plots for the planar 2D plume-edge case: eos vs mass (or any set of runs).

    python3 plumeEdge.py <outdir> label=runs/<a> label=runs/<b> [--times 2e-09,5e-09,1e-08,2e-08]

Needs reconstructed time directories with the cell centres written:
    postProcess -funcs '(writeCellCentres writeCellVolumes)'     (in each run)
Writes
  fields_<label>.png    T, p, |U|, rho/rhoM (rhoM = written rhoM, or sum c.* for eos runs) at the chosen times
  fronts.png            lateral plume / shock front along the wall, height on the axis, max |U|, max T,
                        and the solver's per-step max |rho/rhoM - 1| line, vs time, all runs
  fronts.csv            the front data
"""
import os, re, sys, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
from plumeFront import read_field, times_of, VAP

OUT = sys.argv[1]
args = [a for a in sys.argv[2:] if not a.startswith("--")]
RUNS = [a.split("=", 1) for a in args]
TS = ["2e-09", "5e-09", "1e-08", "2e-08"]
for a in sys.argv[2:]:
    if a.startswith("--times="):
        TS = a.split("=", 1)[1].split(",")
os.makedirs(OUT, exist_ok=True)


def load(run, d):
    td = os.path.join(run, d)
    C = read_field(os.path.join(td, "C"))
    n = len(C)
    g = lambda f: read_field(os.path.join(td, f), n)
    s = dict(x=C[:, 0], y=C[:, 1], T=g("T"), p=g("p"), U=np.sqrt((g("U")**2).sum(1)), rho=g("rho"))
    s["vap"] = sum(v for v in (g(f) for f in VAP) if v is not None)
    rhoM = g("rhoM")
    if rhoM is None:
        cs = [read_field(f, n) for f in glob.glob(os.path.join(td, "c.*")) if re.fullmatch(r"c\.\w+", os.path.basename(f))]
        rhoM = np.maximum(sum(cs), 1e-3)
    s["mis"] = s["rho"]/rhoM - 1
    return s


def grid(s, q):
    xs, ys = np.unique(s["x"]), np.unique(s["y"])
    ix, iy = np.searchsorted(xs, s["x"]), np.searchsorted(ys, s["y"])
    a = np.full((len(ys), len(xs)), np.nan)
    a[iy, ix] = s[q]
    return xs, ys, a


def loglines(run):
    rows, t = [], None
    for l in open(os.path.join(run, "log.solver"), errors="replace"):
        if l.startswith("Time = "):
            t = float(l.split()[2])
        m = re.match(r"momentumDensity \w+: max \|rho/rhoM - 1\| = (\S+), max \|ddt", l)
        if m:
            rows.append((t, float(m.group(1))))
    txt = open(os.path.join(run, "log.solver"), errors="replace").read()
    prob = re.search(r"PROBLEM DETECTED at stage '([^']+)', time ([-\d.eE+]+)", txt)
    return np.array(rows), prob


# ------------------------------------------------------------------ fields
QS = [("T", "T [K]", "inferno", None), ("p", "log10 p [Pa]", "viridis", None),
      ("U", "|U| [m/s]", "magma", None), ("mis", "rho/rhoM - 1", "RdBu_r", 0.3)]
for label, run in RUNS:
    have = [d for d in TS if os.path.isfile(os.path.join(run, d, "C"))]
    if not have:
        continue
    fig, axs = plt.subplots(len(have), 4, figsize=(19, 3.6*len(have)), squeeze=False)
    states = [load(run, d) for d in have]
    # zoom to the disturbed region: 1.25 x the largest shock extent over the plotted times
    lim = 10e-6 + 1.25*max(max(np.abs(s["x"][s["p"] > 2e5]).max(), np.abs(s["y"][s["p"] > 2e5]).max())
                           for s in states)
    for r, (d, s) in enumerate(zip(have, states)):
        for c, (q, lab, cm, sym) in enumerate(QS):
            xs, ys, a = grid(s, q)
            if q == "p":
                a = np.log10(a)
            kw = dict(cmap=cm, shading="auto")
            if sym:
                kw.update(vmin=-sym, vmax=sym)
            ax = axs[r, c]
            im = ax.pcolormesh(xs*1e6, ys*1e6, a, **kw)
            # vapour boundary (sum alpha vapour = 0.5) and shock (p = 2e5)
            _, _, va = grid(s, "vap"); _, _, pa = grid(s, "p")
            ax.contour(xs*1e6, ys*1e6, va, [0.5], colors="w" if q != "mis" else "k", linewidths=0.8)
            ax.contour(xs*1e6, ys*1e6, pa, [2e5], colors="c", linewidths=0.6, linestyles="--")
            ax.set_xlim(-lim*1e6, lim*1e6); ax.set_ylim(-lim*1e6, 0)
            ax.set_aspect("equal"); fig.colorbar(im, ax=ax, shrink=0.85, label=lab)
            ax.set_title("%s  t = %.0f ns" % (lab, float(d)*1e9), fontsize=9)
            if r == len(have) - 1:
                ax.set_xlabel("x [um]")
            if c == 0:
                ax.set_ylabel("y [um] (wall at 0)")
    fig.suptitle("%s: solid = vapour fraction 0.5, dashed = p 2e5 Pa" % label)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fields_%s.png" % label), dpi=95); plt.close(fig)

# ------------------------------------------------------------------ fronts vs time
fr = {}
lines = ["label,time,plume_x_wall,shock_x_wall,plume_h_axis,shock_h_axis,max_U,max_T,max_abs_mismatch"]
for label, run in RUNS:
    rows = []
    for t, d in times_of(run):
        if t <= 0 or not os.path.isfile(os.path.join(run, d, "C")):
            continue
        s = load(run, d)
        dy = np.diff(np.unique(s["y"]))[0]
        wall = s["y"] > -dy          # first gas row
        axis = np.abs(s["x"]) < dy   # the two columns next to x = 0
        f = lambda m, q: q[m].max() if m.any() else np.nan
        row = (t, f(wall & (s["vap"] > 0.01), np.abs(s["x"])), f(wall & (s["p"] > 2e5), np.abs(s["x"])),
               f(axis & (s["vap"] > 0.01), -s["y"]), f(axis & (s["p"] > 2e5), -s["y"]),
               s["U"].max(), s["T"].max(), np.abs(s["mis"]).max())
        rows.append(row)
        lines.append(label + "," + ",".join("%.6g" % v for v in row))
    fr[label] = (np.array(rows), *loglines(run))
open(os.path.join(OUT, "fronts.csv"), "w", newline="\n").write("\n".join(lines) + "\n")

fig, axs = plt.subplots(2, 3, figsize=(17, 9))
cols = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#17becf", "#7f7f7f"]
for (label, (a, ll, prob)), col in zip(fr.items(), cols):
    if len(a) == 0:
        continue
    t = a[:, 0]*1e9
    st = " (stop %.2f ns)" % (float(prob.group(2))*1e9) if prob else ""
    axs[0, 0].plot(t, a[:, 1]*1e6, "o-", color=col, ms=3, label="%s vapour%s" % (label, st))
    axs[0, 0].plot(t, a[:, 2]*1e6, "s--", color=col, ms=3, mfc="none", label="%s shock" % label)
    axs[0, 1].plot(t, a[:, 3]*1e6, "o-", color=col, ms=3, label="%s vapour" % label)
    axs[0, 1].plot(t, a[:, 4]*1e6, "s--", color=col, ms=3, mfc="none", label="%s shock" % label)
    axs[0, 2].plot(t, a[:, 5], "o-", color=col, ms=3, label=label)
    axs[1, 0].plot(t, a[:, 6], "o-", color=col, ms=3, label=label)
    if len(ll):
        axs[1, 1].semilogy(ll[:, 0]*1e9, ll[:, 1], "-", color=col, lw=0.8, label=label)
    # lateral front speed along the wall
    if len(a) > 2:
        v = np.gradient(a[:, 1], a[:, 0])
        axs[1, 2].plot(t, v, "o-", color=col, ms=3, label="%s vapour" % label)
        axs[1, 2].plot(t, np.gradient(a[:, 2], a[:, 0]), "s--", color=col, ms=3, mfc="none", label="%s shock" % label)
tt = np.linspace(1, max((a[-1, 0]*1e9 for a, _, _ in fr.values() if len(a)), default=20), 100)
for ax in (axs[0, 0], axs[0, 1]):
    ax.plot(tt, (2*0.212/1.6)**0.25*np.sqrt(tt*1e-9)*1e6, "k:", lw=1, label="half-cylinder blast, E_l 0.21 J/m")
axs[0, 0].set_ylabel("lateral front along the wall [um] (first gas row)")
axs[0, 1].set_ylabel("height on the axis [um]")
axs[0, 2].set_ylabel("max |U| [m/s]"); axs[1, 0].set_ylabel("max T [K]")
axs[1, 1].set_ylabel("solver line: max |rho/rhoM - 1| per step")
axs[1, 2].set_ylabel("lateral front speed [m/s]")
for ax in axs.flat:
    ax.set_xlabel("t [ns]"); ax.grid(alpha=0.3); ax.legend(fontsize=7)
fig.suptitle("2D plume edge (r0 30 um, 1000 bar, 12 kK Fe vapour into argon 1 bar): fronts vs time")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fronts.png"), dpi=110); plt.close(fig)
print("\n".join(lines))
