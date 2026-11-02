"""Ablation 1D profiles at 60 ns (A, Af vs exact; B) for V0-V4."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = r"\\wsl.localhost\Ubuntu\home\g05228sr"
OUT = sys.argv[1]
M = H + r"\shockTube_mm\runs"
sys.path.insert(0, H + r"\shockTube_v5"); sys.argv = sys.argv[:1]
import summarise, riemann
VAR = [("V0", "full"), ("V1", "full_mm"), ("V2", "full_nc+mm"), ("V3", "full_o10+nc+mm"), ("V4", "contErr_o10+nc+mm")]
COL = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
for c in ("A", "Af", "B"):
    fig, axs = plt.subplots(2, 2, figsize=(14, 8), sharex=True)
    if c != "B":
        sol = riemann.solution("A"); xe = np.linspace(-200e-6, 400e-6, 2001)
        ex = np.array([(x,) + riemann.sample(sol, x, 6e-8) for x in xe])
        for ax, ec in zip(axs.flat, (3, 2, 4, 1)):
            ax.plot(ex[:, 0]*1e6, ex[:, ec], "k-", lw=1.2, label="exact")
    for (k, s), col in zip(VAR, COL):
        pp = os.path.join(M, "%s_%s" % (c, s), "postProcessing", "sampleDict")
        ts = sorted(os.listdir(pp), key=float)
        tn = "6e-08" if "6e-08" in ts else None
        if c == "B" and k == "V0":
            tn = "2.6e-08"
        if tn is None or not os.path.isdir(os.path.join(pp, tn)):
            continue
        d = summarise.read_sets(os.path.join(pp, tn))
        ls = ":" if k == "V3" else "-"
        for ax, q in zip(axs.flat, ("p", "U", "T", "rho")):
            if q in d:
                y = d[q] if q != "U" else d["U"]
                ax.plot(np.array(d["x"])*1e6, y, ls, color=col, lw=1.3 if k != "V3" else 2.2,
                        label="%s%s" % (k, " (26 ns, stopped)" if tn != "6e-08" else ""))
    for ax, lab in zip(axs.flat, ("p [Pa]", "U_x [m/s]", "T [K]", "rho [kg/m3]")):
        ax.set_ylabel(lab); ax.grid(alpha=0.3)
    for ax in axs[1]:
        ax.set_xlabel("x [um]")
    axs[0, 0].legend(fontsize=8)
    fig.suptitle("Case %s at 60 ns: ablation V0 full / V1 +mm / V2 +nc / V3 +o10 / V4 contErr (V3 dotted under V4)" % c)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "profiles_%s.png" % c), dpi=110); plt.close(fig)
print("ok", list(summarise.read_sets(os.path.join(M, "A_full", "postProcessing", "sampleDict", "6e-08")).keys()))
