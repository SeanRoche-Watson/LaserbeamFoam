"""Ablation V0-V4 (1D A/Af/B, 2D plume edge): tables.

    python ablation.py <outdir>
Needs <outdir>/energy_table.txt (energy_check.py on links A_V0 ...) and <outdir>/2d/fronts.csv (plumeEdge.py).
"""
import os, re, sys, statistics
import numpy as np

H = r"\\wsl.localhost\Ubuntu\home\g05228sr"
OUT = sys.argv[1]
M = H + r"\shockTube_mm\runs"
P = H + r"\plumeEdge2D\runs"
VAR = [("V0", "full", ""), ("V1", "full", "mm"), ("V2", "full", "nc+mm"), ("V3", "full", "o10+nc+mm"),
       ("V4", "contErr", "o10+nc+mm")]
DESC = {"V0": "full (eos, cap on, resCtrl)", "V1": "full:mm", "V2": "full:nc+mm", "V3": "full:o10+nc+mm",
        "V4": "contErr:o10+nc+mm"}
run = lambda base, c, mode, v: os.path.join(base, "%s_%s%s" % (c, mode, "_" + v if v else ""))


def loginfo(r, tend):
    txt = open(os.path.join(r, "log.solver"), errors="replace").read()
    prob = re.search(r"PROBLEM DETECTED at stage '([^']+)', time ([-\d.eE+]+)", txt)
    ok = re.search(r"^End\s*$", txt, re.M)
    status = ("STOP %.2f ns" % (float(prob.group(2))*1e9) if prob else ("%g ns" % (tend*1e9) if ok else "incomplete"))
    ct = re.findall(r"ClockTime = (\d+)", txt)
    rows, t = [], None
    for l in txt.split("\n"):
        if l.startswith("Time = "):
            t = float(l.split()[2])
        m = re.match(r"momentumDensity \w+: max \|rho/rhoM - 1\| = (\S+), max \|ddt\(rhoM\) \+ div\(rhoPhi\)\| = (\S+)", l)
        if m:
            rows.append((t, float(m.group(1)), float(m.group(2))))
    a = np.array(rows)
    if prob and len(a):
        a = a[a[:, 0] <= float(prob.group(2))]
    nsteps = len(re.findall(r"^deltaT = ", txt, re.M))
    nouter = len(re.findall(r"^PIMPLE: iteration", txt, re.M))
    return dict(status=status, wall=int(ct[-1]) if ct else None, mis=a[:, 1].max() if len(a) else np.nan,
                mis_med=np.median(a[:, 1]) if len(a) else np.nan, cont=a[:, 2].max() if len(a) else np.nan,
                steps=nsteps, outer=nouter/max(nsteps, 1), prob=prob)


# energy drift: Table 1 "last good" column
drift = {}
et = open(os.path.join(OUT, "energy_table.txt")).read()
t1 = et[et.index("Table 1:"):et.index("Table 2")]
for l in t1.split("\n"):
    m = re.match(r"(\w+_V\d)\s+(.*?)\s{2,}.*?([+-]\d+\.\d+) \(\s*([\d.]+)\)", l)
    if m:
        drift[m.group(1)] = (float(m.group(3)), float(m.group(4)))

sys.path.insert(0, H + r"\shockTube_v5"); sys.argv = sys.argv[:1]
import summarise

L = ["Ablation V0-V4, cc-mass-momentum (33e3963f build), writeInterval 1e-9 (1D) / 1e-9 (2D)", ""]
for k, _, _ in VAR:
    L.append("  %s  %s" % (k, DESC[k]))
L += ["  resCtrl = residualControl on (default, ~2 outer iterations); o10 = 10 outer correctors, residualControl off;",
      "  nc = gasPackingCap false; mm = momentumDensity mass. Wall time: serial 1D runs, 4 at a time on 12 cores",
      "  (V4 and B V0/V1 reused from 27 Oct, run 2-4 at a time); 2D: 2 ranks each, 4 at a time.", ""]
L.append("max mis = max over steps of the solver line max |rho/rhoM - 1|; max cont = max |ddt(rhoM)+div(rhoPhi)| [kg/m3/s]")
L.append("(both up to any stop); drift = E_tot drift with transported masses (energy_check.py Table 1) at the last good time.")
L.append("")

f = lambda v, s: "-" if v is None or (isinstance(v, float) and np.isnan(v)) else s % v
EX = ("exact", "-84.5", "4579", "2.48e6", "86.6", "2.48e6", "1029", "2138", "5.57")
for c in ("A", "Af", "B"):
    L.append("=== Case %s (60 ns)" % c)
    h = "%-4s %-20s %-12s %6s %6s %7s %8s %9s %8s | %7s %7s %8s %8s %8s %7s %6s %7s"
    L.append(h % ("", "settings", "status", "wall s", "steps", "outer/s", "max mis", "max cont", "drift %",
                  "head", "T vap", "p vap", "shock", "p*", "u*", "T air", "rho air"))
    L.append("-"*160)
    if c != "B":
        L.append(h % (("", "exact", "", "", "", "", "", "", "") + EX[1:]))
    for k, mode, v in VAR:
        r = run(M, c, mode, v)
        if not os.path.isdir(r):
            L.append("%-4s missing %s" % (k, r)); continue
        li = loginfo(r, 6e-8)
        d = drift.get("%s_%s" % (c, k))
        ds = "%+.3f" % d[0] + ("" if d[1] >= 59.9 else "@%.0f" % d[1]) if d else "-"
        pp = os.path.join(r, "postProcessing", "sampleDict")
        ts = sorted(os.listdir(pp), key=float) if os.path.isdir(pp) else []
        if li["prob"]:
            ts = [t for t in ts if float(t) < float(li["prob"].group(2)) - 1e-15]
        cols = ["-"]*8
        if ts:
            tn = "6e-08" if "6e-08" in ts else ts[-1]
            dd = summarise.read_sets(os.path.join(pp, tn)); m = summarise.metrics(dd)
            x, p, T, a = dd["x"], dd["p"], dd["T"], dd["alpha.air"]
            p0 = 0.99e7 if c != "B" else 0.99e8
            head = min((x[i] for i in range(len(x)) if p[i] < p0), default=None)
            vap = [i for i in range(len(x)) if a[i] < 0.01 and -10e-6 < x[i] < m["contact"] - 5e-6]
            Tvp = statistics.median(T[i] for i in vap) if vap else None
            pvp = statistics.median(p[i] for i in vap) if vap else None
            cols = [f(head and head*1e6, "%.1f"), f(Tvp, "%.0f"), f(pvp, "%.4g"), f(m.get("shock") and m["shock"]*1e6, "%.1f"),
                    f(m.get("pStar"), "%.4g"), f(m.get("uStar"), "%.0f"), f(m.get("Tair"), "%.0f"), f(m.get("rhoAir"), "%.3f")]
            if tn != "6e-08":
                cols[0] = "@%.0fns %s" % (float(tn)*1e9, cols[0])
        L.append(h % (k, DESC[k][:20], li["status"], f(li["wall"], "%d"), li["steps"], "%.1f" % li["outer"], f(li["mis"], "%.3f"),
                      f(li["cont"], "%.2e"), ds, *cols))
    L.append("")

# ------------------------------------------------------------ V3 vs V4
L.append("=== V3 (full:o10+nc+mm) vs V4 (contErr:o10+nc+mm): max relative difference over all written times")
FLD = ["p", "U", "T", "rho", "c.air", "c.metal1vapour", "alpha.air"]


def rd(fn):
    txt = open(fn).read()
    s = txt[txt.index("internalField") + 13:]
    m = re.match(r"\s*uniform\s+(\([^)]*\)|\S+?);", s)
    if m:
        return np.array([float(x) for x in m.group(1).strip("()").split()])
    m = re.match(r"\s*nonuniform\s+List<\w+>\s*(\d+)\s*\(", s)
    body = s[m.end():]
    body = body[:body.index("\n)")]
    return np.array(body.replace("(", " ").replace(")", " ").split(), float)


for c, base in (("A", M), ("Af", M), ("B", M), ("P", P)):
    r3, r4 = run(base, c, "full", "o10+nc+mm"), run(base, c, "contErr", "o10+nc+mm")
    if not (os.path.isdir(r3) and os.path.isdir(r4)):
        L.append("  %s: missing" % c); continue
    isT = lambda d: re.fullmatch(r"[\d.e+-]+", d) and float(d) > 0
    ts = sorted(set(filter(isT, os.listdir(r3))) & set(filter(isT, os.listdir(r4))), key=float)
    worst, nsame, ntot = {}, 0, 0
    for t in ts:
        for q in FLD:
            a, b = os.path.join(r3, t, q), os.path.join(r4, t, q)
            if not (os.path.isfile(a) and os.path.isfile(b)):
                continue
            ntot += 1
            if open(a, "rb").read() == open(b, "rb").read():
                nsame += 1; continue
            x, y = rd(a), rd(b)
            e = np.abs(x - y).max()/max(np.abs(x).max(), 1e-300)
            if e > worst.get(q, (0,))[0]:
                worst[q] = (e, t)
    L.append("  %-3s %d times, %d/%d field files byte-identical; worst: %s" % (
        c, len(ts), nsame, ntot, ", ".join("%s %.1e (%.0f ns)" % (q, e, float(t)*1e9) for q, (e, t) in worst.items()) or "none"))
L.append("")

# ------------------------------------------------------------ 2D
L.append("=== 2D plume edge P (160x100, 0-20 ns): fronts along the wall (first gas row) [um], vapour (sum alpha_vap > 0.01) / shock (p > 2e5)")
fr = {}
fc = os.path.join(OUT, "2d", "fronts.csv")
if os.path.isfile(fc):
    for l in open(fc).read().split("\n")[1:]:
        if l:
            s = l.split(","); fr.setdefault(s[0], {})[round(float(s[1])*1e9)] = [float(v) for v in s[2:]]
h2 = "%-4s %-20s %-12s %6s %6s %7s %8s %9s | %-13s %-13s %-13s | %7s %7s"
L.append(h2 % ("", "settings", "status", "wall s", "steps", "outer/s", "max mis", "max cont", "5 ns vap/shk", "10 ns vap/shk",
               "20 ns vap/shk", "maxU 20", "mis@20"))
L.append("-"*150)
for k, mode, v in VAR:
    r = run(P, "P", mode, v)
    if not os.path.isdir(r):
        L.append("%-4s missing" % k); continue
    li = loginfo(r, 2e-8)
    g = fr.get(k, {})
    cell = lambda t: "%5.1f / %5.1f" % (g[t][0]*1e6, g[t][1]*1e6) if t in g else "-"
    L.append(h2 % (k, DESC[k][:20], li["status"], f(li["wall"], "%d"), li["steps"], "%.1f" % li["outer"], f(li["mis"], "%.3f"),
                   f(li["cont"], "%.2e"), cell(5), cell(10), cell(20),
                   "%.0f" % g[20][4] if 20 in g else "-", "%.3f" % g[20][6] if 20 in g else "-"))
L.append("")
txt = "\n".join(L)
print(txt)
open(os.path.join(OUT, "ablation_tables.txt"), "w", newline="\n").write(txt + "\n")
