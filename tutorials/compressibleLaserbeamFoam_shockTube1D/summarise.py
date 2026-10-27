#!/usr/bin/env python3
"""Summarise shock-tube runs:  python3 summarise.py runs/*

For each run: how far it got, whether the health check stopped it, and at
each sampled time the shock/contact positions and plateau values. Case A
runs are compared with the exact Riemann solution (riemann.py).
"""
import sys, os, glob, re, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import riemann

def read_sets(d):
    data = {}
    for fn in glob.glob(os.path.join(d, "line_*.xy")):
        names = os.path.basename(fn)[5:-3].split("_")
        rows = [l.split() for l in open(fn) if l.strip() and not l.startswith("#")]
        if not rows: continue
        x = [float(r[0]) for r in rows]
        data["x"] = x
        if names == ["U"]:
            data["U"] = [float(r[1]) for r in rows]          # x-component
        else:
            for k, n in enumerate(names):
                data[n] = [float(r[k + 1]) for r in rows]
    return data

def metrics(d):
    x, p, a = d["x"], d["p"], d.get("alpha.air")
    ish = max((i for i in range(len(x)) if p[i] > 1.5e5 and x[i] > 0), default=None)
    icon = next((i for i in range(len(x)) if a and a[i] > 0.5), None)
    m = dict(maxT=max(d["T"]), minp=min(p), maxU=max(abs(u) for u in d.get("U", [0])))
    if ish is not None: m["shock"] = x[ish]
    if icon is not None: m["contact"] = x[icon]
    if ish is not None and icon is not None and ish - icon > 4:
        mid = range(icon + 2, ish - 1)
        m["pStar"] = statistics.median(p[i] for i in mid)
        m["Tair"] = statistics.median(d["T"][i] for i in mid)
        if "U" in d: m["uStar"] = statistics.median(d["U"][i] for i in mid)
        if "rho" in d: m["rhoAir"] = statistics.median(d["rho"][i] for i in mid)
        lv = range(max(icon - 6, 0), max(icon - 1, 0))
        if len(lv): m["Tvap"] = statistics.median(d["T"][i] for i in lv)
    return m

def fmt(v, f="%.4g"):
    return "-" if v is None else f % v

for run in sys.argv[1:]:
    name = os.path.basename(run.rstrip("/"))
    case = name.split("_")[0].rstrip("f")
    log = os.path.join(run, "log.solver")
    txt = open(log, errors="replace").read() if os.path.exists(log) else ""
    times = re.findall(r"^Time = (\S+)", txt, re.M)
    mode = re.search(r"kineticEnergyTerms (\w+)", txt)
    prob = re.search(r"PROBLEM DETECTED at stage '([^']+)', time (\S+)", txt)
    fatal = "FOAM FATAL" in txt and not prob
    print("=" * 78)
    print("%s   mode in log: %s   last time: %s   %s" % (
        name, mode.group(1) if mode else "NOT REPORTED (old binary?)",
        times[-1] if times else "-",
        ("STOPPED: %s at %s" % prob.groups()) if prob else ("FATAL ERROR, see log" if fatal else
         ("finished" if "End" in txt.split("\n")[-20:] or re.search(r"^End\s*$", txt, re.M) else "running/incomplete"))))
    ex = riemann.solution(case) if case == "A" else None
    print("  %-9s %-10s %-10s %-10s %-9s %-8s %-8s %-9s %-9s %-9s" %
          ("t [ns]", "shock[um]", "contact", "p* [Pa]", "u* [m/s]", "T air", "T vap", "max T", "min p", "max|U|"))
    sd = sorted(glob.glob(os.path.join(run, "postProcessing", "sampleDict", "*")),
                key=lambda s: float(os.path.basename(s)))
    for tdir in sd:
        t = float(os.path.basename(tdir))
        if t <= 0: continue
        d = read_sets(tdir)
        if "p" not in d: continue
        m = metrics(d)
        print("  %-9s %-10s %-10s %-10s %-9s %-8s %-8s %-9s %-9s %-9s" % (
            "%.1f" % (t*1e9), fmt(m.get("shock") and m["shock"]*1e6, "%.1f"),
            fmt(m.get("contact") and m["contact"]*1e6, "%.1f"), fmt(m.get("pStar")),
            fmt(m.get("uStar"), "%.0f"), fmt(m.get("Tair"), "%.0f"), fmt(m.get("Tvap"), "%.0f"),
            fmt(m["maxT"], "%.0f"), fmt(m["minp"]), fmt(m["maxU"], "%.0f")))
        if ex:
            print("  %-9s %-10s %-10s %-10s %-9s %-8s %-8s" % (
                "  exact", "%.1f" % (ex["S"]*t*1e6), "%.1f" % (ex["us"]*t*1e6),
                "%.4g" % ex["ps"], "%.0f" % ex["us"], "%.0f" % ex["TRs"], "%.0f" % ex["TLs"]))
