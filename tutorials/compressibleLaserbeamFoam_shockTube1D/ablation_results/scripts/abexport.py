"""Small CSVs behind the ablation tables and plots.

    python abexport.py <outdir>

  csv/series_<case>_<V>.csv        per time step, from log.solver
  csv/contact_<case>_<V>.csv     1D, per write: contact position/velocity, max |U|, min T
  csv/fields_P_<V>.csv           2D, per write: max |U|, min T, max |rho/rhoM - 1| from the written fields
  csv/profile_<case>_<V>.csv     1D profile at 60 ns (B V0: 26 ns, its last write before the stop); A/Af with exact
"""
import os, re, sys
import numpy as np

H = r"\\wsl.localhost\Ubuntu\home\g05228sr"
OUT = os.path.join(sys.argv[1], "csv")
os.makedirs(OUT, exist_ok=True)
M = H + r"\shockTube_mm\runs"
P = H + r"\plumeEdge2D\runs"
VAR = [("V0", "full"), ("V1", "full_mm"), ("V2", "full_nc+mm"), ("V3", "full_o10+nc+mm"), ("V4", "contErr_o10+nc+mm")]
sys.path.insert(0, H + r"\shockTube_v5")
sys.path.insert(0, H + r"\laserbeamfoam\tutorials\tools")
sys.argv = sys.argv[:1]
import summarise, riemann
from plumeFront import read_field, times_of, VAP


def write(name, head, rows):
    with open(os.path.join(OUT, name), "w", newline="\n") as f:
        f.write(head + "\n")
        for r in rows:
            f.write(",".join("" if v is None else (v if isinstance(v, str) else "%.8g" % v) for v in r) + "\n")


def stop_time(run):
    txt = open(os.path.join(run, "log.solver"), errors="replace").read()
    m = re.search(r"PROBLEM DETECTED at stage '[^']+', time ([-\d.eE+]+)", txt)
    return float(m.group(1)) if m else None


def log_series(run):
    """One row per step: last max(U)/min(T) of the step, largest max(U) over its correctors, the momentumDensity line."""
    rows, t, dt, us, ts, mis, cont = [], None, None, [], [], None, None

    def flush():
        if t is not None and us:
            rows.append((t, dt, us[-1], max(us), ts[-1] if ts else None, mis, cont))
    for l in open(os.path.join(run, "log.solver"), errors="replace"):
        if l.startswith("Time = "):
            flush()
            t, us, ts, mis, cont = float(l.split()[2]), [], [], None, None
        elif l.startswith("deltaT = "):
            dt = float(l.split()[2])
        elif l.startswith("max(U) "):
            try:
                us.append(float(l.split()[1]))
            except (ValueError, IndexError):
                pass
        elif l.startswith("min(T) "):
            try:
                ts.append(float(l.split()[1]))
            except (ValueError, IndexError):
                pass
        else:
            m = re.match(r"momentumDensity \w+: max \|rho/rhoM - 1\| = (\S+), max \|ddt\(rhoM\) \+ div\(rhoPhi\)\| = (\S+)", l)
            if m:
                mis, cont = float(m.group(1)), float(m.group(2))
    flush()
    return rows


LOGH = "t_s,deltaT_s,maxU_final_m_s,maxU_step_m_s,minT_K,max_abs_rho_over_rhoM_minus1,max_contErr_kg_m3_s"
for k, s in VAR:
    for c, base in (("A", M), ("Af", M), ("B", M), ("P", P)):
        run = os.path.join(base, "%s_%s" % (c, s))
        rows = log_series(run)
        write("series_%s_%s.csv" % (c, k), LOGH, rows)
        print(c, k, len(rows), "steps; last t", rows[-1][0], "maxU", rows[-1][3], "minT", rows[-1][4])


def cross(x, a, y):
    """Interpolate y at the first alpha.air = 0.5 crossing."""
    i = next((j for j in range(len(x) - 1) if a[j] < 0.5 <= a[j + 1]), None)
    if i is None:
        return None, None
    w = (0.5 - a[i])/(a[i + 1] - a[i])
    return x[i] + w*(x[i + 1] - x[i]), y[i] + w*(y[i + 1] - y[i])


for c in ("A", "Af", "B"):
    for k, s in VAR:
        run = os.path.join(M, "%s_%s" % (c, s))
        pp = os.path.join(run, "postProcessing", "sampleDict")
        st = stop_time(run)
        rows = []
        for tn in sorted(os.listdir(pp), key=float):
            t = float(tn)
            if t == 0:
                continue
            d = summarise.read_sets(os.path.join(pp, tn))
            x, a, U, T, p = (np.array(d[q], float) for q in ("x", "alpha.air", "U", "T", "p"))
            xc, Uc = cross(x, a, U)
            _, Tc = cross(x, a, T)
            vm = a < 0.5
            rows.append([t, xc, Uc, None, np.abs(U).max(), T.min(), T[vm].min() if vm.any() else None, Tc, p.max(), "health dump at stop" if st and t >= st - 1e-15 else ""])
        # contact speed from the written positions (central differences, one-sided at the ends)
        ts = np.array([r[0] for r in rows]); xs = np.array([r[1] if r[1] is not None else np.nan for r in rows])
        v = np.gradient(xs, ts)
        for r, q in zip(rows, v):
            r[3] = q
        write("contact_%s_%s.csv" % (c, k), "t_s,x_contact_m,U_contact_m_s,dx_contact_dt_m_s,maxU_m_s,minT_K,minT_vapour_K,T_contact_K,maxp_Pa,note", rows)
        # 60 ns profile
        tn = "6e-08" if "6e-08" in os.listdir(pp) else "2.6e-08"
        d = summarise.read_sets(os.path.join(pp, tn))
        x = list(d["x"])
        cols = [("p_Pa", "p"), ("U_m_s", "U"), ("T_K", "T"), ("rho_kg_m3", "rho"), ("alpha_air", "alpha.air"),
                ("c_air", "c.air"), ("c_vapour", "c.metal1vapour")]
        head = "x_m," + ",".join(n for n, _ in cols)
        pr = [[x[i]] + [d[q][i] for _, q in cols] for i in range(len(x))]
        if c != "B":
            sol = riemann.solution("A")
            head += ",p_exact_Pa,U_exact_m_s,T_exact_K,rho_exact_kg_m3"
            for r, xx in zip(pr, x):
                e = riemann.sample(sol, xx, 6e-8)
                r += [e[2], e[1], e[3], e[0]]
        write("profile_%s_%s_t%s.csv" % (c, k, tn), head, pr)

# 2D, per write, from the fields
for k, s in VAR:
    run = os.path.join(P, "P_%s" % s)
    st = stop_time(run)
    rows = []
    for t, d in times_of(run):
        td = os.path.join(run, d)
        if t <= 0 or not os.path.isfile(os.path.join(td, "C")):
            continue
        n = len(read_field(os.path.join(td, "C")))
        g = lambda f: read_field(os.path.join(td, f), n)
        rhoM = g("rhoM")
        if rhoM is None:
            rhoM = np.maximum(g("c.air") + g("c.metal1vapour"), 1e-3)
        vap = sum(v for v in (g(f) for f in VAP) if v is not None)
        rows.append((t, np.sqrt((g("U")**2).sum(1)).max(), g("T").min(), g("T")[vap > 0.5].min() if (vap > 0.5).any() else None, g("T").max(), g("p").max(),
                     np.abs(g("rho")/rhoM - 1).max(), "health dump at stop" if st and t >= st - 1e-15 else ""))
    write("fields_P_%s.csv" % k, "t_s,maxU_m_s,minT_K,minT_vapour_K,maxT_K,maxp_Pa,max_abs_rho_over_rhoM_minus1,note", rows)
    print("P", k, "field writes", len(rows), "last", rows[-1][:3])
