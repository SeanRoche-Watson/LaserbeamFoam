#!/usr/bin/env python3
"""Mass and energy conservation check for the closed 1D shock tube.

    python3 energy_check.py runs/*          (run from the shockTube case directory)

The tube is closed (wall ends, U = 0, zeroGradient T), laser power 0, no liquid,
so total mass and total energy must stay constant. For every written time of
every run this reads the ASCII cell fields T, U, rho, c.air, c.metal1vapour,
alpha.air, alpha.metal1vapour and computes

    M_i     = sum(c_i V)
    E_int   = sum((c_air e_air(T) + c_vap e_vap(T)) V)
    E_kin   = sum((c_air + c_vap) |U|^2/2 V)
    E_tot   = E_int + E_kin
    E_EOS   = sum(rho (Y_air e_air + Y_vap e_vap + |U|^2/2) V),  Y_i = c_i/sum(c)

with e_i(T) = int_298.15^T Cv_i dT', Cv_i = Cp_i - R_i, Cp_i from the janaf
coefficients in constant/thermophysicalProperties.<phase> (integrated exactly).

t = 0: the solver writes no rho/c fields at the start time, so c_i is rebuilt
from 0/ (alpha_i, p, T) with the perfect-gas EOS, c_i = alpha_i p/(R_i T),
which is how the solver initialises them; this is checked against the analytic
two-state value.

Writes energy_table.txt and energy_<run>.csv (one row per time) next to runs/.
"""
import os, re, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # riemann.py

RU = 8314.47            # J/(kmol K)
TSTD = 298.15
PHASES = {"air": "air", "vap": "metal1vapour"}
CASES = {"A": (1e7, 8000.0), "B": (1e8, 12000.0)}   # left pL, TL; right 1e5 Pa, 300 K
PR, TR = 1e5, 300.0
XMIN, XMAX, XDIA = -200e-6, 400e-6, 0.0
CHECK_TIMES = [10e-9, 20e-9, 30e-9, 40e-9, 50e-9, 60e-9]


# ------------------------------------------------------------------ janaf ---
class Janaf:
    def __init__(self, fn):
        txt = open(fn).read()
        num = lambda k: float(re.search(r"\b%s\s+([-\d.eE+]+)\s*;" % k, txt).group(1))
        coeffs = lambda k: [float(v) for v in re.search(r"\b%s\s*\(([^)]*)\)" % k, txt).group(1).split()]
        self.W = num("molWeight")
        self.R = RU/self.W
        self.Tc = num("Tcommon")
        self.lo = coeffs("lowCpCoeffs")
        self.hi = coeffs("highCpCoeffs")

    @staticmethod
    def cpR(a, T):
        return a[0] + T*(a[1] + T*(a[2] + T*(a[3] + T*a[4])))

    @staticmethod
    def intCpR(a, T):        # antiderivative of Cp/R
        return T*(a[0] + T*(a[1]/2 + T*(a[2]/3 + T*(a[3]/4 + T*a[4]/5))))

    def Cv(self, T):
        a = self.lo if T < self.Tc else self.hi
        return self.R*(self.cpR(a, T) - 1.0)

    def e(self, T):
        """int_TSTD^T Cv dT' [J/kg], piecewise exact (assumes TSTD < Tcommon)."""
        lo, hi = self.lo, self.hi
        if T <= self.Tc:
            hR = self.intCpR(lo, T) - self.intCpR(lo, TSTD)
        else:
            hR = (self.intCpR(lo, self.Tc) - self.intCpR(lo, TSTD)
                  + self.intCpR(hi, T) - self.intCpR(hi, self.Tc))
        return self.R*(hR - (T - TSTD))

    def checks(self):
        Tc = self.Tc
        out = {}
        out["Cp/R low,high at Tcommon"] = (self.cpR(self.lo, Tc), self.cpR(self.hi, Tc))
        # e(T) evaluated purely with each branch at Tcommon: must coincide
        e_lo = self.R*(self.intCpR(self.lo, Tc) - self.intCpR(self.lo, TSTD) - (Tc - TSTD))
        eps = 1e-6*Tc
        out["e(Tc-eps), e(Tc), e(Tc+eps)"] = (self.e(Tc - eps), e_lo, self.e(Tc + eps))
        # janaf absolute enthalpy (with a6) continuity, as OpenFOAM evaluates Ha
        ha = lambda a: self.R*(self.intCpR(a, Tc) + a[5])
        out["Ha/R low,high at Tcommon (with a6) [J/kg]"] = (ha(self.lo), ha(self.hi))
        # quadrature cross-check of e(T) (composite Simpson, 20000 intervals)
        for T in (8000.0, 12000.0, 20000.0):
            n = 20000; h = (T - TSTD)/n
            s = self.Cv(TSTD) + self.Cv(T)
            s += sum((4 if k % 2 else 2)*self.Cv(TSTD + k*h) for k in range(1, n))
            out["e(%g) exact vs Simpson" % T] = (self.e(T), s*h/3)
        return out


# ---------------------------------------------------------------- fields ----
def read_field(fn, n, vector=False):
    txt = open(fn).read()
    i = txt.index("internalField")
    s = txt[i + len("internalField"):]
    m = re.match(r"\s*uniform\s+(\([^)]*\)|[-\d.eE+]+)", s)
    if m:
        v = m.group(1)
        v = [float(x) for x in v.strip("()").split()] if vector else float(v)
        return [v]*n
    m = re.match(r"\s*nonuniform\s+List<\w+>\s*(\d+)\s*", s)
    if not m:
        raise ValueError("cannot parse internalField in " + fn)
    cnt = int(m.group(1)); s = s[m.end():]
    if s.startswith("{"):        # compact form N{value}
        v = re.match(r"\{\s*(\([^)]*\)|[-\d.eE+]+)\s*\}", s).group(1)
        v = [float(x) for x in v.strip("()").split()] if vector else float(v)
        return [v]*cnt
    body = s[s.index("(") + 1:]
    if vector:
        vals = [[float(x) for x in t.split()] for t in re.findall(r"\(([^()]*)\)", body[:body.index(";")])[:cnt]]
    else:
        vals = [float(x) for x in body[:body.index(")")].split()]
    if len(vals) != cnt:
        raise ValueError("%s: expected %d values, got %d" % (fn, cnt, len(vals)))
    return vals


def mesh_info(run):
    """Cell count from blockMeshDict and the cell volume, checked against the points."""
    bm = open(os.path.join(run, "system", "blockMeshDict")).read()
    nx = int(re.search(r"hex\s*\([^)]*\)\s*\(\s*(\d+)\s+1\s+1\s*\)", bm).group(1))
    dx = (XMAX - XMIN)/nx
    V = dx*(3.125e-6)**2
    pts = open(os.path.join(run, "constant", "polyMesh", "points")).read()
    P = [[float(x) for x in t.split()] for t in re.findall(r"\(([-\d.eE+]+ [-\d.eE+]+ [-\d.eE+]+)\)", pts)]
    xs = sorted(set(round(p[0], 12) for p in P))
    ly = max(p[1] for p in P) - min(p[1] for p in P)
    lz = max(p[2] for p in P) - min(p[2] for p in P)
    dxs = [b - a for a, b in zip(xs, xs[1:])]
    Vmesh = (max(dxs) * ly * lz, min(dxs) * ly * lz)
    ok = (len(xs) == nx + 1 and abs(Vmesh[0] - V) < 1e-9*V and abs(Vmesh[1] - V) < 1e-9*V)
    return nx, V, ok, "points: %d x-planes, dx %.6g..%.6g um, ly %.6g um, lz %.6g um" % (
        len(xs), min(dxs)*1e6, max(dxs)*1e6, ly*1e6, lz*1e6)


def time_dirs(run):
    out = []
    for d in os.listdir(run):
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.exists(os.path.join(run, d, "T")):
            out.append((t, d))
    return sorted(out)


def integrals(run, tname, n, V, th, t0=False):
    g = lambda f, vec=False: read_field(os.path.join(run, tname, f), n, vec)
    T = g("T"); U = g("U", True)
    if t0:      # rebuild the conserved masses from the start fields
        p = g("p"); aA = g("alpha.air"); aV = g("alpha.metal1vapour")
        cA = [a*pp/(th["air"].R*TT) for a, pp, TT in zip(aA, p, T)]
        cV = [a*pp/(th["vap"].R*TT) for a, pp, TT in zip(aV, p, T)]
        rho = [a + b for a, b in zip(cA, cV)]
    else:
        cA = g("c.air"); cV = g("c.metal1vapour"); rho = g("rho")
    MA = MV = Ei = Ek = Ei_eos = Ek_eos = 0.0
    for k in range(n):
        eA = th["air"].e(T[k]); eV = th["vap"].e(T[k])
        u2 = sum(c*c for c in U[k])
        cs = cA[k] + cV[k]
        MA += cA[k]*V; MV += cV[k]*V
        Ei += (cA[k]*eA + cV[k]*eV)*V
        Ek += cs*0.5*u2*V
        YA, YV = (cA[k]/cs, cV[k]/cs) if cs > 0 else (0.0, 0.0)
        Ei_eos += rho[k]*(YA*eA + YV*eV)*V
        Ek_eos += rho[k]*0.5*u2*V
    return dict(MA=MA, MV=MV, Ei=Ei, Ek=Ek, Et=Ei + Ek, Ei_eos=Ei_eos, Ek_eos=Ek_eos,
                Et_eos=Ei_eos + Ek_eos, Msum=MA + MV, Meos=sum(rho)*V)


def analytic(case, th, nx, V):
    pL, TL = CASES[case]
    nL = round((XDIA - XMIN)/(XMAX - XMIN)*nx); nR = nx - nL
    rL = pL/(th["vap"].R*TL); rR = PR/(th["air"].R*TR)
    return dict(MV=rL*nL*V, MA=rR*nR*V, Et=(rL*th["vap"].e(TL)*nL + rR*th["air"].e(TR)*nR)*V)


def exact_split(case, t, th, nx, V, E0, s):
    """Case A only (all T < 10 kK, Cp/R = 2.5, so riemann.py is exact): exact E_kin(t)
    from the Riemann solution; E_int_exact = E0 - E_kin_exact since E_tot is conserved.
    Returns the errors (E_int - E_int_exact, E_kin - E_kin_exact) as fractions of E0."""
    if case != "A" or t <= 0:
        return None
    try:
        import riemann
    except ImportError:
        return None
    sol = riemann.solution("A")
    A = V/((XMAX - XMIN)/nx)            # cross-section
    n = 60000; dx = (XMAX - XMIN)/n
    Ek = 0.0
    for k in range(n):
        rho, u, p, T = riemann.sample(sol, XMIN + (k + 0.5)*dx, t)
        Ek += 0.5*rho*u*u*dx
    Ek *= A
    return ((s["Ei"] - (E0 - Ek))/E0, (s["Ek"] - Ek)/E0, Ek/E0)


def log_info(run):
    fn = os.path.join(run, "log.solver")
    txt = open(fn, errors="replace").read() if os.path.exists(fn) else ""
    times = re.findall(r"^Time = (\S+)", txt, re.M)
    prob = re.search(r"PROBLEM DETECTED at stage '([^']+)', time ([-\d.eE+]+)", txt)
    status = ("STOP %.2fns %s" % (float(prob.group(2))*1e9, prob.group(1).replace("after ", "")) if prob else
              ("finished" if re.search(r"^End\s*$", txt, re.M) else "incomplete"))
    drops = [float(v) for v in re.findall(r"max KE cooling dropped \(dissipative\) = (\S+) K", txt)]
    return dict(last=times[-1] if times else "-", status=status, drops=drops)


# ------------------------------------------------------------------- main ---
def main(runs):
    lines_chk, rows, csvs = [], [], {}
    th_cache = {}
    for run in runs:
        name = os.path.basename(run.rstrip("/"))
        case = name.split("_")[0].rstrip("f")
        cp = os.path.join(run, "constant")
        key = tuple(open(os.path.join(cp, "thermophysicalProperties." + PHASES[p])).read() for p in PHASES)
        if key not in th_cache:
            th_cache[key] = {p: Janaf(os.path.join(cp, "thermophysicalProperties." + PHASES[p])) for p in PHASES}
            if len(th_cache) == 1:
                lines_chk.append("janaf checks (from %s):" % name)
                for p, j in th_cache[key].items():
                    lines_chk.append("  %s (%s): molWeight %g, R %.6g J/kgK, Tcommon %g" % (p, PHASES[p], j.W, j.R, j.Tc))
                    for k, v in j.checks().items():
                        lines_chk.append("      %-44s %s" % (k, "  ".join("%.10g" % x for x in v)))
        th = th_cache[key]
        nx, V, mesh_ok, mesh_msg = mesh_info(run)
        td = time_dirs(run)
        series = []
        for t, d in td:
            series.append((t, integrals(run, d, nx, V, th, t0=(t == 0.0))))
        csvs[name] = series
        an = analytic(case, th, nx, V)
        s0 = series[0][1]
        lines_chk.append("%-22s NX %d, V %.6e m3 (mesh check %s; %s)" % (name, nx, V, "OK" if mesh_ok else "MISMATCH", mesh_msg))
        lines_chk.append("%-22s t=0: E_kin %.3e J; E_tot %.10e J vs analytic %.10e J (rel %.1e); "
                         "M_air rel %.1e, M_vap rel %.1e" % (
                             "", s0["Ek"], s0["Et"], an["Et"], (s0["Et"] - an["Et"])/an["Et"],
                             (s0["MA"] - an["MA"])/an["MA"], (s0["MV"] - an["MV"])/an["MV"]))
        rel = lambda k, s: (s[k] - s0[k])/s0[k]
        def at(tq):
            for t, s in series:
                if abs(t - tq) < 1e-12:
                    return s
            return None
        mA = max((rel("MA", s) for _, s in series), key=abs)
        mV = max((rel("MV", s) for _, s in series), key=abs)
        # relative to the first WRITTEN time: removes the 8-digit rounding offset
        # between the exact t = 0 rebuild and the writePrecision-8 files
        s1 = series[1][1] if len(series) > 1 else s0
        mA1 = max(((s["MA"] - s1["MA"])/s1["MA"] for _, s in series[1:]), key=abs, default=0.0)
        mV1 = max(((s["MV"] - s1["MV"])/s1["MV"] for _, s in series[1:]), key=abs, default=0.0)
        mEos = max(((s["Meos"] - s["Msum"])/s["Msum"] for _, s in series), key=abs)
        dE = [(t, rel("Et", s)) for t, s in series]
        dEe = [(t, (s["Et_eos"] - s0["Et_eos"])/s0["Et_eos"]) for t, s in series]
        tmax, emax = max(dE, key=lambda p: abs(p[1]))
        tmaxe, emaxe = max(dEe, key=lambda p: abs(p[1]))
        tl, sl = series[-1]
        li = log_info(run)
        # stopped run: its last write is the health-check dump of the broken state
        good = series[-2] if li["status"].startswith("STOP") and len(series) > 2 else series[-1]
        rows.append(dict(name=name, last=li["last"], status=li["status"], tlast=tl,
                         mA=mA, mV=mV, mA1=mA1, mV1=mV1, mEos=mEos,
                         tgood=good[0], dEgood=(good[1]["Et"] - s0["Et"])/s0["Et"],
                         dEgoode=(good[1]["Et_eos"] - s0["Et_eos"])/s0["Et_eos"],
                         dE={tq: (at(tq) and (at(tq)["Et"] - s0["Et"])/s0["Et"]) for tq in CHECK_TIMES},
                         dEe={tq: (at(tq) and (at(tq)["Et_eos"] - s0["Et_eos"])/s0["Et_eos"]) for tq in CHECK_TIMES},
                         emax=emax, tmax=tmax, emaxe=emaxe, tmaxe=tmaxe,
                         dEi=(good[1]["Ei"] - s0["Ei"])/s0["Et"], dEk=(good[1]["Ek"] - s0["Ek"])/s0["Et"],
                         dEie=(good[1]["Ei_eos"] - s0["Ei_eos"])/s0["Et_eos"],
                         dEke=(good[1]["Ek_eos"] - s0["Ek_eos"])/s0["Et_eos"],
                         exact=exact_split(case, good[0], th, nx, V, s0["Et"], good[1]),
                         drops=li["drops"], E0=s0["Et"]))
        with open("energy_%s.csv" % name, "w") as f:
            f.write("t,M_air,M_vap,E_int,E_kin,E_tot,E_int_EOS,E_kin_EOS,E_tot_EOS,M_EOS\n")
            for t, s in series:
                f.write("%.6e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n" % (
                    t, s["MA"], s["MV"], s["Ei"], s["Ek"], s["Et"], s["Ei_eos"], s["Ek_eos"], s["Et_eos"], s["Meos"]))

    pc = lambda v: "-" if v is None else "%+.3f" % (100*v)
    out = []
    out.append("Energy and mass conservation, closed 1D shock tube (all drifts relative to t = 0)")
    out.append("")
    h ="%-22s %-21s" % ("run", "how far") + "".join("%9s" % ("%dns" % round(t*1e9)) for t in CHECK_TIMES) \
        + "%10s %6s %15s %9s %9s" % ("max|dE|", "at ns", "last good (ns)", "dE_int", "dE_kin")
    for tab, kE, kmax, kt, kg, ki, kk in (
            ("Table 1: relative drift of E_tot = E_int + E_kin with TRANSPORTED masses c_i  [percent]",
             "dE", "emax", "tmax", "dEgood", "dEi", "dEk"),
            ("Table 2: same with the EOS density, E_tot_EOS = sum(rho (Y_air e_air + Y_vap e_vap + |U|^2/2) V)  [percent]",
             "dEe", "emaxe", "tmaxe", "dEgoode", "dEie", "dEke")):
        out += [tab, h, "-"*len(h)]
        for r in rows:
            out.append("%-22s %-21s" % (r["name"], r["status"])
                       + "".join("%9s" % pc(r[kE][t]) for t in CHECK_TIMES)
                       + "%10s %6.1f %8s (%4.1f) %9s %9s" % (pc(r[kmax]), r[kt]*1e9, pc(r[kg]), r["tgood"]*1e9,
                                                             pc(r[ki]), pc(r[kk])))
        out.append("")
    out.append("max|dE|: over all written times (for a STOPped run this is the health-check dump of the broken state).")
    out.append("last good: drift at the last written time before the stop (= the last time for finished runs).")
    out.append("dE_int, dE_kin: change of E_int and E_kin from t = 0 to the last good time, as % of E_tot(0).")
    out.append("  (In a shock tube E_int is physically converted into E_kin, so these two are large and opposite")
    out.append("   even for an exact solution; their SUM is the drift. Table 1b separates the error for case A.)")
    out.append("")
    out.append("Table 1b (case A / Af): error against the exact Riemann solution at the last good time, % of E_tot(0)")
    hb = "%-22s %8s %16s %16s %16s %16s" % ("run", "t [ns]", "E_kin exact", "E_kin solver", "err E_int", "err E_kin")
    out += [hb, "-"*len(hb)]
    for r in rows:
        if r["exact"]:
            ei, ek, ekx = r["exact"]
            out.append("%-22s %8.1f %16s %16s %16s %16s" % (r["name"], r["tgood"]*1e9, pc(ekx), pc(ekx + ek), pc(ei), pc(ek)))
    out.append("err E_int = E_int - (E_tot(0) - E_kin_exact); err E_kin = E_kin - E_kin_exact; they add up to the drift.")
    out.append("")
    out.append("Table 3: mass, and dissipative-mode dropped KE cooling")
    h3 = "%-22s %11s %11s %11s %11s %13s   %s" % ("run", "dM_air(0)", "dM_vap(0)", "dM_air(1ns)", "dM_vap(1ns)",
                                               "rhoEOS-sumc", "dissipative: largest dropped [K]; T solves > 100 K / total")
    out += [h3, "-"*len(h3)]
    for r in rows:
        d = r["drops"]
        ds = ("%.4g; %d / %d" % (max(d), sum(v > 100 for v in d), len(d))) if d else "-"
        out.append("%-22s %11.2e %11.2e %11.2e %11.2e %13.2e   %s" % (
            r["name"], r["mA"], r["mV"], r["mA1"], r["mV1"], r["mEos"], ds))
    out.append("dM_i(0): largest signed (M_i(t) - M_i(0))/M_i(0) over all written times; M_i(0) is the exact EOS rebuild,")
    out.append("         so this includes a ~1e-8 offset from writePrecision 8 in the written c_i files.")
    out.append("dM_i(1ns): the same relative to the first written time (1 ns), i.e. file-to-file; ~1e-8 = output round-off.")
    out.append("max(rhoEOS-sumc): largest (sum(rho V) - sum(c V))/sum(c V), the gap between the two density definitions.")
    out.append("")
    out += lines_chk
    txt = "\n".join(out)
    print(txt)
    open("energy_table.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    runs = [r for r in sys.argv[1:] if os.path.isdir(r)]
    if not runs:
        sys.exit(__doc__)
    main(runs)
