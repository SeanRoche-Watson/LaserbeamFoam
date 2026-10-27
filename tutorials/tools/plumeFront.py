#!/usr/bin/env python3
"""Plume and shock fronts of a compressibleLaserbeamFoam LIBS case, one row per time.

    python3 plumeFront.py <case> [options]

Works on a reconstructed case (time directories in <case>/), ascii or binary, with or
without AMR (a polyMesh per time directory). Cell centres and volumes come from the
files C and V in each time directory. If they are missing, the script runs
    postProcess -funcs '(writeCellCentres writeCellVolumes)' -time <list>
(OpenFOAM must be sourced) and deletes the files it wrote at the end unless --keep-geometry.

Geometry convention (as in the 3D cases): the gas is at y < 0, the metal surface at
y = ySurf (default 0), the laser axis is the y line through (x0, z0) (default 0, 0).
A one-cell-thick mesh in z is treated as planar 2D (r = |x - x0|).

Per time it reports
  plume_h_axis / plume_h_max  max height above the surface (ySurf - y) of cells with
                              sum alpha.metal{1..4}vapour > --vap, on the axis
                              (cells within one cell size of it) / anywhere
  plume_R_eq                  hemisphere-equivalent radius (3V/2pi)^(1/3) of those cells
                              (2D: half-disc sqrt(2A/pi))
  plume_r_surf                lateral extent: max r of those cells in the first gas cell row
  shock_*                     the same for p > --pshock (default 2e5 Pa)
  max_U, max_T, max_p         over all cells
  m_vap                       total vapour mass sum c.metal{1..4}vapour * V (kg; 2D: for
                              the mesh thickness)
  V_mismatch                  volume where |rho/rhoM - 1| > --mismatch (0.1), with rhoM the
                              written rhoM field if present, else sum of all c.* fields
Output: <case>/plumeFront.csv and plumeFront.png (log-log plume and shock height vs time,
with a Sedov-Taylor line if --E (3D, J) or --El (2D, J/m) is given).
"""
import argparse, glob, gzip, os, re, shutil, struct, subprocess, sys
import numpy as np

VAP = ["alpha.metal%dvapour" % i for i in range(1, 5)]
CVAP = ["c.metal%dvapour" % i for i in range(1, 5)]


# ------------------------------------------------------------------ field reader
def _open(path):
    if os.path.exists(path):
        return open(path, "rb").read()
    if os.path.exists(path + ".gz"):
        return gzip.open(path + ".gz", "rb").read()
    return None


def _header(buf):
    m = re.search(rb"FoamFile\s*\{(.*?)\}", buf, re.S)
    h = m.group(1) if m else b""
    fmt = re.search(rb"format\s+(\w+)\s*;", h)
    binary = bool(fmt and fmt.group(1) == b"binary")
    arch = re.search(rb'arch\s+"([^"]*)"', h)
    scalar = 8
    if arch:
        s = re.search(rb"scalar=(\d+)", arch.group(1))
        if s:
            scalar = int(s.group(1))//8
    return binary, scalar, (m.end() if m else 0)


def _list(buf, pos, ncomp, binary, scalar):
    """Parse an OpenFOAM list starting at buf[pos:] ('N(' ... ')' or 'N{v}')."""
    m = re.compile(rb"\s*(\d+)\s*([({])").match(buf, pos)
    if not m:
        raise ValueError("no list at %d" % pos)
    n = int(m.group(1))
    if m.group(2) == b"{":
        e = buf.index(b"}", m.end())
        v = np.array(buf[m.end():e].replace(b"(", b" ").replace(b")", b" ").split(), float)
        return np.tile(v, (n, 1)).squeeze() if ncomp > 1 else np.full(n, v[0])
    if binary:
        a = np.frombuffer(buf, dtype="<f%d" % scalar, count=n*ncomp, offset=m.end()).astype(float)
    else:
        e = buf.index(b"\n)", m.end()) if ncomp > 1 else buf.index(b")", m.end())
        body = buf[m.end():e]
        if ncomp > 1:
            body = body.replace(b"(", b" ").replace(b")", b" ")
        a = np.array(body.split(), float)
    return a.reshape(n, ncomp) if ncomp > 1 else a


def read_field(path, n=None):
    """internalField of a volScalarField / volVectorField file (None if absent)."""
    buf = _open(path)
    if buf is None:
        return None
    binary, scalar, h = _header(buf)
    i = buf.index(b"internalField", h) + len(b"internalField")
    ncomp = 3 if re.search(rb"class\s+vol(Vector)Field", buf[:h]) else 1
    m = re.compile(rb"\s*uniform\s+").match(buf, i)
    if m:
        e = buf.index(b";", m.end())
        v = np.array(buf[m.end():e].replace(b"(", b" ").replace(b")", b" ").split(), float)
        if n is None:
            raise ValueError("uniform field %s needs the cell count" % path)
        return np.tile(v, (n, 1)) if ncomp > 1 else np.full(n, v[0])
    m = re.compile(rb"\s*nonuniform\s+List<\w+>").match(buf, i)
    return _list(buf, m.end(), ncomp, binary, scalar)


def read_points_z(path):
    buf = _open(path)
    binary, scalar, h = _header(buf)
    # the list starts on its own line after the banner comment: "N\n("
    m = re.compile(rb"\n(\d+)\s*\n\(").search(buf, h)
    p = _list(buf, m.start() + 1, 3, binary, scalar)
    return p[:, 2].min(), p[:, 2].max()


# ------------------------------------------------------------------ helpers
def times_of(case):
    out = []
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.isdir(os.path.join(case, d)):
            out.append((t, d))
    return sorted(out)


def ensure_geometry(case, tdirs, keep):
    missing = [d for _, d in tdirs
               if _open(os.path.join(case, d, "C")) is None or _open(os.path.join(case, d, "V")) is None]
    if not missing:
        return []
    if shutil.which("postProcess") is None:
        sys.exit("C/V missing in %d time directories and postProcess not found: source OpenFOAM first" % len(missing))
    print("writing cell centres/volumes for %d times with postProcess ..." % len(missing), flush=True)
    for k in range(0, len(missing), 50):
        chunk = missing[k:k + 50]
        cmd = ["postProcess", "-case", case, "-funcs", "(writeCellCentres writeCellVolumes)",
               "-time", ",".join(chunk)]
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if r.returncode != 0:
            sys.exit("postProcess failed:\n" + r.stdout[-3000:])
    if keep:
        return []
    return [os.path.join(case, d, f) for d in missing for f in ("C", "Cx", "Cy", "Cz", "V")]


def mesh_thickness(case, tdir):
    for p in (os.path.join(case, tdir, "polyMesh", "points"), os.path.join(case, "constant", "polyMesh", "points")):
        if _open(p) is not None:
            z0, z1 = read_points_z(p)
            return z1 - z0
    raise SystemExit("no polyMesh/points found")


def front(mask, h, r, axis, first_row, Vc, dim, thick):
    if not mask.any():
        return np.nan, np.nan, 0.0, np.nan
    h_max = h[mask].max()
    h_axis = h[mask & axis].max() if (mask & axis).any() else np.nan
    vol = Vc[mask].sum()
    if dim == 3:
        Req = (3*vol/(2*np.pi))**(1/3)
    else:
        Req = np.sqrt(2*(vol/thick)/np.pi)
    r_s = r[mask & first_row].max() if (mask & first_row).any() else np.nan
    return h_axis, h_max, Req, r_s


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--axis", nargs=2, type=float, default=[0.0, 0.0], metavar=("X0", "Z0"))
    ap.add_argument("--surface", type=float, default=0.0, help="surface y (gas at y < ySurf)")
    ap.add_argument("--vap", type=float, default=0.01, help="vapour volume-fraction threshold")
    ap.add_argument("--pshock", type=float, default=2e5, help="shock threshold [Pa]")
    ap.add_argument("--mismatch", type=float, default=0.1, help="|rho/rhoM - 1| threshold")
    ap.add_argument("--tmin", type=float, default=0.0)
    ap.add_argument("--tmax", type=float, default=np.inf)
    ap.add_argument("--E", type=float, help="3D: pulse energy [J] for a Sedov-Taylor line")
    ap.add_argument("--El", type=float, help="2D: energy per unit length [J/m] for a cylindrical blast line")
    ap.add_argument("--rho0", type=float, default=1.6, help="ambient density [kg/m3] (argon 1 bar 300 K: 1.6)")
    ap.add_argument("--out", default=None, help="CSV path (default <case>/plumeFront.csv)")
    ap.add_argument("--keep-geometry", action="store_true", help="keep the C/V files written by postProcess")
    a = ap.parse_args()

    case = a.case.rstrip("/")
    tdirs = [(t, d) for t, d in times_of(case) if a.tmin <= t <= a.tmax]
    if not tdirs:
        sys.exit("no time directories in %s (reconstruct the case first)" % case)
    written = ensure_geometry(case, tdirs, a.keep_geometry)
    x0, z0 = a.axis
    cols = ["time", "dim", "ncells",
            "plume_h_axis", "plume_h_max", "plume_R_eq", "plume_r_surf",
            "shock_h_axis", "shock_h_max", "shock_R_eq", "shock_r_surf",
            "max_U", "max_T", "max_p", "m_vap", "V_mismatch", "rhoM_source"]
    rows = []
    try:
        for t, d in tdirs:
            td = os.path.join(case, d)
            C = read_field(os.path.join(td, "C"))
            n = len(C)
            Vc = read_field(os.path.join(td, "V"), n)
            g = lambda f: read_field(os.path.join(td, f), n)
            dim = 2 if np.ptp(C[:, 2]) < 1e-12 else 3
            thick = mesh_thickness(case, d) if dim == 2 else 1.0
            dx = (Vc/thick)**0.5 if dim == 2 else Vc**(1/3)
            h = a.surface - C[:, 1]
            r = np.abs(C[:, 0] - x0) if dim == 2 else np.hypot(C[:, 0] - x0, C[:, 2] - z0)
            axis = r < dx
            first_row = (h > 0) & (h < dx)
            vap = np.zeros(n)
            for f in VAP:
                v = g(f)
                if v is not None:
                    vap += v
            p, T, U = g("p"), g("T"), g("U")
            P = front(vap > a.vap, h, r, axis, first_row, Vc, dim, thick)
            S = front((p > a.pshock) & (h > 0), h, r, axis, first_row, Vc, dim, thick)
            cv = [g(f) for f in CVAP]
            m_vap = sum((c*Vc).sum() for c in cv if c is not None) if any(c is not None for c in cv) else np.nan
            rho, rhoM, src = g("rho"), g("rhoM"), "rhoM"
            if rhoM is None:
                cs = [read_field(f, n) for f in glob.glob(os.path.join(td, "c.*"))
                      if re.fullmatch(r"c\.[A-Za-z0-9]+", os.path.basename(f))]
                rhoM, src = (np.maximum(sum(cs), 1e-3), "sum_c") if cs else (None, "none")
            Vm = Vc[np.abs(rho/rhoM - 1) > a.mismatch].sum() if (rho is not None and rhoM is not None) else np.nan
            rows.append([t, dim, n, *P, *S, np.sqrt((U**2).sum(1)).max(), T.max(), p.max(), m_vap, Vm, src])
            print("t = %-10s plume h %.3g m (axis %.3g), shock h %.3g m, lateral %.3g m, max|U| %.4g, max T %.5g"
                  % (d, P[1], P[0], S[1], P[3], rows[-1][11], rows[-1][12]), flush=True)
    finally:
        for f in written:
            if os.path.exists(f):
                os.remove(f)

    out = a.out or os.path.join(case, "plumeFront.csv")
    with open(out, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for row in rows:
            fh.write(",".join(x if isinstance(x, str) else ("%d" % x if isinstance(x, (int, np.integer)) else "%.6g" % x)
                              for x in row) + "\n")
    print("wrote", out)
    plot(rows, a, os.path.splitext(out)[0] + ".png", case)


def plot(rows, a, png, case):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: no plot")
        return
    R = np.array([[r[0], r[3], r[4], r[5], r[7], r[8], r[9]] for r in rows], float)
    R = R[R[:, 0] > 0]
    if len(R) == 0:
        return
    t = R[:, 0]*1e9
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(t, R[:, 2]*1e6, "o-", color="#d62728", ms=3, label="plume height (max)")
    ax.loglog(t, R[:, 1]*1e6, "s--", color="#d62728", ms=3, mfc="none", label="plume height (axis)")
    ax.loglog(t, R[:, 3]*1e6, ":", color="#d62728", label="plume R_eq")
    ax.loglog(t, R[:, 5]*1e6, "o-", color="#1f77b4", ms=3, label="shock height (max, p > %.0e)" % a.pshock)
    ax.loglog(t, R[:, 4]*1e6, "s--", color="#1f77b4", ms=3, mfc="none", label="shock height (axis)")
    ax.loglog(t, R[:, 6]*1e6, ":", color="#1f77b4", label="shock R_eq")
    tt = np.logspace(np.log10(t.min()), np.log10(t.max()), 50)*1e-9
    if a.E:
        # blast on a surface: the energy fills a hemisphere, i.e. a sphere of energy 2E
        ax.loglog(tt*1e9, (2*a.E/a.rho0)**0.2*tt**0.4*1e6, "k-", lw=1,
                  label="Sedov-Taylor (2E/rho0)^1/5 t^2/5, E = %g J" % a.E)
    if a.El:
        ax.loglog(tt*1e9, (2*a.El/a.rho0)**0.25*tt**0.5*1e6, "k-", lw=1,
                  label="cylindrical blast (2E_l/rho0)^1/4 t^1/2, E_l = %g J/m" % a.El)
    ax.set_xlabel("t [ns]"); ax.set_ylabel("height above surface [um]")
    ax.set_title(os.path.basename(os.path.abspath(case)), fontsize=10)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(png, dpi=120)
    print("wrote", png)


if __name__ == "__main__":
    main()
