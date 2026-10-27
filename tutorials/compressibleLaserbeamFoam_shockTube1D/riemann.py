#!/usr/bin/env python3
"""Exact Riemann solution for the shock tube (two ideal gases, gamma = 5/3).

Valid for case A (everything stays below 10 kK, where the janaf files give
Cp/R = 2.5 exactly). For case B the left state is above 10 kK, where Cv rises,
so the numbers are only indicative.

  python3 riemann.py A            -> wave speeds and star state
  python3 riemann.py A 6e-8 > f   -> exact profile x rho u p T at t
"""
import sys, math

Ru = 8.314462618
CASES = {
    "A": dict(pL=1e7, TL=8000.0),
    "B": dict(pL=1e8, TL=12000.0),
}
ML, MR = 0.05585, 0.039948       # Fe vapour (metal1vapour), "air" phase (argon)
pR, TR = 1e5, 300.0
G = 5.0/3.0

def state(p, T, M):
    R = Ru/M
    rho = p/(R*T)
    return dict(p=p, T=T, rho=rho, R=R, c=math.sqrt(G*p/rho))

def f(p, s):
    if p > s["p"]:
        A = 2.0/((G + 1)*s["rho"]); Bc = (G - 1)/(G + 1)*s["p"]
        return (p - s["p"])*math.sqrt(A/(p + Bc))
    return 2*s["c"]/(G - 1)*((p/s["p"])**((G - 1)/(2*G)) - 1)

def star(L, R):
    lo, hi = 1e-6, 1e13
    for _ in range(300):
        mid = math.sqrt(lo*hi)
        if f(mid, L) + f(mid, R) > 0: hi = mid
        else: lo = mid
    ps = math.sqrt(lo*hi)
    return ps, 0.5*(f(ps, R) - f(ps, L))

def solution(case):
    L = state(CASES[case]["pL"], CASES[case]["TL"], ML)
    R = state(pR, TR, MR)
    ps, us = star(L, R)
    g1 = (G - 1)/(G + 1)
    rhoRs = R["rho"]*(ps/R["p"] + g1)/(g1*ps/R["p"] + 1)
    S = R["c"]*math.sqrt((G + 1)/(2*G)*ps/R["p"] + (G - 1)/(2*G))
    rhoLs = L["rho"]*(ps/L["p"])**(1/G)
    cLs = math.sqrt(G*ps/rhoLs)
    return dict(L=L, R=R, ps=ps, us=us, S=S, rhoRs=rhoRs, rhoLs=rhoLs,
                TRs=ps/(rhoRs*R["R"]), TLs=ps/(rhoLs*L["R"]),
                headL=-L["c"], tailL=us - cLs)

def sample(sol, x, t):
    L, R = sol["L"], sol["R"]; xi = x/t
    if xi < sol["headL"]:
        return L["rho"], 0.0, L["p"], L["T"]
    if xi < sol["tailL"]:
        u = 2/(G + 1)*(L["c"] + xi)
        c = L["c"] - (G - 1)/2*u
        rho = L["rho"]*(c/L["c"])**(2/(G - 1))
        p = L["p"]*(rho/L["rho"])**G
        return rho, u, p, p/(rho*L["R"])
    if xi < sol["us"]:
        return sol["rhoLs"], sol["us"], sol["ps"], sol["TLs"]
    if xi < sol["S"]:
        return sol["rhoRs"], sol["us"], sol["ps"], sol["TRs"]
    return R["rho"], 0.0, R["p"], R["T"]

if __name__ == "__main__":
    case = sys.argv[1] if len(sys.argv) > 1 else "A"
    s = solution(case)
    if len(sys.argv) > 2:
        t = float(sys.argv[2])
        n = 2000
        for i in range(n + 1):
            x = -2e-4 + 6e-4*i/n
            print("%.6e %.6e %.6e %.6e %.6e" % ((x,) + sample(s, x, t)))
    else:
        print("case %s: left rho %.4g kg/m3, c %.0f m/s" % (case, s["L"]["rho"], s["L"]["c"]))
        print("  p* = %.4g Pa, u* = %.0f m/s, shock %.0f m/s" % (s["ps"], s["us"], s["S"]))
        print("  shocked air: rho %.4g, T %.0f K;  expanded vapour: rho %.4g, T %.0f K"
              % (s["rhoRs"], s["TRs"], s["rhoLs"], s["TLs"]))
        print("  rarefaction head %.0f m/s, tail %.0f m/s" % (s["headL"], s["tailL"]))
