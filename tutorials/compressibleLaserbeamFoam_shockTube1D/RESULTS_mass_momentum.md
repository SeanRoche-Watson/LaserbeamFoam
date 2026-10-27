# momentumDensity eos | mass (branch cc-mass-momentum), 27 Oct 2026

## Status: compiled and tested (28 Oct 2026); results at the end of this file
It compiled cleanly on CSD3 and in WSL (ESI v2506), with no compile fixes; the only warning is an unused
`rhoMom` in createFields.H. The design notes below were written before the first build.

## The option
`system/fvSolution`, PIMPLE: `momentumDensity eos;` (default, original) or `momentumDensity mass;`.
- `rhoM = max(sum_i c_i, 1e-3)`, a volScalarField built from the registered `c.<phase>` fields (rhoM.H).
  It is always built (for the diagnostic line) and written only in mass mode. It is never read.
- `rhoMom` is a reference to `rho` (eos) or `rhoM` (mass). The equations use rhoMom, so eos runs the
  original expressions on the same object.
- Per time step: rhoM is recomputed, then `rhoM.oldTime() = rhoM`, right after `rho.oldTime() = rho` and
  before `mixture.solve()`. rhoM is recomputed again straight after `mixture.solve()`.
- Log line at the end of every step (both modes):
  `momentumDensity <mode>: max |rho/rhoM - 1| = ..., max |ddt(rhoM) + div(rhoPhi)| = ... kg/m3/s`.
- Budget block (mass mode): `rho` is rhoM, plus an extra line `EOS rho ..., rho/rhoM = ...`. The momentum check
  uses rhoMom.

## Is rhoPhi_ the flux that updates c_i? Yes (checked in multiphaseMixtureThermo.C)
- `rhoPhi_` is zeroed (l. 2446). Only `rhoPhi_ += Fi` (l. 2899) adds to it, and that is the same `Fi` passed to
  `MULES::explicitSolve(geometricOneField(), c_i, Fi, zeroField(), SuMass_i)`. Fi includes the limited group
  share and the Fickian part. The 5-argument explicitSolve does no limiting:
  c_i = c_i.oldTime() + dt*(SuMass_i - div(Fi)), with boundary faces included in both.
- c_i.oldTime() == c0_ == c_i at the first solveAlphas call of the step (l. 1757-1768). Nothing else changes
  c_i internal values between steps, so sum_i c_i at the start of the step = sum_i c_i.oldTime().
- SuMass is paired (l. 2188-2189: liquid -r, vapour +r), so sum_i SuMass_i = 0.
- nAlphaSubCycles > 1 is a fatal error (l. 1245), so there is no sub-cycle accumulation.
- Departures from exact: the purge (|c| < cPurgeTol = 1e-10 set to 0, after the update), the 1e-3 floor, and
  round-off. So `max |ddt(rhoM)+div(rhoPhi)|` should be ~0 in mass mode, apart from floored or purged cells.

## Where rho is used: changed / not changed (mass mode)
| place | term | changed? |
|---|---|---|
| UEqn | `fvm::ddt(rho,U)` | rhoM |
| UEqn | `ghf*snGrad(rho)` (momentum predictor) | rhoM |
| UEqn | `divDevRhoReff(U)` | no: laminar mu = rho*nu, rho cancels; turbulence model keeps EOS rho |
| pEqn | `interpolate(rho*rAU)*ddtCorr(U,phi)` | rhoM (current level, see below) |
| pEqn | `ghf*snGrad(rho)` in phig | rhoM |
| pEqn | rAU = 1/UEqn.A() | follows implicitly (A now holds rhoM/dt), so laplacian(rAUf,p_rgh) is the rhoM one |
| pEqn | `p = p_rgh + rho*gh`, `p_rgh = p - rho*gh` | rhoM, so -grad p_rgh - gh grad rhoM = -grad p + rhoM g (consistent body force). With g ~ 10 m/s2 on a 1e-4 m domain this is < 1 Pa. |
| pEqn | p_rghEqnComps (phase rho, psi), dgdt, correctRho, `rho = mixture.rho()` | no (compressibility / EOS) |
| TEqn | `fvm::ddt(rho,T)`, `Sp(ddt(rho)+div(rhoPhi),T)` | rhoM |
| TEqn | contErrT, `ddt(rho,K)` in the work term | rhoM (so contErr ~ 0 and full ~ contErr in mass mode) |
| TEqn | K*contErr / dissipative diagnostics `/rho` | rhoM |
| TEqn | budget rho0, keB, momentum-check rho | rhoM |
| TEqn | TRHS `ddt(rho,epsilon1)+div(rhoPhi,epsilon1)` (latent heat of melting) | **no**: not in the spec; only acts where meltMask > 0 (condensed). For full consistency it arguably should be rhoM too. Not changed without a test. |
| TEqn | mass_dot | no: rCv*L*r, no density |
| solverHealth, setDeltaT | rho | no (EOS rho is still the thermodynamic density) |

## Points I am not sure about
- **ddtCorr time level.** With Euler, `fvm::ddt(rhoM,U)` puts rhoM.oldTime()*U.oldTime()/dt into H. So
  strictly, the factor that replaces flux(U.oldTime()) by phi.oldTime() is rhoM.oldTime()*rAU. The original code
  uses the current rho (interFoam-family convention), and I kept that with the current rhoM. The two differ by
  one step's mass change. That is small unless contErr is large, and contErr is ~0 in mass mode by construction.
  Switching to `rhoM.oldTime()` would be a one-line test.
- **Mixed density in the pressure equation.** The momentum inertia (rAU) now uses rhoM. The compressibility
  part still uses the phase EOS densities. This is intended, but it makes the pressure equation a hybrid.
  Where rho and rhoM differ 10-45x (3D), the pressure Laplacian changes by that factor. That is the point of
  the test, but it may change stability either way.
- **ddt schemes.** Only Euler was considered. backward would need rhoM.oldTime().oldTime(). The
  storeOldTimes mechanism provides it, but its value is the floored sum c of two steps back. That is fine,
  but untested.
- **p_rgh boundary conditions that look up `rho`** (prghPressure, prghTotalPressure) still use the EOS rho.
  The shock tube and the planned 2D case use waveTransmissive / zeroGradient, so this does not matter there.

## Self-review: time levels, restart, parallel, AMR
- **Time levels.** rhoM.oldTime() is assigned explicitly once per step from the c fields. It is not
  copied from the previous step's rhoM, so it is exact regardless of mapping. During the outer correctors
  rhoM is constant (c_i only change at the first corrector).
- **Restart.** No file is needed. On the first step after restart, rhoM.oldTime() comes from the c fields read
  from the start time. Those are exactly what solveAlphas uses as c0_, to the c fields' write precision.
- **Parallel.** The sum includes processor-patch values, which are the c fields' halo values after their
  correctBoundaryConditions. The diagnostics use gMax.
- **AMR.** rhoM is a registered field, so dynamicRefineFvMesh maps it and its old time. But it is recomputed
  from the mapped c fields after mesh.update() and before use, so the mapping of rhoM itself does not matter.
  Limitation: with moveMeshOuterCorrectors, a mesh change on a later corrector is not followed by a rhoM
  recompute. rhoPhi is stale then as well, since mixture.solve() only runs on the first corrector.
- **eos bit-identity.** The eos path evaluates the same expressions on the same field (rhoMom is rho). In pEqn
  the eos branch keeps `mixture.rho()*gh`. The added work (rhoM, a registered temp `sumC`, the ddt(rhoM)
  diagnostic) has no side effects on solver fields. Bit-identity still has to be confirmed by a run
  (compare log.solver numbers / time directories against mixture-cv-option f73ba8b2).

## Tests to run (section 3 of the task, nothing run yet)
Build: `cd src && ./Allwmake`, then `cd applications/solvers/compressibleLaserbeamFoam && ./Allwmake`.
Without libROUNDSchemes, add the `vl` variant. CONTEXT.md: div(phi,alpha) only moves beta, so gas results don't change.
```
./Allrun A:contErr:o10+nc+vl B:contErr:o10+nc+vl                     # a) eos, vs CONTEXT.md numbers
./Allrun A:contErr:o10+nc+vl+mm Af:contErr:o10+nc+vl+mm B:contErr:o10+nc+vl+mm   # b) mass
python3 energy_check.py runs/*; python3 summarise.py runs/*
./Allrun B:full:o10+vl B:full:o10+vl+mm                                      # c) stress, gasPackingCap true (default)
```
d) The 2D plume-edge case (tutorials/compressibleLaserbeamFoam_plumeEdge2D) is not created yet. It was part of
the run phase.


## Results (WSL, 28 Oct 2026; tables and plots in shockTube_work/mm_results/)

### Bit-identity of eos
A:contErr:o10+nc and B:contErr:o10+nc (and +ebi) on this branch with the default `momentumDensity eos` are
**bit-identical** to the packing2 runs (43bae9b1, writeInterval 1e-9). Every written field file at every write
matches (4779 files per run). The logs differ only in the VTK output paths.

### mass vs eos in 1D (contErr, o10, gasPackingCap false)
| run (60 ns) | head | T vap | p vap | shock | p* | u* | T air | rho air | energy drift |
|---|---|---|---|---|---|---|---|---|---|
| exact | -84.5 | 4579 | 2.48e6 | 86.6 | 2.48e6 | 1029 | 2138 | 5.57 | 0 |
| A eos | -95.3 | 4794 | 2.709e6 | 95.3 | 2.518e6 | 1030 | 2567 | 4.67 | +0.49% |
| A mass | -95.3 | 4793 | 2.703e6 | 95.3 | 2.517e6 | 1032 | 2518 | 4.77 | +0.39% |
| Af eos | -87.9 | 4649 | 2.526e6 | 89.5 | 2.466e6 | 1022 | 2265 | 5.19 | +0.16% |
| Af mass | -87.9 | 4642 | 2.509e6 | 89.5 | 2.501e6 | 1030 | 2259 | 5.28 | +0.22% |
| B eos | | | | | | | | | -0.57% |
| B mass | | | | | | | | | -0.22% |

- **Size of the change:** mass changes the A/Af plateaus by <= 2%. Over the tube, L1 eos-vs-mass differences
  are 0.1-0.7% in p, U, T and rho. The larger pointwise differences (2-10%) are confined to the shock-front
  cell, a one-cell shift.
- **No sign of a bug.** In Af the drift grows (0.16 -> 0.22%). The E_int error is unchanged (+0.41%) and
  E_kin moves closer to exact (-0.24 -> -0.19%), so the extra drift is kinetic energy that eos was missing.
- **Mass and c_i:** mass is conserved to round-off, and no c_i is ever negative.
- **The new log line:**
  - `max |rho/rhoM - 1|`: per step it is typically 2.5-3% and at most 7-9%, always at the shock-front cell.
    It is the same in eos and mass; it does not grow.
  - `max |ddt(rhoM)+div(rhoPhi)|`: median 3e-3 to 6e-2 kg/m3/s, max 8-26 kg/m3/s. Times dt this is at most
    1e-10 kg/m3 = cPurgeTol, so it is the purge plus round-off; it peaks right after each write, when dt is
    smallest. The line is mode-independent by construction (rhoM is always sum c), which confirms that
    rhoPhi is the flux that moves the c_i.
- **B, contErr o10:** mass removes the cold dip (eos 4600 K vs mass 5800 K at the contact) and the velocity
  spike at the contact. The contact plateau is ~2100 m/s (Riemann u* ~2000).

### Stress test: B, kineticEnergyTerms full, gasPackingCap true
| setting | eos | mass |
|---|---|---|
| resCtrl (default) | **stops 26.70 ns**; mismatch rises 0.1 -> 0.45 from 10 to 25 ns; contact runs away to ~5000 m/s | 60 ns, mismatch <= 0.14 |
| o5nr | **stops 10.15 ns**; mismatch spikes to 1.4 at 9.5 ns | 60 ns, mismatch <= 0.18 |
| o10 | 60 ns (mismatch <= 0.13) | 60 ns (<= 0.17) |

So where the EOS density and the transported mass drift apart, mass stops the spurious acceleration, as it
did in the 3D 10 mJ restart. The cap-on shock is still wrong in both modes (a T/U spike at the front), as before.

### 2D plume edge (tutorials/compressibleLaserbeamFoam_plumeEdge2D), contErr o10 nc eb, 5 ranks, 0-20 ns
- **Both runs reach 20 ns** (eos and mm, 66 min each). No health-check stop, and no energy-budget block
  (no cell |dT| >= 200 K per step).
- **The two modes are practically the same:**
  - Lateral vapour front along the wall: eos 42 / 61 / 86 um, mass 42 / 61 / 83 um at 2 / 10 / 20 ns.
  - Shock 89 um at 20 ns in both.
  - max |U| 1930-2020 m/s in both.
- **Mismatch:** max |rho/rhoM - 1| per step is ~0.1 at start-up, then 0.04-0.05.
- **Geometry:** the lateral extent equals the axial height (a near-hemicylinder). There is no wall jet and
  no fast dense-vapour ring like the 3D failure.
- **This case does not reproduce the 3D mismatch.** Its initial vapour is at rest and in pressure
  equilibrium within the disc; the 3D plume at 5 ns is laser-driven, with large rho/rhoM (10-45x).
- **Next 2D test to make it bite:** start from a mapped 3D plume state, or give the disc an outward velocity
  and a non-equilibrium (rho != sum c) start.

### Tools
- `tutorials/tools/plumeFront.py`: plume and shock height (axis / max / hemisphere-equivalent), lateral
  extent at the surface, max |U| / T / p, vapour mass and the rho/rhoM-mismatch volume, per time.
  - Input: ascii or binary, gz, AMR (polyMesh per time).
  - Output: a CSV and a log-log plot with an optional Sedov-Taylor line.
  - Tested on 2D runs (ascii = binary) and a synthetic 3D hemisphere (R_eq 30.4 um for r = 30 um).
- `LITERATURE_plume_front.md`: published fronts and the blast-wave estimate for 10 mJ in 1 bar argon. The
  shock should be at ~250-280 / 330-370 / 450-540 um at 10 / 20 / 50 ns, with the vapour front somewhat behind.
