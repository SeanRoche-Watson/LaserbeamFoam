# momentumDensity eos | mass (branch cc-mass-momentum), 27 Oct 2026

## Status: implemented and self-reviewed, NOT compiled or run
The cloud environment used for this could not run OpenFOAM. dl.openfoam.com is blocked by its network proxy
(HTTP 403), so the ESI apt repo could not be added. The Ubuntu archive only has OpenFOAM v1912. The v2506
source tarball could be downloaded from SourceForge, but compiling it was not permitted there. **So nothing below
has been compiled.** Expect small compile fixes on the first wmake.

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
