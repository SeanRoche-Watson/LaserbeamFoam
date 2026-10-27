# Context for continuing the compressibleLaserbeamFoam shock-tube investigation

Read this first. It summarises everything done so far (27 Oct 2026).

## The problem
3D 10 mJ LIBS runs (CSD3, 12.8M cells, 160 ranks) fail in the vapour plume after the 5 ns pulse,
at 5.1-5.5 ns. Every failure: a mixed vapour/air cell is cooled to the T floor (250 K), its EOS density
goes to the 1e-3 floor, then |U| -> 1e4-1e5 m/s and T -> 1e6 K (health check stops the run).
Restarts from the same time are deterministic.

## Solver branches (github.com/SeanRoche-Watson/laserbeamfoam)
- 755eb84d  crash diagnostics + write sync (production build on CSD3 so far)
- closure-every-outer 7845b15b: thermophysicalProperties closureEveryOuter (off by default). Tested in 3D: made it worse.
- energy-ke-option 00f190b3: PIMPLE kineticEnergyTerms full | contErr | none (+ K*contErr diagnostic line)
- mixture-cv-option (current, 21b71f16): + mixtureCv volume|mass (thermophysicalProperties),
  + DIAGNOSTICS energyBudget / energyBudgetMinDT / energyBudgetPoint / energyBudgetEveryIter,
  + kineticEnergyTerms dissipative.
All new behaviour is behind options; defaults reproduce the original solver.

## Shock-tube results so far (see the folders in this directory)
- Case A (100 bar, 8000 K, exact Riemann solution valid): `full` wrong at shocks and non-convergent (shock 28% fast at 768 cells);
  `none` loses KE (shock 26% slow); `contErr` gets the shock within 4%.
- p-T coupling: with residualControl only ~2 outer iterations run; rarefaction too slow and expanded vapour ~1000 K too hot.
  o5nr (nOuterCorrectors 5, residualControl off, all relaxation factors 1) fixes it (768 cells: head -87.9 vs -84.5 um,
  T vap 4598 vs 4579 K). nr alone with residualControl on is harmful.
- Remaining A error: p* in shocked air ~25% low, shock ~8% fast (shock capturing), u* fine.
- Energy conservation (closed tube): contErr + o5nr conserves to 0.19% (192) / 0.08% (768 cells). dissipative creates energy
  (up to 32%, worse with resolution): stopgap only. Mass conserved to output round-off.
- Case B (1000 bar, 12000 K, the 3D plume state): full stops at 26 ns; contErr runs but forms a cold pocket at the contact;
  contErr o5nr stops at ~5 ns; dissipative o5nr reaches 60 ns.
- Energy budget in the failing B cells: the KINETIC-ENERGY term cools by 1.4e4-1.6e4 K/step (|U| ~5400 m/s, K ~1.5e7 J/kg,
  Cv ~250 J/kg/K); pressure work 1-2 orders smaller. So kinetic energy is gained without matching pressure work:
  a numerical velocity overshoot at the vapour/air contact.
- mixtureCv mass: no effect. In 1D, EOS rho and sum(c_i) agree to 8e-4 (in 3D they differ 10-45x).

## Current hypothesis and next tests
The B (and 3D) failure is a stability problem at the miscible vapour/air contact, not the energy equation.
Suspects: alpha scheme div(phi,alpha) Gauss ROUNDAplus01 (sharpening), hierarchical alpha recovery in solveAlphas,
pressure non-equilibrium at a contact between gases of different molar mass. Also check janaf (flat-Cv test).
Recommended production settings so far: kineticEnergyTerms contErr; o5nr.

## Contact investigation (contact_results/), 27 Oct
- div(phi,alpha) only transports beta (liquid fraction); gas species use a hard-coded upwind group flux. Alpha schemes have no effect on gas.
- Flat Cv (no janaf knee) delays the B failure ~3x (15.9 vs 4.75 ns) but it still fails the same way. The 10 kK knee makes the
  outer loop oscillate (Cv 287 <-> 1152) in the step before failure.
- Velocity overshoot develops first in a ~50/50 mixed cell (2x contact speed), then a pressure spike in the next cell, which expands to the rho floor.
- Shocked air heats 400 -> 9600 K at CONSTANT density 1.6 instead of compressing ~4x.
- Group flux limiter: caps gas inflow at cPackTol (1.02) x gas EOS density at current p,T; fires every step, blocks gas inflow completely.
  Prime suspect: it prevents gas compression at shocks/contacts (also likely cause of the low p* / fast shock in case A).

## New option (mixture-cv-option @ 43bae9b1)
thermophysicalProperties: gasPackingCap true (default) | false. false = no packing cap for the GAS group (inflow cap and flux-density cap);
positivity outflow limit and condensed-group cap stay. Log prints the cell where the gas inflow cap is strongest.
In the shock tube (no condensed phases) `cPackTol 1e3;` removes only the inflow cap (not the flux-density cap); a partial test without rebuild.

## gasPackingCap results (packing_results/, packing2_results/), 27 Oct
- Default (cap on) at 43bae9b1 is bitwise identical to 21b71f16. The option lives in libmultiphaseVapMixtureThermo.so.
- With the cap on, the inflow cap bites every step at the shock-front cell (A: x = shock, factor ~0.37; B: often 0).
- gasPackingCap false: B contErr reaches 60 ns with NO cold pocket (min T vap 3300-4100 K, contact |U| ~2340 m/s).
  B o5 and o10 (relaxation 0.5 kept) also reach 60 ns, no pocket, contact |U| 2065-2150 m/s (Riemann u* ~2000),
  energy drift -0.6%. B o5nr (no relaxation) still fails at 26.85 ns: the KE term cools a pure-vapour cell
  (23.4 um) 8856 -> 250 K, then the KE term in that floored cell grows each outer corrector (+7e3 -> +8.9e4 K)
  and |U| -> 2.2e5. With relaxation the outer loop converges monotonically. So: relaxation, not the number of
  outer correctors, is what B cannot do without. cPackTol 1e3 (inflow cap only) + o5nr: 34 ns; flat Cv + cap off + o5nr: 18 ns.
- Case A/Af with cap off: p* fixed (2.39-2.52e6 vs 2.48e6; cap on 1.87-1.95e6). Best: contErr o10, Af:
  head -87.9, T vap 4649, p vap 2.53e6, shock 89.5, p* 2.47e6, u* 1022 (exact -84.5/4579/2.48e6/86.6/2.48e6/1029);
  energy drift +0.16% (Af) / +0.49% (A) at 60 ns, slowly growing. o5 cap off: drift 1.5-1.6% growing; default resCtrl: 1.3-1.5%.
  Cost: o10 ~2x o5nr wall time (Af 670 s vs 247 s).
- "Hot/light shocked air" is mostly a measurement artefact: summarise.py's air region includes mixed contact cells and
  shock-ramp cells. In pure-air cells (Af o10): T 2181 K (+2%), rho 5.40 (-3%). Entropy excess is largest next to the
  contact and that zone moves with the widening contact (not with fixed parcels): contact-zone heating, not start-up.
  Behind it a small uniform excess (+12 J/kg/K o10, +25 resCtrl, of a 354 J/kg/K jump) that does not grow: weak systematic
  shock heating. Mass conserved to round-off, c_i never negative, in all cap-off runs.
- Candidate production settings now: kineticEnergyTerms contErr, gasPackingCap false, nOuterCorrectors 5-10 with
  residualControl off and relaxation 0.5 (not o5nr). To test in 3D.

## 3D, 27 Oct evening (mixture-cv-option @ 43bae9b1 / 40cd28c7)
- 3D 10 mJ restart from 5 ns, contErr + gasPackingCap false + o10: failed 5.15 ns. Outer loop diverged in one step
  (K*contErr diagnostic 3.6e3 -> 3e10 K over correctors) in a 7 kbar vapour cell compressed 37% in one step.
- Same + maxCo 0.05 (dt 8.3e-12): a RING of cells ~14 um above the surface, ~40 um off axis, dense vapour (~37 kg/m3)
  at 7-8 km/s (above the ~5 km/s maximum expansion speed of Fe vapour from 12 kK) keeps accelerating; KE term
  -3000 to -9300 K/step, pressure work +150 K, contErr ~ -1e11 kg/m3/s; cells reach the 250 K floor by 5.025 ns.
  So in 3D the EOS density vs transported mass mismatch (absent in 1D) drives spurious acceleration.
- 40cd28c7 adds a "momentum check" line to each [budget] block (rho dU/dt vs |grad p|, rho/sum c).
- Next: option momentumDensity eos|mass (use sum c_i in momentum, KE and T terms), test 1D + new 2D plume-edge case.
