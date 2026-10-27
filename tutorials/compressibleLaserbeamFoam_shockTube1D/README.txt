Shock tube for compressibleLaserbeamFoam
=========================================

1D tube, 600 um long (x = -200..400 um), 3.125 um cells, diaphragm at x = 0.
Left: pure Fe vapour (metal1vapour). Right: the "air" phase at 1 bar, 300 K.
No liquid, no laser (power 0), no AMR. Same thermophysical files, schemes and
solver settings as the 50mj_CSD3 case (fvSolution PIMPLE incl. residualControl),
plus DIAGNOSTICS with stopOnViolation.

Case A  100 bar, 8000 K   exact Riemann solution valid (all < 10 kK)
Case B  1000 bar, 12000 K  the 3D plume state before the 5.133 ns failure

Run on CSD3 (needs the energy-ke-option build):
    source ~/rds/rds-ukaea-ap002-mOlK9qn0PlQ/ir-roch3/OpenFoamv2506/of2506.sh
    cd shockTube && ./Allrun
Then paste summary.txt back.

What to look for
  A: does each mode reproduce the exact shock speed, p*, u*, T air, T vapour?
     Expected at 60 ns: shock 86.6 um, contact 61.7 um, p* 2.48e6 Pa,
     u* 1029 m/s, T air 2138 K, T vapour 4579 K.
  B: does full blow up / hit the health check, and do contErr or none not?
     log.solver also prints "max |K*contErr| as dT per step" every T solve.

Variants (third field):  ./Allrun A:contErr:o5 A:contErr:o5nr A:contErr:o10
  o5 / o10  nOuterCorrectors 5 / 10 with residualControl off
  nr        all under-relaxation factors 1
  o5nr      both
Local (WSL) fixes from the first run are now in the template: AMR libs removed
from controlDict, sampleDict wrapped in functions{} with readFields and midPoint.
v3: variant 'mc' = mixtureCv mass (needs the mixture-cv-option build, c2205c53).
    Variants combine with '+':  ./Allrun B:contErr:mc B:contErr:o5nr+mc
v4: variant 'eb' = DIAGNOSTICS energyBudget (per-term T budget of the cell with the
    largest |dT|, printed when |dT| >= 200 K). Needs mixture-cv-option @ f55a0c0e.
    grep -A16 "\[budget\] t =" runs/<name>/log.solver | tail -60
v5: mode 'dissipative' (kineticEnergyTerms dissipative) and variant 'ebi' (budget every
    outer corrector). Needs mixture-cv-option @ 21b71f16.
