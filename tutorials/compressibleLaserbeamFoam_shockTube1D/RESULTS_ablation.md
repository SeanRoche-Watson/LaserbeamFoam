# Ablation: which fixes are needed for stability, which only for accuracy (28 Oct 2026)

Branch cc-mass-momentum (33e3963f solver build, no source changes). Full tables: shockTube_work/ablation_results/
(ablation_tables.txt, energy_table.txt, profiles_A/Af/B.png, 2d/fields_V*.png, 2d/fronts.png, 2d/fronts.csv).

| | settings | what it adds |
|---|---|---|
| V0 | `full` | original: eos, gasPackingCap on, residualControl (2 outer iterations per step) |
| V1 | `full:mm` | momentumDensity mass |
| V2 | `full:nc+mm` | + gasPackingCap false |
| V3 | `full:o10+nc+mm` | + 10 outer correctors, residualControl off |
| V4 | `contErr:o10+nc+mm` | kineticEnergyTerms contErr (current production settings) |

- 1D: writeInterval 1e-9, serial, run 4 at a time on 12 cores. Reused from 27 Oct: all V4 runs, B V0 (B_full) and B V1 (B_full_mm).
- 2D plume edge: 0-20 ns, writeInterval 1e-9, 2 ranks each, 4 at a time.
- The time step is the same within a case (maxDeltaT-limited), so every difference in cost is outer iterations.
- **Wall times are only good to ~±30%:** runs shared the machine with different neighbours. V3 and V4 do identical
  work and still differ by 12-50%.

## 1D

| case | | status | wall s | max \|rho/rhoM-1\| | max cont. err | energy drift | shock um | p* Pa | u* m/s | T vap K | rho air |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A | exact | | | | | 0 | 86.6 | 2.48e6 | 1029 | 4579 | 5.57 |
| A | V0 | 60 ns | 48 | 0.062 | 14 | +2.70% | 89.1* | 1.81e6* | 977 | 5483 | 4.04 |
| A | V1 | 60 ns | 43 | 0.065 | 14 | +2.16% | 79.7* | 2.80e6 | 898 | 5514 | 4.73 |
| A | V2 | 60 ns | 62 | 0.070 | 14 | +1.61% | 79.7* | 2.95e6 | 864 | 5582 | 4.79 |
| A | V3 | 60 ns | 168 | 0.073 | 26 | +0.39% | 95.3 | 2.52e6 | 1032 | 4793 | 4.77 |
| A | V4 | 60 ns | 251 | 0.073 | 26 | +0.39% | 95.3 | 2.52e6 | 1032 | 4793 | 4.77 |
| Af | V0 | 60 ns | 427 | 0.071 | 8 | +4.86% | 111.3* | 1.15e6* | 1122* | 5313 | 3.43 |
| Af | V1 | 60 ns | 263 | 0.056 | 8 | +2.31% | 82.4 | 2.76e6 | 894 | 5419 | 4.54 |
| Af | V2 | 60 ns | 362 | 0.063 | 8 | +1.50% | 77.7 | 2.92e6 | 853 | 5512 | 4.89 |
| Af | V3 | 60 ns | 1051 | 0.071 | 8 | +0.22% | 89.5 | 2.50e6 | 1030 | 4642 | 5.28 |
| Af | V4 | 60 ns | 930 | 0.071 | 8 | +0.22% | 89.5 | 2.50e6 | 1030 | 4642 | 5.28 |
| B | V0 | **STOP 26.70 ns** | 36 | **0.445** | 2.5 | +11.6% @26 ns | | | contact 4560 (Riemann ~2000) | | |
| B | V1 | 60 ns | 83 | 0.142 | 2.2 | +1.39% | 151.6 | 1.13e7 | 2079 | 6756 | 9.56 |
| B | V2 | 60 ns | 86 | 0.091 | 7.0 | +2.98% | 167.2 | 1.46e7 | 1966 | 7205 | 6.26 |
| B | V3 | 60 ns | 209 | 0.091 | 7.2 | -0.23% | 182.8 | 9.82e6 | 2113 | 6459 | 6.07 |
| B | V4 | 60 ns | 251 | 0.091 | 7.2 | -0.23% | 182.8 | 9.82e6 | 2113 | 6459 | 6.07 |

- Continuity error is max |ddt(rhoM)+div(rhoPhi)| in kg/m3/s.
- Energy drift is E_tot with transported masses, at the last good time.
- **\* Unreliable plateau metrics:** with 2 outer iterations, V0-V2 have a p/U/T/rho spike at the shock, so the
  automatic shock/p*/u* picks are not trustworthy there (see profiles_A.png). The robust signals are the profiles
  and the energy drift:
  - V0-V2 have a shock that lags the exact one, a contact plateau p that is 15-20% high, and u* 10-15% low.
- **The energy error is mostly an E_int excess** (Table 1b of energy_table.txt): +4.1 to +4.8% of E_tot(0) for
  V0-V2, against +0.4 to +0.9% for V3/V4.
- No negative c_i in any run.

## 2D plume edge (P, 160x100, 1000 bar / 12 kK half-disc on the wall)

| | status | wall s | max \|rho/rhoM-1\| | max cont. err | wall front vapour / shock um, 5 ns | 10 ns | 20 ns | max\|U\| 20 ns |
|---|---|---|---|---|---|---|---|---|
| V0 | **STOP 5.90 ns**: \|U\| 31.8 km/s in 4 cells, T at the 250 K floor | 704 | 0.31 (4.25 at the stop) | 6.7 | 45 / 45 | - | - | - |
| V1 | 20 ns | 1447 | 0.18 | 3.0 | 42 / 45 | 52 / 55 | 73 / 77 | 2111 |
| V2 | 20 ns | 1465 | 0.10 | 6.9 | 48 / 48 | 58 / 61 | 80 / 83 | 1787 |
| V3 | 20 ns | 2276 | 0.12 | 7.2 | 48 / 52 | 61 / 64 | 83 / 86 | 1931 |
| V4 | 20 ns | 1595 | 0.12 | 7.2 | 48 / 52 | 61 / 64 | 83 / 86 | 1931 |

- Fronts are measured in the first gas row: vapour where sum alpha_vap > 0.01, shock where p > 2e5.
- V0's blow-up is at (7.8, -42) um, at the plume top near the axis. It is the first time this 2D case reproduces
  the 3D-type failure: the mismatch climbs to 0.3 before it.
- The previous 2D runs had +eb and 5 ranks and did not include V0, so this was not seen before.
- Earlier V4-equivalent run (+eb, 5 ranks): 42 / 61 / 83 um at 2 / 10 / 20 ns, consistent with V4 here.

## V3 vs V4 (full vs contErr under mass)
They agree to round-off:

| case | byte-identical field files | worst max relative difference |
|---|---|---|
| A | 415 of 420 | 1e-15 |
| B | 370 of 420 | 5e-14 |
| Af | 385 of 420 | 1e-7, in U near zero |
| 2D | 18 of 140 | up to 8e-8 |

- Energy budgets are identical to the printed digits.
- This is expected: under mass, ddt(rhoM)+div(rhoPhi) is at the purge level, so the contErr correction term vanishes
  and the two modes differ only in round-off.

## Answers
- **(a) Is mass alone (V1) enough for stability?** Yes, in both 1D and 2D.
  - V0 fails in B (26.7 ns) and in 2D (5.9 ns); V1 runs every case to the end.
  - The mismatch stays bounded: 1D <= 0.14, 2D <= 0.18.
- **(b) What does each later variant add, and at what cost?**
  - **nc (V2):** a modest accuracy gain at no cost (same steps and iterations).
    - It halves the 2D mismatch (0.18 -> 0.10) and moves the 2D wall front from 73 to 80 um at 20 ns
      (the V3/V4 value is 83).
    - In 1D the gain is mixed: energy drift is better in A/Af (2.2 -> 1.6%, 2.3 -> 1.5%) but worse in B
      (1.4 -> 3.0%, with a hot shocked-air band of ~11.5 kK).
  - **o10 (V3):** the big accuracy step.
    - Energy drift falls to 0.2-0.4%.
    - A/Af p*, u* and T vap come back to within ~2% (u*) to ~5% (T vap) of exact, the shock spikes go away, and the
      2D fronts settle.
    - Cost: 5x the outer iterations. Measured wall time is ~2.5-4x in 1D (same-batch V1 vs V3: A 43 -> 168 s,
      Af 263 -> 1051 s, B 83 -> 209 s) and ~1.6x in 2D.
  - **contErr (V4):** no change in accuracy or cost under mass.
- **(c) Does V3 match V4 to round-off?** Yes (see above).

## Recommendation
- **Required for stability:** mass.
- **Required for accuracy:** o10.
- **Optional:** nc (helps the 2D plume, mixed in 1D B) and contErr (redundant under mass, harmless).

## Notes added 2 Nov 2026
- **Data:** tables, plots and CSVs are in `ablation_results/` next to this file; see its README.
- **2D vapour front depends on the threshold.** The vapour fronts above use sum alpha_vap > 0.01, the leading tail of a
  contact smeared numerically over about 18-27 cells (1-99% width of alpha.air in 1D; no compression on vapour/air).
  The small gap to the shock is probably mostly smearing, and the V1-V3 trend may partly reflect the tail. Shock
  positions (p > 2e5) are unaffected. `tools/plumeFront.py` uses the same 0.01 threshold. Not yet re-measured at 0.5.
- **B V0 contact velocity.** The table's 4560 m/s is summarise.py's u*. At alpha.air = 0.5 the contact moves at
  3.7-4.3 km/s from 10 to 26 ns (V1: 2.0-2.2 km/s). Max |U| reaches 32 km/s in the step that stops the run.
