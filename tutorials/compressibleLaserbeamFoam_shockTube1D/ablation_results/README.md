# Ablation V0-V4: data behind RESULTS_ablation.md

The raw runs (`runs/`, ~2.7 GB) are gitignored and stay on the WSL machine (`~/shockTube_mm/runs`, `~/plumeEdge2D/runs`).
This folder holds the tables, plots and small CSVs extracted from them.

Variants: V0 `full`, V1 `full:mm`, V2 `full:nc+mm`, V3 `full:o10+nc+mm`, V4 `contErr:o10+nc+mm`.
Cases: A, Af (case A at 4x resolution), B, and P (2D plume edge).

## Tables and plots
- `ablation_tables.txt`: per-case tables, V3-vs-V4 comparison, 2D fronts.
- `energy_table.txt`: full energy_check.py output for the 1D runs.
- `profiles_A.png`, `profiles_Af.png`, `profiles_B.png`: V0-V4 at 60 ns (B V0 at 26 ns).
- `2d/fields_V0.png` ... `2d/fields_V4.png`, `2d/fronts.png`, `2d/fronts.csv`.
- `scripts/`: the scripts that made all of this.

## CSVs (`csv/`)
| file | content | use it for |
|---|---|---|
| `contact_<A,Af,B>_<V>.csv` | per 1 ns write: contact position and velocity at alpha.air = 0.5, max \|U\|, min T overall and in the vapour | B V0 vs V1: the runaway and the fix |
| `fields_P_<V>.csv` | per 1 ns write, from the 2D fields: max \|U\|, min T overall and in the vapour, max T, max p, max \|rho/rhoM - 1\| | 2D V0 vs V1 |
| `series_<A,Af,B,P>_<V>.csv` | per time step, from log.solver: deltaT, max(U), min(T), mismatch, continuity error | time histories between writes |
| `profile_<A,Af,B>_<V>_t<time>.csv` | sampled profile at 60 ns (B V0 at 26 ns); A and Af include the exact Riemann solution | A/Af V0 vs V3 against exact |

Reading them:
- **Use `minT_vapour_K` for the cold dip.** `minT_K` and the solver's `min(T)` line are tube-wide, so they stay at the
  300 K air until the 250 K floor is hit at the stop. Vapour means alpha.air < 0.5 in 1D and sum alpha_vap > 0.5 in 2D.
- **Rows noted `health dump at stop`** are the fields written at the failure, not a normal write.
- **In `series_*`,** `maxU_final` is the last max(U) of the step and `maxU_step` the largest over its correctors.
- **`dx_contact_dt_m_s`** comes from the written positions and is smooth; it should track `U_contact_m_s`.

## What the CSVs show
| 1D B | t | contact U (alpha.air 0.5) | max \|U\| | min T vapour |
|---|---|---|---|---|
| V0 | 10 ns | 4294 m/s | 4338 | 1755 K |
| V0 | 26 ns | 3715 | 5140 | 1571 |
| V0 | 26.7 ns (stop) | 3782 | 31971 | 250 (floor) |
| V1 | 10 ns | 2232 | 2277 | 6814 |
| V1 | 26 ns | 2076 | 2323 | 6076 |
| V1 | 60 ns | 2044 | 2384 | 5769 |

| 2D P | t | max \|U\| | min T vapour | max \|rho/rhoM-1\| |
|---|---|---|---|---|
| V0 | 2 ns | 2495 m/s | 8507 K | 0.101 |
| V0 | 5 ns | 3168 | 3314 | 0.100 |
| V0 | 5.9 ns (stop) | 31822 | 250 (floor) | 4.25 |
| V1 | 2 ns | 2041 | 9309 | 0.043 |
| V1 | 5 ns | 2187 | 6951 | 0.084 |
| V1 | 20 ns | 2111 | 5331 | 0.104 |

## Guard log lines at ~1.9 ns: not found
No log line from a guard firing at about 1.9 ns exists in any of this data.
- **Searched:** all 540 `log.*` files in the WSL home (excluding paraview and the Python environments), case-insensitive
  for "guard" and including logs grep treats as binary; the docs in this repo; the shockTube_work text files.
  Zero matches.
- **`[health]` lines:** none before 3 ns in any log. In these runs they appear only at the stops
  (B V0 26.7 ns, P V0 5.90 ns).
- **"guard" in the solver source:** only `cMinMelt` (TEqn.H) and `coarsenGuard` (update.H). Neither prints anything.
- **Not searched:** the logs on CSD3, which would be the place to look if the guard fired in a 3D run.
- **Not a substitute:** 2D V0 does print `group flux limiter: gasPackingCap on ...` around 1.9 ns, but that line fires
  on nearly every step from about 0.01 ns, so it does not mark an event at 1.9 ns.
