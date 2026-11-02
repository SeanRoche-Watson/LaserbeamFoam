# Ablation V0-V4 data
Summary: ../RESULTS_ablation.md. Raw runs are not in git.

- `ablation_tables.txt`, `energy_table.txt`: tables and energy_check.py output.
- `profiles_{A,Af,B}.png`, `2d/`: plots and `2d/fronts.csv`.
- `csv/log_<case>_<V>.csv`: per step from log.solver (max U, min T, mismatch, continuity error).
- `csv/contact_<case>_<V>.csv`: 1D per write: contact position and velocity, max |U|, min T.
- `csv/fields_P_<V>.csv`: 2D per write, from the fields.
- `csv/profile_<case>_<V>_t*.csv`: 1D profiles at 60 ns, with exact for A/Af. B V0 is at 26 ns, its last write before the stop.
- Rows noted "health dump at stop" are the fields at the failure.

Caveat: the 2D vapour front uses vapour fraction > 0.01, the tail of the numerically smeared contact, so its gap to
the shock is probably mostly smearing. Shock fronts (p > 2e5) are unaffected.

No guard log lines at ~1.9 ns are included: none were found in the local logs searched.
