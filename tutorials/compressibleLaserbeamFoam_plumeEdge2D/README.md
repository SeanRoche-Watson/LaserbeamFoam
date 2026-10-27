# Planar 2D plume edge for compressibleLaserbeamFoam

A half-disc of Fe vapour (r = 30 um, 1000 bar, 12000 K; the 3D plume state at ~5 ns) sits on a no-slip,
adiabatic wall at y = 0. It expands into the "air" phase (argon properties) at 1 bar, 300 K. The gas is at
y < 0, as in the 3D cases.
- **Mesh:** 160 x 100 cells of 3.125 um (the finest 3D AMR level), 500 x 312.5 um, closed walls.
- **Unlike the 1D shock tube,** the flow turns along the wall at the plume edge, where the 3D runs fail. The
  EOS density and the transported mass (sum c_i) are free to disagree there.
- **Settings:** template constant/ and system/ (schemes, PIMPLE, DIAGNOSTICS) are those of
  ../compressibleLaserbeamFoam_shockTube1D/template. No laser and no AMR.

    source <OpenFOAM v2506 bashrc>
    NP=5 ./Allrun P:contErr:o10+nc+eb P:contErr:o10+nc+eb+mm     # case:mode:variants, see the Allrun header
    END=2e-8 NP=5 ./Allrun ...                                    # end time (default 6e-8), write interval WRITE (2e-9)

Post-processing (needs numpy and matplotlib):

    (cd runs/<run> && postProcess -funcs '(writeCellCentres writeCellVolumes)')
    python3 plumeEdge.py <outdir> eos=runs/P_contErr_o10+nc+eb mass=runs/P_contErr_o10+nc+eb+mm
    python3 ../tools/plumeFront.py runs/<run> --El 0.21            # fronts CSV + log-log plot

Cost (WSL, 12-core laptop): about 480 s per simulated ns on 5 ranks with o10 at maxDeltaT 5e-11, so
20 ns takes ~2.7 h.

Note: the repository .gitignore ignores `0*`, so `template/0.orig` is added with `git add -f`.
