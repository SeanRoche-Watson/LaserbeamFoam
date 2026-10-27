# 1D vapour/air shock tube for compressibleLaserbeamFoam

Investigation case for the LIBS plume failures (see CONTEXT.md for the full history and results).

    source <OpenFOAM v2506 bashrc>
    ./Allrun                       # A and B in full/contErr/none
    ./Allrun B:contErr:o10+eb      # case:mode:variants (see README.txt and Allrun header)
    python3 summarise.py runs/*    # table vs the exact Riemann solution (case A)
    python3 energy_check.py runs/* # mass/energy conservation (closed tube)

Notes: needs libROUNDSchemes (fvSchemes uses ROUNDAplus01 for beta; set e.g. Gauss vanLeer01 if absent).
Options gasPackingCap / mixtureCv live in libmultiphaseVapMixtureThermo.so, so `strings` on the solver binary
will not show them. Set gasPackingCap false in constant/thermophysicalProperties for the recommended settings.

Variants nc (gasPackingCap false), mm (momentumDensity mass, cc-mass-momentum build) and vl (vanLeer01 alpha
scheme, no libROUNDSchemes), e.g. `./Allrun A:contErr:o10+nc+vl A:contErr:o10+nc+vl+mm`.
