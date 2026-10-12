{
    // Temperature-dependent electrical resistivity of the metal, used by the
    // laser ray tracer for the Fresnel absorptivity in interface cells.
    //
    // The gas has no free electrons, so the resistivity the laser should see
    // in a mixed (interface) cell is that of the metal. Weighting by
    // alpha.metal lowers rho_e - and hence the absorptivity - by an amount
    // that depends on where the interface sits in the cell (mesh dependent).
    // It is kept only as an option (resistivityAlphaWeighting yes) to
    // reproduce the original PR #140 behaviour.

    scalarField& eResI = electrical_resistivity;
    const scalarField& TI = T;
    const scalarField& alpha1I = alpha1;

    forAll(eResI, celli)
    {
        const scalar rhoE = electricalResistivity_m->value(TI[celli]);

        if (resistivityAlphaWeighting)
        {
            eResI[celli] = max(min(alpha1I[celli], 1.0), 0.0)*rhoE;
        }
        else
        {
            eResI[celli] = rhoE;
        }
    }

    electrical_resistivity.correctBoundaryConditions();

    Info<< "electrical_resistivity min/max: "
        << min(electrical_resistivity).value() << " "
        << max(electrical_resistivity).value() << nl;
}


// end
