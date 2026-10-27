{

  //https://www.cfd-online.com/Forums/openfoam-solving/251323-teqn-compressibleinterfoam.html





    int iter = 0;


    dimensionedScalar residual("residual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar meanResidual("meanResidual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar residual2("residual2",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar meanResidual2("meanResidual2",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));





    TRHS*=0.0;//=TRHSzero;

    // ---- Kinetic-energy / pressure-work terms: option ----------------------
    // PIMPLE { kineticEnergyTerms full | contErr | none; }
    //   full    (default) original: div(phi,p) + ddt(rho,K) + div(rhoPhi,K)
    //   contErr full minus contErr*K, contErr = ddt(rho) + div(rhoPhi), as in
    //           ESI compressibleInterFoam. With rho from the EOS and rhoPhi
    //           from the c_i transport, contErr is not ~0 where the volume
    //           closure is off, and the full form turns K*contErr into heat.
    //   none    internal-energy form: p*div(U) only (ESI totalInternalEnergy
    //           false). No kinetic-energy terms at all.
    //   dissipative
    //           contErr form split as p*div(U) + D, with
    //           D = U.grad(p) + rho DK/Dt (= the contErr work term - p*div(U)).
    //           Physically D = -(dissipation) <= 0, i.e. it can only heat.
    //           D > 0 (kinetic energy gained beyond the pressure work, e.g. a
    //           numerical velocity overshoot at a vapour/air contact) would
    //           cool the gas; it is dropped. D < 0 (shock heating) is kept.
    //           Work term = p*div(U) + min(D, 0).
    const word keMode
    (
        pimple.dict().getOrDefault<word>("kineticEnergyTerms", "full")
    );
    if
    (
        keMode != "full" && keMode != "contErr" && keMode != "none"
     && keMode != "dissipative"
    )
    {
        FatalErrorInFunction
            << "kineticEnergyTerms " << keMode
            << " unknown, use full, contErr, none or dissipative"
            << exit(FatalError);
    }

    const volScalarField contErrT
    (
        "contErrT",
        fvc::ddt(rho) + fvc::div(mixture.rhoPhi())
    );

    tmp<volScalarField> tWork;
    if (keMode == "none")
    {
        tWork = p*fvc::div(fvc::absolute(phi, U));
    }
    else
    {
        tWork =
            fvc::div(fvc::absolute(phi, U), p)
          + fvc::ddt(rho, K) + fvc::div(mixture.rhoPhi(), K);

        if (keMode == "contErr" || keMode == "dissipative")
        {
            tWork.ref() -= contErrT*K;
        }
    }

    // dissipative: keep only the heating part of D = work - p*div(U)
    tmp<volScalarField> tKeDrop;   // the dropped (cooling) part, >= 0
    const bool keDissipative = (keMode == "dissipative");
    if (keDissipative)
    {
        const volScalarField pDivU(p*fvc::div(fvc::absolute(phi, U)));
        const volScalarField D(tWork() - pDivU);
        const dimensionedScalar zeroD(D.dimensions(), Zero);
        tKeDrop = max(D, zeroD);
        tWork = pDivU + min(D, zeroD);
    }
    const volScalarField& workTerm = tWork();

    // diagnostic, any mode: the T change per step that K*contErr alone
    // would give (the part the full form puts into T and contErr removes)
    {
        const volScalarField rCvDiag(mixture.rCv());
        const scalarField dTce
        (
            mag(contErrT.primitiveField()*K.primitiveField())
           *rCvDiag.primitiveField()*runTime.deltaTValue()
           /max(rho.primitiveField(), scalar(1e-3))
        );
        scalar dMax = -1;
        label iMax = -1;
        forAll(dTce, celli)
        {
            if (dTce[celli] > dMax) { dMax = dTce[celli]; iMax = celli; }
        }
        const scalar gMaxV = returnReduce(dMax, maxOp<scalar>());
        Info<< "TEqn: kineticEnergyTerms " << keMode
            << ", mixtureCv " << mixture.mixtureCvMode()
            << ", max |K*contErr| as dT per step = " << gMaxV << " K";
        if (keDissipative)
        {
            const scalarField dropDT
            (
                tKeDrop().primitiveField()*rCvDiag.primitiveField()
               *runTime.deltaTValue()
               /max(rho.primitiveField(), scalar(1e-3))
            );
            Info<< ", max KE cooling dropped (dissipative) = "
                << gMax(dropDT) << " K";
        }
        Info<< endl;
        if (iMax >= 0 && dMax == gMaxV && gMaxV > 100)
        {
            Pout<< "TEqn:   at " << mesh.C()[iMax]
                << ", T = " << T[iMax] << ", rho = " << rho[iMax]
                << ", |U| = " << mag(U[iMax]) << endl;
        }
    }

    // no spurious mushy drag in gas
    const dimensionedScalar cMinMelt   ("cMinMelt",    dimless, 1e-4); // div guard
    const dimensionedScalar cThreshMelt("cThreshMelt", dimless, 1e-2); // mask ramp

    const volScalarField rCond(1.0/max(condensate, cMinMelt));
    const volScalarField LatentHeatCond(LatentHeat*rCond);  // -> metal L_f
    const volScalarField TSolidusCond  (TSolidus*rCond);    // -> metal T_sol
    const volScalarField TLiquidusCond (TLiquidus*rCond);   // -> metal T_liq

    // Apparent-heat-capacity coefficient
    // LatentHeat is now condensed-only (see update.H), so its  0 in
    // pure-gas cells where it used to be the gas dummy (~1).  
    const volScalarField meltCoeff
    (
        epsilonRel*mixture.Cv()
       /max(LatentHeatCond, dimensionedScalar("LHfloor", LatentHeatCond.dimensions(), 1.0))
    );

    // 1 where metal present, ramps to 0 in the gas
    const volScalarField meltMask
    (
        min(max(condensate/cThreshMelt, scalar(0)), scalar(1))
    );
    // ------------------------------------------------------------------------

    do{
iter++;
epsilon1.storePrevIter();



TRHS = meltMask*LatentHeat*
(
    fvc::ddt(rho, epsilon1)
  + fvc::div(mixture.rhoPhi(), epsilon1)
);

TRHS.correctBoundaryConditions();


    const volScalarField& Joulesource =
        mthd.valid()
      ? mthd->Joulesource()
      : zeroJoulesource;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(mixture.rhoPhi(), T)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(mixture.rhoPhi()), T)
      - fvm::laplacian(mixture.alphaEff(turbulence->mut()), T)
      + workTerm*mixture.rCv()

        ==
          mixture.rCv()*laser.deposition()
        + mixture.rCv()*Joulesource
        - mixture.rCv()*TRHS
        -mass_dot
    );

    TEqn.relax();
    TEqn.solve();



 {
        const dimensionedScalar Tfloor("Tfloor", dimTemperature, 250.0);
        // ceiling to match the floor, kinetic heating spikes have no other
        // bound. default 1e6 = off, set TMax in PIMPLE dict to arm it
        const scalar Tceil = pimple.dict().getOrDefault<scalar>("TMax", 1e6);
        scalarField& TI = T.primitiveFieldRef();
        label nClip = 0, nClipHi = 0;
        forAll(TI, celli)
        {
            if (TI[celli] < Tfloor.value()) { TI[celli] = Tfloor.value(); ++nClip; }
            if (TI[celli] > Tceil) { TI[celli] = Tceil; ++nClipHi; }
        }
        reduce(nClip, sumOp<label>());
        reduce(nClipHi, sumOp<label>());
        if (nClip) Info<< "TEqn: floored " << nClip << " cell(s)" << endl;
        if (nClipHi) Info<< "TEqn: ceilinged " << nClipHi << " cell(s)" << endl;
        T.correctBoundaryConditions();
    }


Tcorr = (TLiquidusCond - TSolidusCond)*epsilon1 + TSolidusCond;
Tcorr.correctBoundaryConditions();

// Masked, normalised melt-fraction update. In metal (meltMask=1) this is
// identical to the original; in (near-)gas (meltMask->0) epsilon1 -> 1.
epsilon1 =
    meltMask*max(min(epsilon1 + meltCoeff*(T - Tcorr), scalar(1)), scalar(0))
  + (scalar(1) - meltMask);


{
    const wordList e1MaskNames
    (
        mesh.lookupObject<IOdictionary>("thermophysicalProperties")
            .lookupOrDefault<wordList>("phaseChangeMaskPatches", wordList())
    );

    label nE1Mask = 0;
    scalarField& e1I = epsilon1.primitiveFieldRef();

    for (const word& pn : e1MaskNames)
    {
        const label patchi = mesh.boundaryMesh().findPatchID(pn);
        if (patchi < 0) continue;   // typo already warned about in solveAlphas

        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        forAll(fc, bfacei)
        {
            if (e1I[fc[bfacei]] < 1.0)
            {
                e1I[fc[bfacei]] = 1.0;
                ++nE1Mask;
            }
        }
    }
    reduce(nE1Mask, sumOp<label>());
    if (nE1Mask)
    {
        Info<< "solidification masked in " << nE1Mask
            << " boundary cell(s)" << endl;
    }
}
epsilon1.correctBoundaryConditions();


residual = gMax(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField()));
meanResidual = gSum(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField())*mesh.V())/gSum(mesh.V().field());


        Info << "Correcting epsilon1, mean residual = " << meanResidual.value()
             << ", max residual = " << residual.value()
             << endl;


}while((iter < minTCorr || (residual > epsilonTol)) && iter <= maxTCorr );


    // ---- Energy budget of one cell (DIAGNOSTICS { energyBudget true; }) ----
    // Each term of the T equation expressed as the temperature change it
    // gives over this step, evaluated explicitly with the solver's own fields
    // and schemes after the final-corrector T solve. Euler ddt with the
    // continuity Sp term gives rho0*(T - T0)/dt, so dT_term = term*dt/rho0.
    // "rest" = actual dT minus the sum: implicitness, relaxation, T clipping.
    if
    (
        healthLimits.energyBudget
     && (pimple.finalIter() || healthLimits.energyBudgetEveryIter)
    )
    {
        crashStage::set("TEqn: energy budget diagnostic");

        const scalar dtB = runTime.deltaTValue();
        const volScalarField& T0B = T.oldTime();

        // pick the cell
        label ci = -1;
        scalar key = -GREAT;
        if (healthLimits.energyBudgetHasPoint)
        {
            ci = mesh.findCell(healthLimits.energyBudgetPoint);
            if (ci >= 0) key = mag(T[ci] - T0B[ci]);
        }
        else
        {
            forAll(T, i)
            {
                const scalar d = mag(T[i] - T0B[i]);
                if (d > key) { key = d; ci = i; }
            }
        }
        const scalar gKey = returnReduce(key, maxOp<scalar>());

        if (gKey >= healthLimits.energyBudgetMinDT)
        {
            // collective field evaluations, all processors
            const volScalarField rCvB(mixture.rCv());
            const volScalarField CvMassB(mixture.massWeightedCv());
            const volScalarField rho0B(max(rho.oldTime(), dimensionedScalar(dimDensity, 1e-3)));

            const volScalarField advB
            (
                fvc::div(mixture.rhoPhi(), T) - T*fvc::div(mixture.rhoPhi())
            );
            const volScalarField condB
            (
                fvc::laplacian(mixture.alphaEff(turbulence->mut()), T)
            );
            const volScalarField pwB
            (
                (keMode == "none" || keMode == "dissipative")
              ? volScalarField(p*fvc::div(fvc::absolute(phi, U)))
              : volScalarField(fvc::div(fvc::absolute(phi, U), p))
            );
            const volScalarField keB
            (
                fvc::ddt(rho, K) + fvc::div(mixture.rhoPhi(), K)
            );
            const volScalarField& JouleB =
                mthd.valid() ? mthd->Joulesource() : zeroJoulesource;
            const volScalarField& depB = laser.deposition();

            if (ci >= 0 && key == gKey)
            {
                const scalar f = dtB/rho0B[ci];
                const scalar rc = rCvB[ci];

                const scalar dAdv  = -advB[ci]*f;
                const scalar dCond =  condB[ci]*f;
                const scalar dPw   = -pwB[ci]*rc*f;
                const scalar dCeK  = contErrT[ci]*K[ci]*rc*f;   // +: removed by contErr
                scalar dKe = (keMode == "none") ? 0 : -keB[ci]*rc*f;
                scalar dCeKApplied = (keMode == "contErr") ? dCeK : 0;
                scalar dDrop = 0;
                if (keMode == "dissipative")
                {
                    // applied work = p*div(U) + min(D,0); report the D part
                    // as "kinetic energy" and the dropped cooling separately
                    dKe = -(workTerm[ci] - pwB[ci])*rc*f;
                    dCeKApplied = 0;
                    dDrop = tKeDrop()[ci]*rc*f;
                }
                const scalar dLas  =  depB[ci]*rc*f;
                const scalar dJou  =  JouleB[ci]*rc*f;
                const scalar dMelt = -TRHS[ci]*rc*f;
                const scalar dPc   = -mass_dot[ci]*f;

                const scalar dT = T[ci] - T0B[ci];
                const scalar sum =
                    dAdv + dCond + dPw + dKe + dCeKApplied
                  + dLas + dJou + dMelt + dPc;

                OStringStream ph;
                for (const auto& phs : mixture.phases())
                {
                    if (phs[ci] > 1e-3)
                    {
                        ph  << ' ' << phs.name() << '=' << phs[ci];
                    }
                }

                Pout<< "[budget] t = " << runTime.timeName()
                    << ", outer " << pimple.corr()
                    << ", cell at " << mesh.C()[ci]
                    << (healthLimits.energyBudgetHasPoint ? " (fixed point)" : " (largest |dT|)") << nl
                    << "[budget]   T " << T0B[ci] << " -> " << T[ci]
                    << "  (dT " << dT << " K)" << nl
                    << "[budget]   advection of T        " << dAdv << nl
                    << "[budget]   conduction            " << dCond << nl
                    << "[budget]   pressure work         " << dPw
                    << (keMode == "none" ? "   (p*div(U))" : "   (div(phi,p))") << nl
                    << "[budget]   kinetic energy        " << dKe
                    << (keMode == "none" ? "   (off)" : "")
                    << (keMode == "dissipative" ? "   (heating part of D only)" : "") << nl
                    << "[budget]   KE cooling dropped    " << dDrop
                    << (keMode == "dissipative" ? "   (not applied)" : "   (n/a)") << nl
                    << "[budget]   K*contErr             " << dCeKApplied
                    << "   (size " << dCeK << ", applied only in contErr mode)" << nl
                    << "[budget]   laser                 " << dLas << nl
                    << "[budget]   Joule                 " << dJou << nl
                    << "[budget]   melting (TRHS)        " << dMelt << nl
                    << "[budget]   phase change (mass_dot) " << dPc << nl
                    << "[budget]   sum " << sum << ", rest (implicit/relax/clip) "
                    << dT - sum << nl
                    << "[budget]   rho " << rho0B[ci] << " -> " << rho[ci]
                    << ", p " << p[ci] << ", |U| " << mag(U[ci])
                    << ", K " << K[ci] << ", contErr " << contErrT[ci] << nl
                    << "[budget]   Cv used 1/rCv = " << 1.0/max(rc, VSMALL)
                    << " (mixtureCv " << mixture.mixtureCvMode()
                    << "), mass-weighted Cv = " << CvMassB[ci] << nl
                    << "[budget]   phases:" << ph.str().c_str() << endl;
            }
        }
    }



    mixture.correct();
}