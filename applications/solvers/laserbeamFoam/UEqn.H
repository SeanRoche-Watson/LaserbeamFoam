MRF.correctBoundaryVelocity(U);


if (laser.powderSim())
{
    DC =
        DarcyConstantlarge*Foam::pow((1.0-epsilon1mask),2)
       /(Foam::pow(epsilon1mask, 3) + DarcyConstantsmall);
}
else
{
    DC =
        DarcyConstantlarge*Foam::pow((1.0 - epsilon1), 2)
       /(Foam::pow(epsilon1, 3) + DarcyConstantsmall);
}

DC.correctBoundaryConditions();  // ADDED

gradAlpha = fvc::grad(alpha1);
gradAlpha.correctBoundaryConditions();  // ADDED - CRITICAL!

nHatM = gradAlpha/(mag(gradAlpha) + deltaN);
nHatM.correctBoundaryConditions();  // ADDED

const scalar Tliquidus1Value = Tliquidus1.value();
const scalarField& TI = T;

if (sahooMarangoni)
{
    // Only molten cells carry a Marangoni force (zeroed below), so clamp the
    // evaluation temperature at the liquidus to keep K finite elsewhere.
    scalarField& dsI = dsigmadT.primitiveFieldRef();
    forAll(dsI, celli)
    {
        dsI[celli] = sahooDsigmadT(max(TI[celli], Tliquidus1Value));
    }
    dsigmadT.correctBoundaryConditions();
}

Marangoni =
    epsilon1*dsigmadT*(gradT - nHatM*(nHatM & gradT))*mag(gradAlpha);

vectorField& MarangoniI = Marangoni;
forAll(MarangoniI, celli)
{
    if(TI[celli] < Tliquidus1Value)
    {
        MarangoniI[celli] = vector::zero;
    }
}

Marangoni.correctBoundaryConditions();

// Diagnostic: are the gradients behind the Marangoni force being limited away?
// Reference = unlimited Gauss linear, whatever gradSchemes says.
// "retained" = sum|grad used| / sum|grad reference| over molten interface cells
// (1 = no limiting, ~0 = the limiter has removed it); "zeroed" = fraction of those
// cells where the used gradient is < 1% of the reference.
if (gradientDiagnostics)
{
    const fv::gaussGrad<scalar> unlimitedGrad(mesh);
    const volVectorField gradAlphaRef(unlimitedGrad.calcGrad(alpha1, "gradAlphaRef"));
    const volVectorField gradTRef(unlimitedGrad.calcGrad(T, "gradTRef"));

    scalar usedA = 0, refA = 0, usedT = 0, refT = 0;
    label nCells = 0, nZeroA = 0;
    const scalarField& V = mesh.V();

    forAll(TI, celli)
    {
        if (TI[celli] < Tliquidus1Value) continue;

        // interface cell: alpha changes by > 0.1 across the cell
        const scalar gA = mag(gradAlphaRef[celli]);
        if (gA*Foam::cbrt(V[celli]) < 0.1) continue;

        const vector n = gradAlphaRef[celli]/gA;
        const scalar uA = mag(gradAlpha[celli]);

        usedA += uA;
        refA  += gA;
        usedT += mag(gradT[celli] - n*(n & gradT[celli]));          // tangential part
        refT  += mag(gradTRef[celli] - n*(n & gradTRef[celli]));
        ++nCells;
        if (uA < 0.01*gA) ++nZeroA;
    }

    reduce(usedA, sumOp<scalar>());
    reduce(refA,  sumOp<scalar>());
    reduce(usedT, sumOp<scalar>());
    reduce(refT,  sumOp<scalar>());
    reduce(nCells, sumOp<label>());
    reduce(nZeroA, sumOp<label>());

    Info<< "Marangoni gradient check (molten interface cells " << nCells << "):"
        << " grad(alpha) retained " << usedA/(refA + VSMALL)
        << "  zeroed " << scalar(nZeroA)/max(nCells, label(1))
        << "  | tangential grad(T) retained " << usedT/(refT + VSMALL)
        << endl;
}

if (damperSwitch)
{
    // To account for high density difference
    damper = 2.0*rho/(rho1 + rho2);
    damper.correctBoundaryConditions();  // ADDED
}

pVap = 0.54*p0*Foam::exp(LatentHeatVap*Mm*((T - Tvap)/(R*T*Tvap)));
pVap.correctBoundaryConditions();  // ADDED


if (mthd.valid())
{
    mthd->updateLorentzForce();
}

const volVectorField& JcrossB =
    mthd.valid()
  ? mthd->JcrossB()
  : zeroJcrossB;

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
  + fvm::Sp(DC, U)
  - Marangoni*damper
  - JcrossB*damper
 ==
    fvOptions(rho, U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()*fvc::interpolate(damper)
                + (fvc::interpolate(pVap)
                    *fvc::snGrad(alpha1)
                    *fvc::interpolate(damper))
              - ghf*fvc::snGrad(rho*rhok)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
