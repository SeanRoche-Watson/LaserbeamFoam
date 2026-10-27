

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

if (darcyScalePtr)
{
    const scalar darcyScale =
        Foam::max(darcyScalePtr->value(runTime.value()), scalar(0));

    DC *= darcyScale;

    Info<< "Darcy scale = " << darcyScale << endl;
}


Marangoni = condensate*epsilon1*mixture.MarangoniForce(T);

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
    fvm::ddt(rhoMom, U)   // rho (eos) or rhoM (mass), see createFields
  + fvm::div(mixture.rhoPhi(), U)
  + turbulence->divDevRhoReff(U)
  + fvm::Sp(DC, U)
  + Marangoni
  - JcrossB
);


UEqn.relax();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rhoMom)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    K = 0.5*magSqr(U);
}
