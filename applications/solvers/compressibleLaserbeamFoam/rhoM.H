// rhoM = max(sum_i c_i, rhoMFloor), boundary values included
{
    volScalarField sumC
    (
        IOobject("sumC", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimDensity, Zero)
    );
    for (const phaseModel& ph : mixture.phases())
    {
        sumC += mesh.lookupObject<volScalarField>("c." + ph.name());
    }
    rhoM = max(sumC, rhoMFloor);
}
