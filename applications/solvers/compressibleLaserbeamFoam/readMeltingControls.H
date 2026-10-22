IOdictionary solutionDict
(
    IOobject
    (
        "fvSolution",
        runTime.system(), //folder where the fvSolution file is located
        runTime,
        IOobject::MUST_READ
    )
);



dictionary& meltingDict = solutionDict.subDict("MELTING");

int minTCorr(readLabel(meltingDict.lookup("minTempCorrector")));
int maxTCorr(readLabel(meltingDict.lookup("maxTempCorrector")));

dimensionedScalar epsilonTol(readScalar(meltingDict.lookup("epsilonTolerance")));
dimensionedScalar epsilonRel(readScalar(meltingDict.lookup("epsilonRelaxation")));

// Optional time-dependent multiplier on the Darcy (Carman-Kozeny) momentum
// sink, e.g. to treat the condensed metal as a compressible fluid during a
// ns pulse and restore the mushy-zone damping afterwards:
//
//     MELTING
//     {
//         ...
//         darcyScale table ((0 0) (5.1e-9 0) (6e-9 1));
//     }
//
// Absent -> scale = 1 (original behaviour). Re-read every time step.
autoPtr<Function1<scalar>> darcyScalePtr
(
    Function1<scalar>::NewIfPresent("darcyScale", meltingDict)
);

// Field health checks (optional DIAGNOSTICS dict in fvSolution, see
// solverHealth.H). Re-read every time step.
const Foam::solverHealth::Limits healthLimits
(
    Foam::solverHealth::readLimits(solutionDict)
);
