/*---------------------------------------------------------------------------*\
Description
    Per-stage field health checks for compressibleLaserbeamFoam.

    At each checkpoint in the time loop (after mesh update, after the alpha /
    phase-change solve, after TEqn, after pEqn) the main fields are scanned:

      T, p, p_rgh, rho, |U|               non-finite and out-of-range
      alpha.<phase>                        non-finite and outside [-tol, 1+tol]
      rho of each phase thermo             non-finite and <= rhoMin
      c.<phase> (partial mass)             non-finite and < -cTol

    Counts are reduced over all processors and printed with Info. The
    processor holding the worst value prints its location with Pout, so the
    log says *where* it is, e.g.

      [health] after TEqn: T above 1e+05 in 12 cell(s); min 299.998 max 1.2e+05
      [125] [health]   max T = 1.2e+05 at (1.56e-06 1.56e-06 -1.56e-06)

    If anything is non-finite, all fields are written at the current time
    (writeOnFailure) and the run stops with a FatalError naming the stage
    (abortOnNonFinite), instead of an FPE some steps later with no fields.

    Controlled from system/fvSolution (re-read every time step, all optional):

        DIAGNOSTICS
        {
            enabled           true;   // default true
            checkEvery        1;      // time steps between checks
            verbose           false;  // print min/max even when healthy
            TMin              1;
            TMax              1e6;
            pMin              -1e12;
            pMax              1e12;
            UMax              1e6;
            rhoMin            0;
            alphaTol          1e-6;
            cTol              1e-6;   // kg/m3
            writeOnFailure    true;   // write fields when non-finite found
            abortOnNonFinite  true;
            stopOnViolation   false;  // also write+stop on range violations
        }

\*---------------------------------------------------------------------------*/

#ifndef solverHealth_H
#define solverHealth_H

#include "fvCFD.H"
#include "crashStage.H"
#include <cmath>

namespace Foam
{
namespace solverHealth
{

struct Limits
{
    bool enabled = true;
    label checkEvery = 1;
    bool verbose = false;
    scalar TMin = 1;
    scalar TMax = 1e6;
    scalar pMin = -1e12;
    scalar pMax = 1e12;
    scalar UMax = 1e6;
    scalar rhoMin = 0;
    scalar alphaTol = 1e-6;
    scalar cTol = 1e-6;
    bool writeOnFailure = true;
    bool abortOnNonFinite = true;
    bool stopOnViolation = false;
};


inline Limits readLimits(const dictionary& solutionDict)
{
    Limits L;
    const dictionary* dPtr = solutionDict.findDict("DIAGNOSTICS");
    if (!dPtr)
    {
        return L;
    }
    const dictionary& d = *dPtr;

    L.enabled          = d.getOrDefault<bool>("enabled", L.enabled);
    L.checkEvery       = max(label(1), d.getOrDefault<label>("checkEvery", L.checkEvery));
    L.verbose          = d.getOrDefault<bool>("verbose", L.verbose);
    L.TMin             = d.getOrDefault<scalar>("TMin", L.TMin);
    L.TMax             = d.getOrDefault<scalar>("TMax", L.TMax);
    L.pMin             = d.getOrDefault<scalar>("pMin", L.pMin);
    L.pMax             = d.getOrDefault<scalar>("pMax", L.pMax);
    L.UMax             = d.getOrDefault<scalar>("UMax", L.UMax);
    L.rhoMin           = d.getOrDefault<scalar>("rhoMin", L.rhoMin);
    L.alphaTol         = d.getOrDefault<scalar>("alphaTol", L.alphaTol);
    L.cTol             = d.getOrDefault<scalar>("cTol", L.cTol);
    L.writeOnFailure   = d.getOrDefault<bool>("writeOnFailure", L.writeOnFailure);
    L.abortOnNonFinite = d.getOrDefault<bool>("abortOnNonFinite", L.abortOnNonFinite);
    L.stopOnViolation  = d.getOrDefault<bool>("stopOnViolation", L.stopOnViolation);

    return L;
}


//- Result of checking one field
struct Result
{
    label nNonFinite = 0;
    label nViolation = 0;
};


//- Check a list of cell values (internal field) and boundary values.
//  'name' is what is printed. Collective: call on all processors.
inline Result checkValues
(
    const word& stage,
    const word& name,
    const scalarField& v,
    const label nBoundaryNonFinite,
    const vectorField& C,
    const scalar lo,
    const scalar hi,
    const bool verbose
)
{
    label nBad = 0, nLo = 0, nHi = 0;
    label iBad = -1, iMin = -1, iMax = -1;
    scalar vMin = GREAT, vMax = -GREAT;

    forAll(v, i)
    {
        const scalar x = v[i];
        if (!std::isfinite(x))
        {
            if (iBad < 0) iBad = i;
            ++nBad;
            continue;
        }
        if (x < vMin) { vMin = x; iMin = i; }
        if (x > vMax) { vMax = x; iMax = i; }
        if (x < lo) ++nLo;
        if (x > hi) ++nHi;
    }

    label nBadB = nBoundaryNonFinite;

    const scalar gMinV = returnReduce(vMin, minOp<scalar>());
    const scalar gMaxV = returnReduce(vMax, maxOp<scalar>());

    const label nBadLocal = nBad;
    reduce(nBad, sumOp<label>());
    reduce(nBadB, sumOp<label>());
    reduce(nLo, sumOp<label>());
    reduce(nHi, sumOp<label>());

    if (nBad || nBadB)
    {
        Info<< "[health] " << stage << ": " << name << " NON-FINITE in "
            << nBad << " cell(s) and " << nBadB << " boundary face(s)"
            << endl;

        if (nBadLocal)
        {
            Pout<< "[health]   first non-finite " << name << " on this "
                << "processor at " << C[iBad] << " (" << nBadLocal
                << " here)" << endl;
        }
    }

    if (nLo)
    {
        Info<< "[health] " << stage << ": " << name << " below " << lo
            << " in " << nLo << " cell(s); global min " << gMinV << endl;
        if (iMin >= 0 && vMin == gMinV)
        {
            Pout<< "[health]   min " << name << " = " << vMin << " at "
                << C[iMin] << endl;
        }
    }

    if (nHi)
    {
        Info<< "[health] " << stage << ": " << name << " above " << hi
            << " in " << nHi << " cell(s); global max " << gMaxV << endl;
        if (iMax >= 0 && vMax == gMaxV)
        {
            Pout<< "[health]   max " << name << " = " << vMax << " at "
                << C[iMax] << endl;
        }
    }

    if (verbose && !nBad && !nBadB && !nLo && !nHi)
    {
        Info<< "[health] " << stage << ": " << name << " ok, min "
            << gMinV << " max " << gMaxV << endl;
    }

    Result r;
    r.nNonFinite = nBad + nBadB;
    r.nViolation = nLo + nHi;
    return r;
}


inline label countBoundaryNonFinite(const volScalarField& f)
{
    label n = 0;
    forAll(f.boundaryField(), patchi)
    {
        const scalarField& pf = f.boundaryField()[patchi];
        forAll(pf, facei)
        {
            if (!std::isfinite(pf[facei])) ++n;
        }
    }
    return n;
}


inline Result checkScalar
(
    const word& stage,
    const word& name,
    const volScalarField& f,
    const scalar lo,
    const scalar hi,
    const bool verbose
)
{
    return checkValues
    (
        stage, name, f.primitiveField(), countBoundaryNonFinite(f),
        f.mesh().C().primitiveField(), lo, hi, verbose
    );
}


inline Result checkVector
(
    const word& stage,
    const word& name,
    const volVectorField& f,
    const scalar hi,
    const bool verbose
)
{
    const vectorField& vf = f.primitiveField();
    scalarField m(vf.size());
    forAll(vf, i)
    {
        const vector& u = vf[i];
        m[i] =
            (std::isfinite(u.x()) && std::isfinite(u.y()) && std::isfinite(u.z()))
          ? Foam::mag(u)
          : std::nan("");
    }

    label nB = 0;
    forAll(f.boundaryField(), patchi)
    {
        const vectorField& pf = f.boundaryField()[patchi];
        forAll(pf, facei)
        {
            const vector& u = pf[facei];
            if
            (
                !std::isfinite(u.x())
             || !std::isfinite(u.y())
             || !std::isfinite(u.z())
            )
            {
                ++nB;
            }
        }
    }

    return checkValues
    (
        stage, name, m, nB, f.mesh().C().primitiveField(), -GREAT, hi, verbose
    );
}


//- Full check of the solver state at a named stage. Collective.
template<class MixtureType>
inline void check
(
    const word& stage,
    const Limits& L,
    Time& runTime,
    const fvMesh& mesh,
    const volScalarField& T,
    const volScalarField& p,
    const volScalarField& p_rgh,
    const volScalarField& rho,
    const volVectorField& U,
    const MixtureType& mixture
)
{
    if (!L.enabled || (runTime.timeIndex() % L.checkEvery) != 0)
    {
        return;
    }

    crashStage::set(std::string("health check: ") + stage);

    label nNonFinite = 0;
    label nViolation = 0;

    auto add = [&](const Result& r)
    {
        nNonFinite += r.nNonFinite;
        nViolation += r.nViolation;
    };

    add(checkScalar(stage, "T", T, L.TMin, L.TMax, L.verbose));
    add(checkScalar(stage, "p", p, L.pMin, L.pMax, L.verbose));
    add(checkScalar(stage, "p_rgh", p_rgh, L.pMin, L.pMax, L.verbose));
    add(checkScalar(stage, "rho", rho, L.rhoMin, GREAT, L.verbose));
    add(checkVector(stage, "mag(U)", U, L.UMax, L.verbose));

    for (const auto& ph : mixture.phases())
    {
        const word phName(ph.name());

        add
        (
            checkScalar
            (
                stage, word("alpha." + phName), ph,
                -L.alphaTol, 1 + L.alphaTol, L.verbose
            )
        );

        {
            const tmp<volScalarField> tRhoPh(ph.thermo().rho());
            add
            (
                checkScalar
                (
                    stage, word("rho." + phName), tRhoPh(),
                    L.rhoMin, GREAT, L.verbose
                )
            );
        }

        const volScalarField* cPtr =
            mesh.findObject<volScalarField>(word("c." + phName));
        if (cPtr)
        {
            add
            (
                checkScalar
                (
                    stage, word("c." + phName), *cPtr,
                    -L.cTol, GREAT, L.verbose
                )
            );
        }
    }

    // already globally reduced inside checkValues

    const bool fail =
        nNonFinite > 0 || (L.stopOnViolation && nViolation > 0);

    if (fail)
    {
        Info<< nl << "[health] *** PROBLEM DETECTED at stage '" << stage
            << "', time " << runTime.timeName() << ": " << nNonFinite
            << " non-finite value(s), " << nViolation
            << " range violation(s) ***" << endl;

        if (L.writeOnFailure)
        {
            Info<< "[health] writing all fields at time "
                << runTime.timeName() << " for inspection" << endl;
            crashStage::set(std::string("health check: writing fields after failure at ") + stage);
            runTime.writeNow();
        }

        if (L.abortOnNonFinite || L.stopOnViolation)
        {
            FatalErrorInFunction
                << "Solver state is invalid at stage '" << stage
                << "', time " << runTime.timeName() << nl
                << "    non-finite values: " << nNonFinite << nl
                << "    range violations:  " << nViolation << nl
                << "See the [health] lines above for fields and locations."
                << exit(FatalError);
        }
    }
}

} // End namespace solverHealth
} // End namespace Foam

#endif
