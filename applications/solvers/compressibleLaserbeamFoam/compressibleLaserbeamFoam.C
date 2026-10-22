/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    compressibleLaserbeamFoam

Group
    grpMultiphaseSolvers

Description
    cheeky bad solver for wicked simulation of laser-substrate interactions. A true naughty piece of work

Authors
    Tom Flint UoM
    Roman Bialek UoM
    Philip Cardiff UcD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "multiphaseMixtureThermo.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "Function1.H"
#include "crashStage.H"
#include "solverHealth.H"
#include "CorrectPhi.H"

#include "laserHeatSource.H"
#include "mthdModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for N compressible, non-isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"

    #include "createDynamicFvMesh.H"

    // #include "createMesh.H"
    // #include "createControl.H"
    // #include "createTimeControls.H"

    #include "initContinuityErrs.H"//new
    #include "createDyMControls.H"//new

    // bool bRestartFirstLoop = false; <- declared within createFields.H
    #include "createFields.H"

    // Crash reporting: on SIGFPE/SIGSEGV each processor prints the last
    // stage it reached (see crashStage.H) before OpenFOAM's stack trace.
    crashStage::install(Pstream::myProcNo());
    Info<< "crashDiagnostics: signal reporter installed "
        << "(stage markers + DIAGNOSTICS health checks)" << endl;

    // Field health check at a named stage; see solverHealth.H
    #define HEALTH_CHECK(stageName)                                           \
        solverHealth::check                                                   \
        (                                                                     \
            stageName, healthLimits, runTime, mesh, T, p, p_rgh, rho, U,      \
            mixture                                                           \
        )

    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"



    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();

    turbulence->validate();

    #include "update.H"


    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

        crashStage::set("time loop: reading controls / setDeltaT");

        #include "readMeltingControls.H"
        #include "readDyMControls.H"

        // #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        crashStage::set("time " + runTime.timeName() + ": start of step");

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {

            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                crashStage::set
                (
                    "time " + runTime.timeName()
                  + ": mesh.update (AMR refine/unrefine/balance, field mapping)"
                );
                mesh.update();

                if (mesh.changing())
                {
                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi && !bRestartFirstLoop)
                    {
                        // Prevent corrections on the first loop of a restarted
                        // solver.

                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.esi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    
                    crashStage::set
                    (
                        "time " + runTime.timeName()
                      + ": mixture.correct() after mesh change"
                    );
                    mixture.correct();

                    // AMR remapping is a recurring crash suspect: check the
                    // freshly mapped fields before anything uses them
                    HEALTH_CHECK("after mesh.update (mesh changed)");

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // Disable flag: this permits correctPhi on subsequent loops.
            bRestartFirstLoop = false;

            // Semi-PIMPLE Implementation - Roman Bialek
            if (pimple.firstIter())
            {
                rho.oldTime() = rho;   // set once per time step

                HEALTH_CHECK("start of step, before alpha/phase-change solve");

                crashStage::set
                (
                    "time " + runTime.timeName()
                  + ": mixture.solve (alphas + phase change)"
                );
                vDot = mixture.solve(&mass_dot, &vDotP);
                vDot.correctBoundaryConditions();
                mass_dot.correctBoundaryConditions();
                vDotP.correctBoundaryConditions();


                crashStage::set
                (
                    "time " + runTime.timeName()
                  + ": update.H (condensate, solidus/liquidus, resistivity)"
                );
                #include "update.H"

                HEALTH_CHECK("after alpha/phase-change solve");
            }

            
            crashStage::set("time " + runTime.timeName() + ": rho = mixture.rho()");
            rho = mixture.rho();

            

            // Ray-trace expensive: evaluate once per time step.
            if (pimple.firstIter())
            {
                crashStage::set
                (
                    "time " + runTime.timeName()
                  + ": laser ray tracing / deposition"
                );
                laser.updateDeposition
                (
                    condensateFiltered, n_filtered, electrical_resistivity
                );
            }

            // PIMPLE Section ((U-p)-H-T coupling)

            crashStage::set("time " + runTime.timeName() + ": UEqn");
            #include "UEqn.H"
            if (mthd.valid())
            {
                crashStage::set("time " + runTime.timeName() + ": MTHD solve");
                mthd->solve(phi, U);
            }
            crashStage::set("time " + runTime.timeName() + ": TEqn");
            #include "TEqn.H"

            HEALTH_CHECK("after TEqn");

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                crashStage::set("time " + runTime.timeName() + ": pEqn");
                #include "pEqn.H"
            }

            HEALTH_CHECK("after pEqn");

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        crashStage::set("time " + runTime.timeName() + ": writing");
        mixture.writeOldTimeValues();
        #include "writeOldTimeStorage.H"

        runTime.write();

        // Write ray paths to VTK files
        if (runTime.outputTime())
        {
            laser.writeRayPathsToVTK();
        }

        runTime.printExecutionTime(Info);
    }

    // Write a VTK series file for easy-opening of the ray files
    laser.writeRayPathVTKSeriesFile();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
