/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.
\*---------------------------------------------------------------------------*/

#include "laserHeatSource.H"
#include "fvc.H"
#include "constants.H"
#include "findLocalCell.H"
#include "SortableList.H"
#include "globalIndex.H"

namespace Foam
{

defineTypeNameAndDebug(laserHeatSource, 0);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void laserHeatSource::createInitialRays
(
    List<compactRay>& rays,
    const fvMesh& mesh,
    const vector& currentLaserPosition,
    const scalar laserRadius,
    const label N_sub_divisions,
    const label nRadial,
    const label nAngular,
    const vector& V_incident,
    const scalar Radius_Flavour,
    const scalar Q_cond,
    const scalar beam_radius
) const
{
    DynamicList<vector> initial_points;
    DynamicList<scalar> point_assoc_power;

    const scalarField& yDimI = yDim_;
    const scalar pi = constant::mathematical::pi;

    // --- NEW: Super-Gaussian n=4 Pre-factor calculation ---
    // This ensures total power is conserved for the n=4 profile
    const scalar preFactor = 
        (2.0 * Foam::sqrt(Radius_Flavour) * Q_cond) 
      / (Foam::sqr(beam_radius) * Foam::pow(pi, 1.5));

    if (radialPolarHeatSource())
    {
        Info<< "nRadial: " << nRadial << nl
            << "nAngular: "<< nAngular <<endl;

        const scalar rMax = 1.5*beam_radius;
        const label totalSamples = nRadial * nAngular;
        const label samplesPerProc = totalSamples/Pstream::nProcs();
        const label remainder = totalSamples % Pstream::nProcs();
        const label myRank = Pstream::myProcNo();
        const label startIdx = myRank * samplesPerProc + min(myRank, remainder);
        const label endIdx = startIdx + samplesPerProc + (myRank < remainder ? 1 : 0);
        const label localSamples = endIdx - startIdx;

        List<scalar> radialPoints(nRadial);
        for (label iR = 0; iR < nRadial; ++iR)
        {
            const scalar fraction = scalar(iR + 0.5)/nRadial;
            radialPoints[iR] = rMax * pow(fraction,1.0);
        }

        const point P0 (currentLaserPosition.x(),currentLaserPosition.y(),currentLaserPosition.z());
        const vector V_i(V_incident/(mag(V_incident) + SMALL));

        const vector a = (mag(V_i.z()) < 0.9) ? vector(0, 0, 1) : vector(0, 1, 0);
        vector u = (V_i ^ a);
        u = u/mag(u);
        const vector v = (V_i ^ u);
        const vector perturbation (1e-10,1e-10,1e-10);

        for (label localIdx = 0; localIdx < localSamples; ++localIdx)
        {
            const label globalIdx = startIdx + localIdx;
            const label iTheta = globalIdx/nRadial;
            const label iR = globalIdx % nRadial;

            const scalar theta = 2.0*pi*iTheta/nAngular;
            const scalar r = radialPoints[iR];
            const scalar deltaTheta = 2.0*pi/nAngular;

            scalar deltaR = (iR == 0) ? radialPoints[0] : (radialPoints[iR] - radialPoints[iR - 1]);
            const scalar area = r*deltaR*deltaTheta;

            const scalar x_local = r*cos(theta);
            const scalar y_local = r*sin(theta);
            const vector globalPos = P0 + x_local*u + y_local*v;

            initial_points.append(globalPos + perturbation);

            // --- AMENDED: Super-Gaussian n=4 exponent ---
            point_assoc_power.append
            (
                area * preFactor * Foam::exp
                (
                  - Radius_Flavour * (Foam::pow(r, 4.0) / Foam::pow(beam_radius, 4.0))
                )
            );
        }
    }
    else 
    {
        const vectorField& CI = mesh.C();

        forAll(CI, celli)
        {
            const scalar x_coord = CI[celli].x();
            const scalar z_coord = CI[celli].z();

            const scalar r = sqrt(sqr(x_coord - currentLaserPosition.x()) + sqr(z_coord - currentLaserPosition.z()));

            if (r <= (1.5*beam_radius) && laserBoundary_[celli] > SMALL)
            {
                for (label Ray_j = 0; Ray_j < N_sub_divisions; Ray_j++)
                {
                    for (label Ray_k = 0; Ray_k < N_sub_divisions; Ray_k++)
                    {
                        const point p_1
                        (
                            CI[celli].x() - (yDimI[celli]/2.0) + ((yDimI[celli]/(N_sub_divisions+1))*(Ray_j+1)),
                            CI[celli].y(),
                            CI[celli].z() - (yDimI[celli]/2.0) + ((yDimI[celli]/(N_sub_divisions+1))*(Ray_k+1))
                        );

                        initial_points.append(p_1);

                        // --- AMENDED: Super-Gaussian n=4 exponent ---
                        point_assoc_power.append
                        (
                            sqr(yDimI[celli]/N_sub_divisions) * preFactor * Foam::exp
                            (
                              - Radius_Flavour * (Foam::pow(r, 4.0) / Foam::pow(beam_radius, 4.0))
                            )
                        );
                    }
                }
            }
        }
    }

    // Gather and broadcast logic remains unchanged
    List<pointField> gatheredData(Pstream::nProcs());
    List<scalarField> gatheredData_powers(Pstream::nProcs());
    gatheredData[Pstream::myProcNo()] = initial_points;
    Pstream::gatherList(gatheredData);
    gatheredData_powers[Pstream::myProcNo()] = point_assoc_power;
    Pstream::gatherList(gatheredData_powers);
    Pstream::broadcastList(gatheredData);
    Pstream::broadcastList(gatheredData_powers);

    pointField rayCoords(ListListOps::combine<Field<vector>>(gatheredData, accessOp<Field<vector>>()));
    scalarField rayPowers(ListListOps::combine<Field<scalar>>(gatheredData_powers, accessOp<Field<scalar>>()));

    rays.setSize(rayCoords.size());
    forAll(rays, i)
    {
        rays[i] = compactRay(rayCoords[i], V_incident, rayPowers[i]);
        rays[i].globalRayIndex_ = i;
        rays[i].currentCell_ = mesh.findCell(rayCoords[i]);
        rays[i].path_.append(rayCoords[i]);
    }
}

// ... (Constructors remain unchanged) ...

void laserHeatSource::updateDeposition
(
    const volScalarField& alphaFiltered,
    const volVectorField& nFiltered,
    const volScalarField& resistivity_in,
    const label laserID,
    const vector& currentLaserPosition,
    const scalar currentLaserPower,
    const scalar laserRadius,
    const label N_sub_divisions,
    const label nRadial,
    const label nAngular,
    const vector& V_incident,
    const scalar wavelength,
    const scalar e_num_density,
    const scalar dep_cutoff,
    const scalar Radius_Flavour,
    const Switch useLocalSearch,
    const label maxLocalSearch,
    const scalar rayPowerRelTol,
    const boundBox& globalBB
)
{
    // ... (Setup logic remains unchanged) ...

    while (myCellID != -1)
    {
        // ... (Ray propagation logic) ...

        if (mag(nFilteredI[myCellID]) > 0.5 && alphaFilteredI[myCellID] >= dep_cutoff)
        {
            // --- AMENDED: Absorptivity set to 1.0 for Energy Conservation Check ---
            // Original Fresnel calculations bypassed per supervisor instructions
            scalar absorptivity = 1.0; 

            deposition_[myCellID] += absorptivity * curRay.power_ / VI[myCellID];
            curRay.power_ *= (1.0 - absorptivity);

            const vector dRef = curRay.direction_ - 2.0*(curRay.direction_ & n)*n;
            curRay.direction_ = dRef;
        }
        else if (alphaFilteredI[myCellID] >= dep_cutoff)
        {
            deposition_[myCellID] += curRay.power_/VI[myCellID];
            curRay.power_ = 0.0;
            curRay.path_.append(curRay.position_);
            break;
        }
        curRay.path_.append(curRay.position_);
    }

    // ... (Sync and Path recording logic) ...

    const scalar TotalQ = fvc::domainIntegrate(deposition_).value();
    Info<< "    Total Q deposited: " << TotalQ << endl;
}

} // End namespace Foam
