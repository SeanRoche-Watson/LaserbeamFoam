/*---------------------------------------------------------------------------*\
License
    This file is part of laserbeamFoam.

    laserbeamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    laserbeamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with laserbeamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "laserHeatSource.H"
#include "fvc.H"
#include "constants.H"
#include "SortableList.H"

#include "findLocalCell.H"
#include "IFstream.H"
#include "OSspecific.H"
#include <fstream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace
{
    // Normalised (super-)Gaussian beam intensity [W/m^2]
    //
    //   I(r) = Q*RF^(1/n)/(pi*w^2*Gamma(1 + 1/n)) * exp(-RF*(r/w)^(2n))
    //
    // RF = Radius_Flavour (2 -> w is the 1/e^2 radius),
    // n  = superGaussianOrder.
    // The prefactor makes the integral over the plane equal Q for any n;
    // n = 1 gives the original RF*Q/(pi*w^2)*exp(-RF*r^2/w^2).
    inline scalar beamIntensity
    (
        const scalar r,
        const scalar w,
        const scalar Q,
        const scalar RF,
        const scalar n
    )
    {
        const scalar pi = constant::mathematical::pi;
        const scalar I0 =
            Q*Foam::pow(RF, 1.0/n)/(pi*sqr(w)*std::tgamma(1.0 + 1.0/n));
        return I0*Foam::exp(-RF*Foam::pow(r/w, 2.0*n));
    }

    // Integral over the plane of exp(-(r - R)^2/(2 s^2)):
    //   2 pi [ s^2 exp(-R^2/(2 s^2)) + R s sqrt(pi/2) (1 + erf(R/(sqrt(2) s))) ]
    inline scalar ringArea(const scalar R, const scalar s)
    {
        const scalar pi = constant::mathematical::pi;
        return
            2.0*pi
           *(
                sqr(s)*Foam::exp(-sqr(R)/(2.0*sqr(s)))
              + R*s*Foam::sqrt(pi/2.0)
               *(1.0 + std::erf(R/(Foam::sqrt(2.0)*s)))
            );
    }

    // Multi-term beam intensity [W/m^2]; each term carries fraction*Q
    inline scalar multiTermIntensity
    (
        const scalar r,
        const scalar Q,
        const List<FixedList<scalar, 4>>& terms
    )
    {
        const scalar pi = constant::mathematical::pi;
        scalar I = 0;
        forAll(terms, i)
        {
            const FixedList<scalar, 4>& t = terms[i];
            if (t[0] < 0.5)
            {
                I += t[1]*Q*2.0/(pi*sqr(t[2]))*Foam::exp(-2.0*sqr(r/t[2]));
            }
            else
            {
                I += t[1]*Q/ringArea(t[2], t[3])
                    *Foam::exp(-sqr(r - t[2])/(2.0*sqr(t[3])));
            }
        }
        return I;
    }

    // Piecewise-linear tabulated profile G(r) (zero beyond the last node)
    inline scalar tableValue
    (
        const scalar r,
        const scalarField& R,
        const scalarField& G
    )
    {
        if (r >= R.last()) return 0;
        label i = 0;
        while (i < R.size() - 2 && r > R[i + 1]) ++i;
        const scalar f = (r - R[i])/(R[i + 1] - R[i] + VSMALL);
        return (1 - f)*G[i] + f*G[i + 1];
    }

    // Plane integral of a piecewise-linear G(r): sum 2 pi int r G dr
    inline scalar tableIntegral(const scalarField& R, const scalarField& G)
    {
        const scalar pi = constant::mathematical::pi;
        scalar s = 0;
        for (label i = 0; i < R.size() - 1; ++i)
        {
            const scalar r0 = R[i], r1 = R[i + 1];
            s += 2*pi*(r1 - r0)
                *(G[i]*(2*r0 + r1) + G[i + 1]*(r0 + 2*r1))/6.0;
        }
        return s;
    }

    // Seeding radius that captures each term
    inline scalar multiTermRMax(const List<FixedList<scalar, 4>>& terms)
    {
        scalar rMax = 0;
        forAll(terms, i)
        {
            const FixedList<scalar, 4>& t = terms[i];
            rMax = max(rMax, t[0] < 0.5 ? 1.5*t[2] : t[2] + 3.0*t[3]);
        }
        return rMax;
    }
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laserHeatSource, 0);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


void laserHeatSource::seedRayCloud
(
    Cloud<laserRayParticle>& cloud,
    const fvMesh& mesh,
    const vector& currentLaserPosition,
    const scalar laserRadius,
    const label N_sub_divisions,
    const label nRadial,
    const label nAngular,
    const vector& V_incident,
    const scalar Radius_Flavour,
    const scalar superGaussianOrder,
    const scalar Q_cond,
    const scalar beam_radius,
    label& nTotalRays
) const
{
    // ================================================================
    // Compute ALL ray starting positions and powers identically on
    // every processor (same as the old createInitialRays), but only
    // add a particle to the cloud on the processor that owns the
    // starting cell.
    // ================================================================

    const scalarField& yDimI = yDim_;
    const scalar pi = constant::mathematical::pi;

    // -- Build global list of ray positions and powers --

    DynamicList<vector> initial_points;
    DynamicList<scalar> point_assoc_power;

    // All seed points will be perturbed by a small amount to avoid landing
    // exactly on a face
    // The magnitude of the perturbation is calculated as a small fraction of
    // the smallest cell dimension, as opposed to an absolute value, and we take
    // care to perturb in the ray direction and transverse to the ray direction

    // Ray direction
    const vector d = V_incident/(mag(V_incident) + VSMALL);

    // Pick a stable “arbitrary” axis not parallel to d
    const vector a =
        (mag(d & vector(1,0,0)) < 0.9) ? vector(1,0,0) : vector(0,0,1);

    // Perpendicular unit vector
    vector t = d ^ a;
    t /= (mag(t) + VSMALL);

    // Minimum length scale
    const scalar lCell = cbrt(gMin(mesh.V()));
    const scalar eps = sqrt(SMALL)*lCell;

    // Safe perturbation in the ray direction and transverse ray direction
    const vector perturbation = eps*(d + t);
    Info<< "    Perturbation vector for ray seed points = " << perturbation
        << endl;

    if (radialPolarHeatSource())
    {
        Info<< "superGaussianOrder: " << superGaussianOrder << nl
            << "nRadial: " << nRadial << nl
            << "nAngular: " << nAngular << endl;

        const scalar rMax =
            ctrlR_.size() ? ctrlR_.last()
          : beamTerms_.size() ? multiTermRMax(beamTerms_)
          : 1.5*beam_radius;
        Info<< "beam seeding radius rMax: " << rMax << endl;

        // ================================================================
        // Improved radial discretisation for smoother Gaussian
        //
        // Use cell-edge based radial boundaries to compute exact annulus
        // areas. Each ring iR covers r from r_inner[iR] to r_outer[iR].
        // The sample point is at the area-weighted centroid of the annulus.
        // ================================================================

        // Compute radial boundaries (edges of annular rings)
        List<scalar> radialBoundaries(nRadial + 1);
        for (label iR = 0; iR <= nRadial; ++iR)
        {
            // Uniform spacing of boundaries
            radialBoundaries[iR] = rMax * scalar(iR) / scalar(nRadial);
        }

        // Compute sample points (area-weighted centroid of each annulus)
        // For an annulus from r1 to r2, the centroid is at:
        //   r_centroid = (2/3) * (r2^3 - r1^3) / (r2^2 - r1^2)
        List<scalar> radialPoints(nRadial);
        List<scalar> annulusAreas(nRadial);

        for (label iR = 0; iR < nRadial; ++iR)
        {
            const scalar r_inner = radialBoundaries[iR];
            const scalar r_outer = radialBoundaries[iR + 1];

            // Exact annulus area (full ring, will be divided by nAngular)
            annulusAreas[iR] = pi * (sqr(r_outer) - sqr(r_inner));

            // Area-weighted centroid for sample position
            if (r_inner < SMALL)
            {
                // Central disc: centroid at 2/3 * r_outer
                radialPoints[iR] = (2.0/3.0) * r_outer;
            }
            else
            {
                // Annulus: use exact centroid formula
                radialPoints[iR] =
                    (2.0/3.0)
                   * (pow3(r_outer) - pow3(r_inner))
                   / (sqr(r_outer) - sqr(r_inner));
            }
        }

        const label totalSamples = nRadial * nAngular;
        const label samplesPerProc = totalSamples/Pstream::nProcs();
        const label remainder = totalSamples % Pstream::nProcs();
        const label myRank = Pstream::myProcNo();
        const label startIdx =
            myRank * samplesPerProc + min(myRank, remainder);
        const label endIdx =
            startIdx + samplesPerProc + (myRank < remainder ? 1 : 0);
        const label localSamples = endIdx - startIdx;

        const point P0
        (
            currentLaserPosition.x(),
            currentLaserPosition.y(),
            currentLaserPosition.z()
        );

        const vector V_i(V_incident/(mag(V_incident) + SMALL));

        const vector a =
            (mag(V_i.z()) < 0.9) ? vector(0, 0, 1) : vector(0, 1, 0);
        vector u = (V_i ^ a);
        u = u/mag(u);
        const vector v = (V_i ^ u);

        for (label localIdx = 0; localIdx < localSamples; ++localIdx)
        {
            const label globalIdx = startIdx + localIdx;
            const label iTheta = globalIdx/nRadial;
            const label iR = globalIdx % nRadial;

            const scalar theta = 2.0*pi*iTheta/nAngular;
            const scalar r = radialPoints[iR];

            // Area element for this ray (annulus sector)
            const scalar area = annulusAreas[iR] / scalar(nAngular);

            const scalar x_local = r*cos(theta);
            const scalar y_local = r*sin(theta);

            const vector globalPos = P0 + x_local*u + y_local*v;

            initial_points.append(globalPos + perturbation);

            // Gaussian power distribution using the sample point radius
            point_assoc_power.append
            (
                area
               *(
                    ctrlR_.size()
                  ? Q_cond*tableValue(r, ctrlR_, ctrlG_)
                  : beamTerms_.size()
                  ? multiTermIntensity(r, Q_cond, beamTerms_)
                  : beamIntensity
                    (
                        r, beam_radius, Q_cond, Radius_Flavour,
                        superGaussianOrder
                    )
                )
            );
        }
    }
    else // One ray for each boundary patch face within the laser radius
    {
        const vectorField& CI = mesh.C();

        forAll(CI, celli)
        {
            const scalar x_coord = CI[celli].x();
            const scalar z_coord = CI[celli].z();

            const scalar r =
                sqrt
                (
                    sqr(x_coord - currentLaserPosition.x())
                  + sqr(z_coord - currentLaserPosition.z())
                );

            if
            (
                r <= (ctrlR_.size() ? ctrlR_.last()
                      : beamTerms_.size()
                      ? multiTermRMax(beamTerms_) : 1.5*beam_radius)
             && laserBoundary_[celli] > SMALL
            )
            {
                for (label Ray_j = 0; Ray_j < N_sub_divisions; Ray_j++)
                {
                    for (label Ray_k = 0; Ray_k < N_sub_divisions; Ray_k++)
                    {
                        point p_1
                        (
                            CI[celli].x()
                          - (yDimI[celli]/2.0)
                          + ((yDimI[celli]/(N_sub_divisions+1))*(Ray_j+1)),
                            CI[celli].y(),
                            CI[celli].z()
                          - (yDimI[celli]/2.0)
                          + ((yDimI[celli]/(N_sub_divisions+1))*(Ray_k+1))
                        );

                        // Perturb the point
                        p_1 += perturbation;

                        initial_points.append(p_1);

                        point_assoc_power.append
                        (
                            sqr(yDimI[celli]/N_sub_divisions)
                           *(
                                ctrlR_.size()
                              ? Q_cond*tableValue(r, ctrlR_, ctrlG_)
                              : beamTerms_.size()
                              ? multiTermIntensity(r, Q_cond, beamTerms_)
                              : beamIntensity
                                (
                                    r,
                                    beam_radius,
                                    Q_cond,
                                    Radius_Flavour,
                                    superGaussianOrder
                                )
                            )
                        );
                    }
                }
            }
        }
    }

    // List with size equal to number of processors
    // Gather all initial points/powers to all processors
    List<pointField> gatheredData(Pstream::nProcs());
    List<scalarField> gatheredPowers(Pstream::nProcs());

    // Populate and gather the list onto the master processor.
    // Distibulte the data accross the different processors
    gatheredData[Pstream::myProcNo()] = initial_points;
    Pstream::gatherList(gatheredData);
    Pstream::broadcastList(gatheredData);

    gatheredPowers[Pstream::myProcNo()] = point_assoc_power;
    Pstream::gatherList(gatheredPowers);
    Pstream::broadcastList(gatheredPowers);

    // List of initial points
    pointField rayCoords
    (
        ListListOps::combine<Field<vector>>
        (
            gatheredData,
            accessOp<Field<vector>>()
        )
    );

    scalarField rayPowers
    (
        ListListOps::combine<Field<scalar>>
        (
            gatheredPowers,
            accessOp<Field<scalar>>()
        )
    );


    // Create a list of Ray objects
    nTotalRays = rayCoords.size();

    label seedCellI = -1;

    forAll(rayCoords, i)
    {
        const label cellI = findLocalCell
            (
                rayCoords[i],
                seedCellI,
                mesh,
                100,    // maxLocalSearch
                false   // debug
            );

        if (cellI >= 0)
        {
            // Update seed for next iteration - nearby rays benefit
            seedCellI = cellI;

            // This processor owns this ray's starting cell
            laserRayParticle* pPtr = new laserRayParticle
                (
                    mesh,
                    rayCoords[i],
                    cellI,
                    V_incident,
                    rayPowers[i],
                    0.0,        // dA
                    i           // globalRayIndex
                );

            cloud.addParticle(pPtr);
        }
    }

    Info<< "    Total rays: " << nTotalRays
        << ", local particles: " << cloud.size() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

laserHeatSource::laserHeatSource
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "LaserProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    deposition_
    (
        IOobject
        (
            "Deposition",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("deposition", dimensionSet(1, -1, -3, -0, 0), -1.0)
    ),
    laserBoundary_
    (
        IOobject
        (
            "Laser_boundary",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    errorTrack_
    (
        IOobject
        (
            "errorTrack",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "errorTrack", dimensionSet(0, 0, 0, -0, 0), 0.0
        )
    ),
    rayNumber_
    (
        IOobject
        (
            "rayNumber",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rayNumber", dimensionSet(0, 0, 0, -0, 0), -1.0)
    ),
    rayQ_
    (
        IOobject
        (
            "rayQ",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rayQ", dimensionSet(1, 0, -3, 0, 0), scalar(0.0))
    ),
    yDim_
    (
        IOobject
        (
            "yDim",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("yDim", dimensionSet(0, 1, 0, 0, 0), 1.0)
    ),
    powderSim_(lookupOrDefault<Switch>("PowderSim", false)), // Needed for cases where we have spherical particles (that have jagged edges or cartesian grids) which will cause spurious currents due to these large gradients - we extend the damping term by 1 cell
    radialPolarHeatSource_
    (
        found("radialPolarHeatSource")
      ? Switch(lookup("radialPolarHeatSource"))
      : lookupOrDefault<Switch>("Radial_Polar_HS", true)
    ),
    laserNames_(0),
    laserDicts_(0),
    timeVsLaserPosition_(0),
    timeVsLaserPower_(0),
    rayPaths_(0),
    vtkTimes_(),
    globalBB_(mesh.bounds()),
    beamTerms_(),
    ctrlR_(),
    ctrlG_(),
    ctrlInit_(false),
    ctrlLastUpdate_(-GREAT)
{
    Info<< "radialPolarHeatSource = " << radialPolarHeatSource_ << endl;

    // Calculate global bounding box
    {
        boundBox localBB = mesh.bounds();
        globalBB_ = localBB;

        reduce(globalBB_.min(), minOp<vector>());
        reduce(globalBB_.max(), maxOp<vector>());

        globalBB_.inflate(0.01);
        // Inflate the bounding box by 1% to avoid issues with rays starting just at the mesh boundary -> probably not an issue but worth doing

        Info<< "Scaled global mesh bounding box: " << globalBB_ << endl;
    }

    // Initialise the laser power and position
    if (found("lasers"))
    {
        const PtrList<entry> laserEntries(lookup("lasers"));

        laserNames_.setSize(laserEntries.size());
        laserDicts_.setSize(laserEntries.size());
        timeVsLaserPosition_.setSize(laserEntries.size());
        timeVsLaserPower_.setSize(laserEntries.size());

        if (Pstream::master())
        {
            rayPaths_.setSize(laserEntries.size());
        }

        forAll(laserEntries, laserI)
        {
            laserNames_[laserI] = laserEntries[laserI].keyword();
            Info<< "Reading laser " << laserNames_[laserI] << endl;

            laserDicts_.set
            (
                laserI, new dictionary(laserEntries[laserI].dict())
            );

            timeVsLaserPosition_.set
            (
                laserI,
                new interpolationTable<vector>
                (
                    laserEntries[laserI].dict().subDict("timeVsLaserPosition")
                )
            );

            timeVsLaserPower_.set
            (
                laserI,
                new interpolationTable<scalar>
                (
                    laserEntries[laserI].dict().subDict("timeVsLaserPower")
                )
            );
        }

            // Check that a single laser is not also defined

        if (found("timeVsLaserPosition"))
        {
            FatalErrorInFunction
                << "timeVsLaserPosition should not be defined in the main dict"
                << " if a list of lasers is provided" << exit(FatalError);
        }

        if (found("timeVsLaserPower"))
        {
            FatalErrorInFunction
                << "timeVsLaserPower should not be defined in the main dict"
                << " if a list of lasers is provided" << exit(FatalError);
        }
    }
    else
    {
        // There is no lists of lasers, just one
        // Backward compatibility: single laser specified in main dict
        rayPaths_.setSize(1);
        laserNames_.setSize(1);
        laserDicts_.setSize(1);
        timeVsLaserPosition_.setSize(1);
        timeVsLaserPower_.setSize(1);

        laserNames_[0] = "laser0";

        // Copy the main dict
        laserDicts_.set(0, new dictionary(*this));

        timeVsLaserPosition_.set
        (
            0,
            new interpolationTable<vector>(subDict("timeVsLaserPosition"))
        );

        timeVsLaserPower_.set
        (
            0,
            new interpolationTable<scalar>(subDict("timeVsLaserPower"))
        );
    }

        // Update laserBoundary - only used in old boundary based initialisation
    laserBoundary_ = fvc::average(laserBoundary_);

    if (debug)
    {
        errorTrack_.writeOpt() = IOobject::AUTO_WRITE;
        rayNumber_.writeOpt() = IOobject::AUTO_WRITE;
    }

    // Give errors if the old input format is found
    if (found("HS_bg"))
    {
        FatalErrorInFunction
            << "'HS_bg' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }
    // Give errors if the old input format is found
    if (found("HS_lg"))
    {
        FatalErrorInFunction
            << "'HS_lg' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }
    // Give errors if the old input format is found
    if (found("HS_velocity"))
    {
        FatalErrorInFunction
            << "'HS_velocity' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }
    // Give errors if the old input format is found
    if (found("HS_Q"))
    {
        FatalErrorInFunction
            << "'HS_Q' is deprecated: please instead specify the laser "
            << "power in time via the laserPowereVsTime sub-dict"
            << exit(FatalError);
    }
    // Give errors if the old input format is found
    if (found("elec_resistivity"))
    {
        FatalErrorInFunction
            << "'elec_resistivity' is deprecated: resistivity is now "
            << "passed in from the solver as a field"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void laserHeatSource::updateDeposition
(
    const volScalarField& alphaFiltered,
    const volVectorField& nFiltered,
    const volScalarField& resistivity_in
)
{
    Info<< "Updating deposition" << endl;

    deposition_ *= 0.0;
    laserBoundary_ *= 0.0;
    laserBoundary_ = fvc::average(laserBoundary_);
    errorTrack_ *= 0.0;
    rayNumber_ *= 0.0;
    rayQ_ *= 0.0;

    const scalar time = deposition_.time().value();

    forAll(laserNames_, laserI)
    {
        // Lookup the current laser position and power
        vector currentLaserPosition =
            timeVsLaserPosition_[laserI](time);
        const scalar currentLaserPower =
            timeVsLaserPower_[laserI](time);

        Info<< "Laser: " << laserNames_[laserI] << nl
            << "    mean position = " << currentLaserPosition << nl
            << "    power = " << currentLaserPower << endl;

            // Dict for current laser
        const dictionary& dict = laserDicts_[laserI];

        // If defined, add oscillation to laser position
        if (dict.found("HS_oscAmpX"))
        {
            const scalar oscAmpX(readScalar(dict.lookup("HS_oscAmpX")));
            const scalar oscFreqX(readScalar(dict.lookup("HS_oscFreqX")));
            const scalar pi = constant::mathematical::pi;
            currentLaserPosition[vector::X] +=
                oscAmpX*sin(2*pi*oscFreqX*time);
        }// If defined, add oscillation to laser position

        if (dict.found("HS_oscAmpZ"))
        {
            const scalar oscAmpZ(readScalar(dict.lookup("HS_oscAmpZ")));
            const scalar oscFreqZ(readScalar(dict.lookup("HS_oscFreqZ")));
            const scalar pi = constant::mathematical::pi;
            currentLaserPosition[vector::Z] +=
                oscAmpZ*cos(2*pi*oscFreqZ*time);
        }// If defined, add oscillation to laser position

        Info<< "    position including any oscillation = "
            << currentLaserPosition << endl;

        scalar laserRadius = 0.0;
        if (dict.found("HS_a") && dict.found("laserRadius"))
        {
            FatalErrorInFunction
                << "The laser radius should be specified via 'laserRadius' or "
                << "'HS_a', not both!" << exit(FatalError);
        }

        if (dict.found("HS_a"))
        {
            laserRadius = readScalar(dict.lookup("HS_a"));
        }
        else if (dict.found("laserRadius"))
        {
            laserRadius = readScalar(dict.lookup("laserRadius"));
        }
        else
        {
            FatalErrorInFunction
                << "The laser radius should be specified via 'laserRadius' "
                << "or 'HS_a'"
                << exit(FatalError);
        }

        const label nRadial(dict.lookupOrDefault<label>("nRadial", 5));
        const label nAngular(dict.lookupOrDefault<label>("nAngular", 30));
        const label N_sub_divisions
        (
            dict.lookupOrDefault<label>("N_sub_divisions", 1)
        );
        const vector V_incident(dict.lookup("V_incident"));
        const scalar wavelength(readScalar(dict.lookup("wavelength")));
        const scalar e_num_density
        (
            readScalar(dict.lookup("e_num_density"))
        );
        const scalar dep_cutoff
        (
            dict.lookupOrDefault<scalar>("dep_cutoff", 0.5)
        );
        const scalar Radius_Flavour
        (
            dict.lookupOrDefault<scalar>("Radius_Flavour", 2.0)
        );
        // Super-Gaussian order n: I ~ exp(-Radius_Flavour*(r/w)^(2n)).
        // n = 1 (default) is the original Gaussian; n = 2 is a flat-topped
        // multimode profile. Power is renormalised so the beam still
        // carries the full laser power for any n.
        const scalar superGaussianOrder
        (
            dict.lookupOrDefault<scalar>("superGaussianOrder", 1.0)
        );
        if (superGaussianOrder < 0.5)
        {
            FatalErrorInFunction
                << "superGaussianOrder must be >= 0.5, got "
                << superGaussianOrder << exit(FatalError);
        }
        // Optional multi-term beam profile, e.g.
        //   beamTerms
        //   {
        //       core { type gaussian; fraction 0.67; w 850e-6; }
        //       halo { type ring; fraction 0.33; R 1109e-6; s 334e-6; }
        //   }
        // gaussian: I ~ exp(-2 r^2/w^2) (w = 1/e^2 radius)
        // ring:     I ~ exp(-(r - R)^2/(2 s^2))
        // Each term is normalised to carry fraction*laser power (fractions
        // are rescaled to sum to 1). When present, laserRadius,
        // Radius_Flavour and superGaussianOrder do not set the profile.
        beamTerms_.clear();
        if (dict.found("beamTerms"))
        {
            const dictionary& bt = dict.subDict("beamTerms");
            const wordList names(bt.toc());
            DynamicList<FixedList<scalar, 4>> terms;
            scalar fSum = 0;
            forAll(names, i)
            {
                const dictionary& td = bt.subDict(names[i]);
                const word type(td.lookup("type"));
                FixedList<scalar, 4> t(0.0);
                t[1] = readScalar(td.lookup("fraction"));
                if (type == "gaussian")
                {
                    t[0] = 0; t[2] = readScalar(td.lookup("w"));
                }
                else if (type == "ring")
                {
                    t[0] = 1;
                    t[2] = readScalar(td.lookup("R"));
                    t[3] = readScalar(td.lookup("s"));
                }
                else
                {
                    FatalErrorInFunction
                        << "beamTerms/" << names[i] << ": unknown type "
                        << type << " (use gaussian or ring)"
                        << exit(FatalError);
                }
                if (t[1] < 0 || t[2] < 0 || (t[0] > 0.5 && t[3] <= 0)
                 || (t[0] < 0.5 && t[2] <= 0))
                {
                    FatalErrorInFunction
                        << "beamTerms/" << names[i]
                        << ": fraction/w/R must be >= 0 and w, s > 0"
                        << exit(FatalError);
                }
                fSum += t[1];
                terms.append(t);
            }
            if (fSum <= 0)
            {
                FatalErrorInFunction
                    << "beamTerms: fractions sum to zero" << exit(FatalError);
            }
            forAll(terms, i)
            {
                terms[i][1] /= fSum;
            }
            beamTerms_.transfer(terms);
            Info<< "    multi-term beam (type 0 gaussian w | 1 ring R s): "
                << beamTerms_ << " (fractions summed to " << fSum << ")"
                << endl;
        }

        if (dict.found("beamController"))
        {
            updateBeamController
            (
                dict.subDict("beamController"),
                currentLaserPosition,
                alphaFiltered,
                laserRadius,
                Radius_Flavour,
                superGaussianOrder
            );
        }

        const Switch useLocalSearch
        (
            dict.lookupOrDefault<Switch>("useLocalSearch", true)
        );
        const label maxLocalSearch
        (
            dict.lookupOrDefault<label>("maxLocalSearch", 100)
        );
        const scalar rayPowerRelTol
        (
            dict.lookupOrDefault<scalar>("rayPowerRelTol", 1e-4)
        );
        const label maxRayBounces
        (
            dict.lookupOrDefault<label>("maxRayBounces", 1000)
        );
        const word maxRayBounceAction
        (
            dict.lookupOrDefault<word>("maxRayBounceAction", "deposit")
        );
        const scalar nearZeroAbsorptivityTol
        (
            dict.lookupOrDefault<scalar>("nearZeroAbsorptivityTol", SMALL)
        );

        updateDeposition
        (
            alphaFiltered,
            nFiltered,
            resistivity_in,
            laserI,
            currentLaserPosition,
            currentLaserPower,
            laserRadius,
            N_sub_divisions,
            nRadial,
            nAngular,
            V_incident,
            wavelength,
            e_num_density,
            dep_cutoff,
            Radius_Flavour,
            superGaussianOrder,
            useLocalSearch,
            maxLocalSearch,
            rayPowerRelTol,
            maxRayBounces,
            maxRayBounceAction,
            nearZeroAbsorptivityTol,
            globalBB_
        );
    }

    deposition_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * Beam controller * * * * * * * * * * * * * * //
//
//  Optional sub-dictionary in LaserProperties:
//
//  beamController
//  {
//      active        true;     // false: fixed tabulated profile only
//      profileFile   "...";    // optional: start from (or, if inactive,
//                              // use) a saved profile; otherwise the
//                              // controller starts from the case's own
//                              // beam (beamTerms or laserRadius)
//      resume        true;     // re-read <case>/beamController/latest on
//                              // a restart (default true)
//      nNodes        21;  rMax 2e-3;     // radial nodes 0..rMax
//      Tiso          1660;     // fusion isotherm [K]
//      alphaMin      0.5;      // metal cells: alpha_filtered >= alphaMin
//      surfaceLevel  0;        // original plate top (depthDir coordinate)
//      depthDir      (0 1 0);  // into the plate
//      crossDir      (1 0 0);  // across the track
//      binWidth      50e-6;
//      start 0.1;  end 0.37;  interval 0.005;    // [s]
//      gain          0.05;     // per update, on log G
//      maxStep       0.1;      // cap on |d log G| per update
//      smooth        0.25;     // neighbour smoothing of d log G
//      conservePower false;    // true: shape only, total power fixed
//      maxPowerScale 3;        // cap on total power / set power
//      depthOffset   50e-6;    // regularises the log error near edges
//      targetDepth   ( (x0 D0) (x1 D1) ... );   // |x| [m], depth [m]
//  }
//
//  The measured quantity is the current melt pool projected onto the
//  cross-section: for every |x| bin, the deepest metal cell with
//  T >= Tiso (anywhere along the track). In steady state this is the
//  fusion line the pool leaves behind. Each update compares it with
//  targetDepth and moves the radial node at r ~ |x| by
//      d log G = gain * log((D_target + d0)/(D_sim + d0))
//  where a bin with no melt counts as D = -d0/2. Absorption (Fresnel,
//  notch trapping) is still computed by the ray tracer: the controller
//  only changes the incident beam.
//  Writes <case>/beamController/{history.dat, depth.dat, latest}.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void laserHeatSource::updateBeamController
(
    const dictionary& cd,
    const vector& laserPosition,
    const volScalarField& alphaFiltered,
    const scalar laserRadius,
    const scalar Radius_Flavour,
    const scalar superGaussianOrder
)
{
    const fvMesh& mesh = deposition_.mesh();
    const Time& runTime = mesh.time();
    const scalar t = runTime.value();
    const fileName dir(runTime.globalPath()/"beamController");

    const Switch active(cd.lookupOrDefault<Switch>("active", true));

    // ---- initialise the table once ----
    if (!ctrlInit_)
    {
        ctrlInit_ = true;
        fileName src;
        if (cd.lookupOrDefault<Switch>("resume", true)
         && isFile(dir/"latest") && runTime.startTime().value() > SMALL)
        {
            src = dir/"latest";
        }
        else if (cd.found("profileFile"))
        {
            src = cd.get<fileName>("profileFile");
            src.expand();
        }

        if (!src.empty())
        {
            IFstream is(src);
            if (!is.good())
            {
                FatalErrorInFunction
                    << "beamController: cannot read " << src
                    << exit(FatalError);
            }
            const dictionary pd(is);
            ctrlR_ = pd.get<scalarField>("r");
            ctrlG_ = pd.get<scalarField>("G");
            if (ctrlR_.size() < 2 || ctrlR_.size() != ctrlG_.size())
            {
                FatalErrorInFunction
                    << "beamController: r and G in " << src
                    << " must be the same length (>= 2)" << exit(FatalError);
            }
            Info<< "    beamController: profile read from " << src << endl;
        }
        else
        {
            const label nN(cd.lookupOrDefault<label>("nNodes", 21));
            const scalar rM(cd.lookupOrDefault<scalar>("rMax", 2e-3));
            ctrlR_.setSize(nN);
            ctrlG_.setSize(nN);
            forAll(ctrlR_, i)
            {
                const scalar r = rM*scalar(i)/scalar(nN - 1);
                ctrlR_[i] = r;
                ctrlG_[i] =
                    beamTerms_.size()
                  ? multiTermIntensity(r, 1.0, beamTerms_)
                  : beamIntensity
                    (
                        r, laserRadius, 1.0, Radius_Flavour,
                        superGaussianOrder
                    );
            }
            // the table carries exactly the laser power at start
            ctrlG_ /= tableIntegral(ctrlR_, ctrlG_);
            ctrlG_[nN - 1] = 0;
            ctrlG_ /= tableIntegral(ctrlR_, ctrlG_);
            Info<< "    beamController: profile sampled from the case beam"
                << endl;
        }
        if (Pstream::master()) mkDir(dir);
    }

    if (!active) return;

    const scalar tStart(cd.lookupOrDefault<scalar>("start", 0.1));
    const scalar tEnd(cd.lookupOrDefault<scalar>("end", GREAT));
    const scalar interval(cd.lookupOrDefault<scalar>("interval", 0.005));
    if (t < tStart || t > tEnd || t < ctrlLastUpdate_ + interval - SMALL)
    {
        return;
    }
    ctrlLastUpdate_ = t;

    const scalar Tiso(cd.lookupOrDefault<scalar>("Tiso", 1660.0));
    const scalar alphaMin(cd.lookupOrDefault<scalar>("alphaMin", 0.5));
    const scalar y0(cd.lookupOrDefault<scalar>("surfaceLevel", 0.0));
    vector dDir(cd.lookupOrDefault<vector>("depthDir", vector(0, 1, 0)));
    vector xDir(cd.lookupOrDefault<vector>("crossDir", vector(1, 0, 0)));
    dDir /= mag(dDir);
    xDir /= mag(xDir);
    const scalar dx(cd.lookupOrDefault<scalar>("binWidth", 50e-6));
    const scalar gain(cd.lookupOrDefault<scalar>("gain", 0.05));
    const scalar maxStep(cd.lookupOrDefault<scalar>("maxStep", 0.1));
    const scalar smooth(cd.lookupOrDefault<scalar>("smooth", 0.25));
    const scalar d0(cd.lookupOrDefault<scalar>("depthOffset", 50e-6));
    const Switch conservePower
    (
        cd.lookupOrDefault<Switch>("conservePower", false)
    );
    const List<Tuple2<scalar, scalar>> target
    (
        cd.get<List<Tuple2<scalar, scalar>>>("targetDepth")
    );
    if (target.size() < 2)
    {
        FatalErrorInFunction
            << "beamController: targetDepth needs >= 2 points"
            << exit(FatalError);
    }

    const scalar rM = ctrlR_.last();
    const label nb = label(rM/dx) + 1;

    // ---- measure: projected melt-pool depth per |x| bin ----
    const volScalarField& T = mesh.lookupObject<volScalarField>("T");
    const vectorField& Cc = mesh.C();
    scalarField Dsim(nb, -GREAT);
    forAll(Cc, celli)
    {
        if (T[celli] >= Tiso && alphaFiltered[celli] >= alphaMin)
        {
            const scalar xc = mag((Cc[celli] - laserPosition) & xDir);
            const label b = label(xc/dx);
            if (b < nb)
            {
                Dsim[b] = max(Dsim[b], (Cc[celli] & dDir) - y0);
            }
        }
    }
    forAll(Dsim, b)
    {
        reduce(Dsim[b], maxOp<scalar>());
    }

    // target on the same bins (linear interpolation, 0 = no melt)
    scalarField Dexp(nb, -GREAT);
    const scalar xTmax = target.last().first();
    for (label b = 0; b < nb; ++b)
    {
        const scalar xb = (b + 0.5)*dx;
        if (xb > xTmax) continue;
        label i = 0;
        while (i < target.size() - 2 && xb > target[i + 1].first()) ++i;
        const scalar x0 = target[i].first(), x1 = target[i + 1].first();
        const scalar f = min(max((xb - x0)/(x1 - x0 + VSMALL), 0.0), 1.0);
        Dexp[b] = (1 - f)*target[i].second() + f*target[i + 1].second();
    }

    // ---- error per bin (log ratio), skip bins molten in neither ----
    scalarField err(nb, 0);
    boolList used(nb, false);
    scalar wSim = 0, wExp = 0, rms = 0;
    label nUsed = 0;
    for (label b = 0; b < nb; ++b)
    {
        const bool ms = Dsim[b] > -GREAT/2, me = Dexp[b] > -GREAT/2;
        if (ms) wSim = (b + 1)*dx;
        if (me) wExp = (b + 1)*dx;
        if (!ms && !me) continue;
        const scalar ds = ms ? max(Dsim[b], 0.0) : -0.5*d0;
        const scalar de = me ? max(Dexp[b], 0.0) : -0.5*d0;
        err[b] = Foam::log((de + d0)/(ds + d0));
        used[b] = true;
        rms += sqr(err[b]);
        ++nUsed;
    }
    if (!nUsed || wSim <= 0)
    {
        Info<< "    beamController: no melt pool yet, no update" << endl;
        return;
    }

    // ---- map bin errors to radial nodes ----
    const label nN = ctrlR_.size();
    const scalar ell = max(rM/(nN - 1), dx);
    scalarField dlog(nN, 0);
    forAll(ctrlR_, k)
    {
        scalar sw = 0, se = 0;
        for (label b = 0; b < nb; ++b)
        {
            if (!used[b]) continue;
            const scalar w =
                Foam::exp(-sqr(((b + 0.5)*dx - ctrlR_[k])/ell));
            sw += w;
            se += w*err[b];
        }
        // normalise only where the node is well covered by bins, so nodes
        // far outside both pools are not driven by the edge error
        dlog[k] = gain*se/max(sw, 1.0);
    }
    // neighbour smoothing, then cap
    scalarField sm(dlog);
    for (label k = 0; k < nN; ++k)
    {
        const scalar l = dlog[max(k - 1, 0)], r = dlog[min(k + 1, nN - 1)];
        sm[k] = (1 - smooth)*dlog[k] + 0.5*smooth*(l + r);
    }
    const scalar Gmax = max(ctrlG_);
    forAll(ctrlG_, k)
    {
        if (k == nN - 1) continue;    // last node stays 0
        ctrlG_[k] *= Foam::exp(min(max(sm[k], -maxStep), maxStep));
        ctrlG_[k] = max(ctrlG_[k], 1e-6*Gmax);
    }
    const scalar maxScale
    (
        cd.lookupOrDefault<scalar>("maxPowerScale", 3.0)
    );
    if (conservePower)
    {
        ctrlG_ /= tableIntegral(ctrlR_, ctrlG_);
    }
    else if (tableIntegral(ctrlR_, ctrlG_) > maxScale)
    {
        ctrlG_ *= maxScale/tableIntegral(ctrlR_, ctrlG_);
    }
    const scalar scale = tableIntegral(ctrlR_, ctrlG_);

    Info<< "    beamController t = " << t
        << ": power scale = " << scale
        << ", centre depth sim/target = " << max(Dsim[0], 0.0)
        << " / " << max(Dexp[0], 0.0)
        << ", width sim/target = " << 2*wSim << " / " << 2*wExp
        << ", rms log error = " << Foam::sqrt(rms/nUsed) << endl;

    // ---- write ----
    if (Pstream::master())
    {
        const bool first = !isFile(dir/"history.dat");
        std::ofstream hs((dir/"history.dat").c_str(), std::ios::app);
        std::ofstream ds((dir/"depth.dat").c_str(), std::ios::app);
        hs.precision(8);
        ds.precision(6);
        if (first)
        {
            hs << "# t power_scale width_sim width_exp rms | r nodes [m]:";
            forAll(ctrlR_, k) hs << ' ' << ctrlR_[k];
            hs << "\n# then G at each node [1/m^2] (I = P*G)\n";
            ds << "# t | D_sim per bin (bin centre = (i+0.5)*" << dx
               << " m, -1 = no melt) | then D_target per bin\n";
        }
        hs << t << ' ' << scale << ' ' << 2*wSim << ' ' << 2*wExp << ' '
           << Foam::sqrt(rms/nUsed);
        forAll(ctrlG_, k) hs << ' ' << ctrlG_[k];
        hs << '\n';
        ds << t;
        forAll(Dsim, b) ds << ' ' << (Dsim[b] > -GREAT/2 ? Dsim[b] : -1.0);
        forAll(Dexp, b) ds << ' ' << (Dexp[b] > -GREAT/2 ? Dexp[b] : -1.0);
        ds << '\n';

        OFstream ls(dir/"latest");
        ls << "// beamController profile at t = " << t
           << " (I = P*G, power scale " << scale << ")" << nl
           << "r " << ctrlR_ << token::END_STATEMENT << nl
           << "G " << ctrlG_ << token::END_STATEMENT << nl;
        OFstream ts(dir/("profile_" + runTime.timeName()));
        ts << "r " << ctrlR_ << token::END_STATEMENT << nl
           << "G " << ctrlG_ << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Lagrangian particle-based parallel ray-tracing: updateDeposition
//
//  Algorithm:
//    1. Create a Cloud<laserRayParticle>
//    2. Seed it with ray particles (one per ray, owned by local proc)
//    3. Call cloud.move() - OpenFOAM handles:
//       - Exact face-to-face tracking (no step size, no missed cells)
//       - Automatic transfer at processor boundaries
//       - Hit callbacks for walls, patches, processor patches
//    4. Collect ray paths from finished particles for VTK output
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


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
    const scalar superGaussianOrder,
    const Switch useLocalSearch,
    const label maxLocalSearch,
    const scalar rayPowerRelTol,
    const label maxRayBounces,
    const word& maxRayBounceAction,
    const scalar nearZeroAbsorptivityTol,
    const boundBox& globalBB
)
{
    const fvMesh& mesh = deposition_.mesh();
    const scalar pi = constant::mathematical::pi;
    const scalar beam_radius = laserRadius;

    // Material constants
    const scalar plasma_frequency = Foam::sqrt
    (
        (
            e_num_density
           *constant::electromagnetic::e.value()
           *constant::electromagnetic::e.value()
        )
       /(
           constant::atomic::me.value()
          *constant::electromagnetic::epsilon0.value()
       )
    );
    const scalar angular_frequency =
        2.0*pi*constant::universal::c.value()/wavelength;

    // ==================================================================
    // Step 1: Create the ray particle cloud
    // ==================================================================

    Cloud<laserRayParticle> rayCloud
    (
        mesh,
        "laserRays",
        IDLList<laserRayParticle>()
    );

    // ==================================================================
    // Step 2: Seed the cloud with ray particles
    //
    // seedRayCloud computes ALL ray positions/powers globally, but
    // each processor only adds particles for cells it owns.
    // ==================================================================

    label nTotalRays = 0;

    seedRayCloud
    (
        rayCloud,
        mesh,
        currentLaserPosition,
        laserRadius,
        N_sub_divisions,
        nRadial,
        nAngular,
        V_incident,
        Radius_Flavour,
        superGaussianOrder,
        currentLaserPower,
        beam_radius,
        nTotalRays
    );

    // ==================================================================
    // Step 3: Move the particles through the mesh
    //
    //  Cloud::move() handles everything:
    //    - Face-to-face tracking via trackToAndHitFace()
    //    - Automatic parallel transfer at processor patches
    //    - Callback to hitWallPatch() / hitPatch() at boundaries
    //
    //  The deposition field is written to directly by each particle
    //  during tracking - only to local cells, so no reduce needed.
    // ==================================================================

    DynamicList<label> finishedRayIDs;
    DynamicList<DynamicList<point>> finishedRayPaths;
    label maxBounceTerminations = 0;
    scalar maxBounceResidualPower = 0.0;
    scalar minAbsorptivity = GREAT;
    label nearZeroAbsorptivityCount = 0;
    label maxObservedBounces = 0;
    FixedList<scalar, 4> powerSplit(Zero);

    laserRayParticle::trackingData td
    (
        rayCloud,
        alphaFiltered,
        nFiltered,
        resistivity_in,
        deposition_,
        rayQ_,
        rayNumber_,
        dep_cutoff,
        plasma_frequency,
        angular_frequency,
        rayPowerRelTol,
        maxRayBounces,
        maxRayBounceAction,
        nearZeroAbsorptivityTol,
        maxBounceTerminations,
        maxBounceResidualPower,
        minAbsorptivity,
        nearZeroAbsorptivityCount,
        maxObservedBounces,
        powerSplit,
        finishedRayIDs,
        finishedRayPaths
    );

    // Move all ray particles through the mesh.
    // trackTime is not physically meaningful for rays (instantaneous),
    // but the Cloud::move interface requires it. We use GREAT to ensure
    // the rays are tracked to completion in a single call.
    rayCloud.storeGlobalPositions();
    rayCloud.move(rayCloud, td, GREAT);

    reduce(maxBounceTerminations, sumOp<label>());
    reduce(maxBounceResidualPower, sumOp<scalar>());
    reduce(minAbsorptivity, minOp<scalar>());
    reduce(nearZeroAbsorptivityCount, sumOp<label>());
    reduce(maxObservedBounces, maxOp<label>());
    forAll(powerSplit, i)
    {
        reduce(powerSplit[i], sumOp<scalar>());
    }

    // ==================================================================
    // Step 4: Gather ray path segments to master for VTK output
    //
    // Each ray may produce multiple segments (one per processor it
    // traverses). Segments are stored in trackingData both at
    // processor-boundary crossings and at particle death. We
    // gather all segments and concatenate them per rayID.
    // ==================================================================

    if (Pstream::master())
    {
        rayPaths_[laserID].clear();
        rayPaths_[laserID].setSize(nTotalRays);
    }

    {
        // Gather to master
        List<labelList> allIDs(Pstream::nProcs());
        allIDs[Pstream::myProcNo()] = finishedRayIDs;
        Pstream::gatherList(allIDs);

        List<List<DynamicList<point>>> allPaths(Pstream::nProcs());
        allPaths[Pstream::myProcNo()] = finishedRayPaths;
        Pstream::gatherList(allPaths);

        if (Pstream::master())
        {
            for (label procI = 0; procI < Pstream::nProcs(); ++procI)
            {
                const labelList& ids = allIDs[procI];
                const List<DynamicList<point>>& paths = allPaths[procI];

                forAll(ids, i)
                {
                    const label rayID = ids[i];
                    if (rayID >= 0 && rayID < rayPaths_[laserID].size())
                    {
                        DynamicList<point>& existing =
                            rayPaths_[laserID][rayID];
                        const DynamicList<point>& segment = paths[i];

                        if (existing.empty())
                        {
                            // First segment for this ray
                            existing = segment;
                        }
                        else if (segment.size() > 0)
                        {
                            // Append segment, skipping the first
                            // point if it duplicates the junction
                            label startJ = 0;
                            if
                            (
                                mag(existing.last() - segment[0])
                              < SMALL
                            )
                            {
                                startJ = 1;
                            }
                            for (label j = startJ; j < segment.size(); ++j)
                            {
                                existing.append(segment[j]);
                            }
                        }
                    }
                }
            }
        }
    }

    deposition_.correctBoundaryConditions();

    const scalar TotalQ = fvc::domainIntegrate(deposition_).value();
    Info<< "    Total Q deposited: " << TotalQ << endl;

    // Diagnostic: where this laser's power went (fractions of laser power).
    // 'not absorbed' = reflected out of the domain, discarded, or cut off
    // by rayPowerRelTol.
    {
        const scalar P = mag(currentLaserPower) + VSMALL;
        const scalar absorbed =
            powerSplit[0] + powerSplit[1] + powerSplit[2] + powerSplit[3];

        Info<< "    Laser power split (fraction of " << currentLaserPower
            << " W): fresnel " << powerSplit[0]/P
            << "  bulkDump " << powerSplit[1]/P
            << "  escapeDump " << powerSplit[2]/P
            << "  wall " << powerSplit[3]/P
            << "  totalAbsorbed " << absorbed/P
            << "  notAbsorbed " << 1.0 - absorbed/P
            << endl;
    }

    if (maxBounceTerminations > 0)
    {
        Info<< "    Ray-tracing safety diagnostics:" << nl
            << "        maxRayBounces = " << maxRayBounces << nl
            << "        maxRayBounceAction = " << maxRayBounceAction << nl
            << "        max observed bounces = " << maxObservedBounces << nl
            << "        rays handled by maxRayBounces = "
            << maxBounceTerminations << nl
            << "        residual power handled by maxRayBounces = "
            << maxBounceResidualPower << nl
            << "        residual fraction of current laser power = "
            << maxBounceResidualPower/(mag(currentLaserPower) + VSMALL) << nl
            << "        minimum Fresnel absorptivity encountered = "
            << (minAbsorptivity == GREAT ? scalar(-1) : minAbsorptivity)
            << nl
            << "        near-zero absorptivity interactions = "
            << nearZeroAbsorptivityCount << nl
            << "        near-zero absorptivity tolerance = "
            << nearZeroAbsorptivityTol << endl;
    }
}


void laserHeatSource::writeRayPathsToVTK()
{
    const Time& runTime = deposition_.time();
    if (Pstream::master()) // Create a directory for the VTK files
    {
        fileName vtkDir;
        if (Pstream::parRun())
        {
            vtkDir = runTime.path()/".."/"VTKs";
        }
        else
        {
            vtkDir = runTime.path()/"VTKs";
        }

        mkDir(vtkDir);

        vtkTimes_.insert(runTime.value());

        forAll(laserNames_, laserI)
        {
            const fileName vtkFileName
            (
                vtkDir/"rays_" + laserNames_[laserI] + "_"
              + Foam::name(runTime.value()) + ".vtk"
            );
            Info<< "Writing " << rayPaths_[laserI].size()
                << " ray paths to " << vtkFileName << endl;
            writeRayPathsToVTK(rayPaths_[laserI], vtkFileName);
        }
    }
}


void laserHeatSource::writeRayPathsToVTK
(
    const List<DynamicList<point>>& rays,
    const fileName& filename
)
{
    OFstream file(filename);

    if (!file.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << filename
            << exit(FatalError);
    }

    label totalPoints = 0;
    label totalLines = 0;

    forAll(rays, rayI)
    {
        totalPoints += rays[rayI].size();
        if (rays[rayI].size() > 1)
        {
            totalLines += (rays[rayI].size() - 1);
        }
    }

    file<< "# vtk DataFile Version 3.0" << nl;
    file<< "Multiple ray data" << nl;
    file<< "ASCII" << nl;
    file<< "DATASET POLYDATA" << nl;

    file<< "POINTS " << totalPoints << " float" << nl;
    forAll(rays, rayI)
    {
        const DynamicList<point>& ray = rays[rayI];
        forAll(ray, pointI)
        {
            const point& pt = ray[pointI];
            file<< pt.x() << " " << pt.y() << " " << pt.z() << nl;
        }
    }

    file<< "LINES " << totalLines << " " << (totalLines * 3) << nl;

    label pointOffset = 0;
    forAll(rays, rayI)
    {
        const DynamicList<point>& ray = rays[rayI];

        for (label i = 0; i < ray.size() - 1; i++)
        {
            file<< "2 " << (pointOffset + i) << " "
                << (pointOffset + i + 1) << nl;
        }

        pointOffset += ray.size();
    }
}


void laserHeatSource::writeRayPathVTKSeriesFile() const
{
    const Time& runTime = deposition_.time();
    if (Pstream::master())
    {
        const fileName vtkDir =
            Pstream::parRun()
          ? runTime.path()/".."/"VTKs"
          : runTime.path()/"VTKs";

        const scalarList times = vtkTimes_.toc();

        SortableList<scalar> sortedTimes(times);

        forAll(laserNames_, laserI)
        {
            const word& laserName = laserNames_[laserI];

            const fileName seriesFile =
                vtkDir/"rays_" + laserName + ".vtk.series";
            Info<< "Writing ray path series file: " << seriesFile << nl;

            OFstream os(seriesFile);
            os.precision(12);

            os  << "{\n"
                << "  \"file-series-version\": \"1.0\",\n"
                << "  \"files\": [\n";

            for (label i = 0; i < sortedTimes.size(); ++i)
            {
                const scalar t = sortedTimes[i];

                os  << "    { \"name\": \"rays_"
                    << laserName << "_" << Foam::name(t)
                    << ".vtk\", \"time\": " << t << " }";

                if (i+1 < sortedTimes.size()) os << ",";
                os << "\n";
            }

            os  << "  ]\n"
                << "}\n";
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
