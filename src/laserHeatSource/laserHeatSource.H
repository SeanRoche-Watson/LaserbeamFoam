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

Class
    laserHeatSource

Description
    Calculate the deposition heat source field from a laser using ray-tracing.

    Uses OpenFOAM's Lagrangian particle tracking infrastructure for exact
    cell-to-cell ray propagation with automatic parallel processor-boundary
    transfer.

    Laser parameters are defined via the provided dictionary.

    It is assumed that the laser enters the domain at the boundary with an
    outward-facing unit normal of (0 1 0).

    Laser rays are written as VTKs for visualisation.

Author
    Tom Flint, The University of Manchester.
    Philip Cardiff, University College Dublin.

SourceFiles
    laserHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef laserHeatSource_H
#define laserHeatSource_H

#include "volFields.H"
#include "dictionary.H"
#include "interpolationTable.H"
#include "OFstream.H"
#include "Cloud.H"
#include "laserRayParticle.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                        Class laserHeatSource Declaration
\*---------------------------------------------------------------------------*/


class laserHeatSource
:
    public IOdictionary
{
    // Private data

        // Heat source field
        volScalarField deposition_;

        //- Field indicating the boundary faces where the laser enters the
        //  domain
        volScalarField laserBoundary_;

        //- Field for tracking error codes
        volScalarField errorTrack_;

        //- Field for tracking the latest ray to vist each cell
        volScalarField rayNumber_;

        //- Field for tracking the energy deposited in each cell
        volScalarField rayQ_;

        //- y dimension field
        volScalarField yDim_;

        //- Powder sim flag
        const Switch powderSim_;

        //- New laser discretisation flag
        const Switch radialPolarHeatSource_;

        //- List of the laser names
        wordList laserNames_;

        //- Copy of the laser dicts
        PtrList<dictionary> laserDicts_;

        //- List of the time vs laser position interpolation series
        PtrList<interpolationTable<vector>> timeVsLaserPosition_;

        //- List of the time vs laser power interpolation series
        PtrList<interpolationTable<scalar>> timeVsLaserPower_;

        //- Record the ray paths for each laser
        List<List<DynamicList<point>>> rayPaths_;

        //- Record of times when VTK files were written
        HashSet<scalar> vtkTimes_;

        //- Global bounding box of the mesh
        boundBox globalBB_;

        //- Optional multi-term beam profile for the laser currently being
        //  seeded (empty -> single super-Gaussian of radius laserRadius).
        //  Each entry: (type, fraction, a, b) with
        //    type 0 = gaussian: a = 1/e^2 radius w, b unused
        //    type 1 = ring:     a = ring radius R, b = Gaussian sigma s
        //  Fractions are normalised to sum to 1.
        List<FixedList<scalar, 4>> beamTerms_;

        //- Tabulated radial beam profile (beamController). When set, it
        //  overrides beamTerms / laserRadius: I(r) = P*G(r), G piecewise
        //  linear between the nodes ctrlR_ [m], zero beyond the last node.
        //  G has units 1/m^2; its plane integral is the power scale factor.
        scalarField ctrlR_;
        scalarField ctrlG_;

        //- Controller bookkeeping
        bool ctrlInit_;
        scalar ctrlLastUpdate_;

        //- Fraction of the tabulated laser power that reaches the target
        //  (plasma-shielding transmission). 1 = original behaviour. Set by
        //  the solver each time step via setTransmission().
        scalar transmission_;


    // Private Member Functions

        //- Beam-profile controller / fixed tabulated profile (see the
        //  beamController sub-dictionary documented in laserHeatSource.C)
        void updateBeamController
        (
            const dictionary& dict,
            const vector& laserPosition,
            const volScalarField& alphaFiltered,
            const scalar laserRadius,
            const scalar Radius_Flavour,
            const scalar superGaussianOrder
        );

        //- Seed the ray cloud with initial ray particles.
        //  Creates laserRayParticle objects and adds them to the cloud.
        //  Only creates particles on the processor that owns the starting
        //  cell. Returns the total number of rays across all processors
        //  via nTotalRays.
        void seedRayCloud
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
        ) const;

        //- Disallow default bitwise copy construct
        laserHeatSource(const laserHeatSource&);

        //- Disallow default bitwise assignment
        void operator=(const laserHeatSource&);

public:

    //- Runtime type information
    TypeName("laserHeatSource");

    // Constructors

        //- Construct from the mesh
        laserHeatSource(const fvMesh& mesh);

    // Destructor

        virtual ~laserHeatSource()
        {}


    // Member Functions

        //- Const access to the deposition field
        const volScalarField& deposition() const
        {
            return deposition_;
        }

        //- Non-const access to the yDim field
        volScalarField& yDim()
        {
            return yDim_;
        }

        //- Const access to the powderSim flag
        const Switch& powderSim() const
        {
            return powderSim_;
        }

        //- Const access to the HS discretisation flag
        const Switch& radialPolarHeatSource() const
        {
            return radialPolarHeatSource_;
        }

        //- Update the deposition field for all lasers
        //  The function takes three arguments
        //  The filtered volume of metal field  (alphaFiltered), the normal
        //  vector field nFiltered and the scalar resistivity field
        void updateDeposition
        (
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in
        );

        //- Update the deposition for the given laser parameters
        void updateDeposition
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
        );

        //- Set the transmitted fraction of the laser power (clipped to
        //  [0, 1]); applied to every laser at the next updateDeposition
        void setTransmission(const scalar tau)
        {
            transmission_ = Foam::min(Foam::max(tau, scalar(0)), scalar(1));
        }

        //- Current transmitted fraction of the laser power
        scalar transmission() const
        {
            return transmission_;
        }

        //- Const access to the global bounding box
        const boundBox& globalBB() const
        {
            return globalBB_;
        }

        //- Write the given ray paths to a VTK file with the given name
        void writeRayPathsToVTK
        (
            const List<DynamicList<point>>& rays,
            const fileName& filename
        );

        //- Write rayPath_ to a VTK file with the given name
        void writeRayPathsToVTK();

        //- Write a "vtk.series" file
        void writeRayPathVTKSeriesFile() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
