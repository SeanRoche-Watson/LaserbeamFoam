#include "createRDeltaT.H"
#include "Function1.H" // added Function1 header file

const dimensionedScalar deltaN = 1e-8/pow(average(mesh.V()), 1.0/3.0);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Temperature\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const word interfaceTrackingScheme
(
    transportProperties.lookup("interfaceTrackingScheme")
);

if
(
    interfaceTrackingScheme != "MULES"
 && interfaceTrackingScheme != "isoAdvector"
)
{
    FatalErrorIn("laserbeamFoam.C")
        << "Unknown interface tracking scheme "
        << interfaceTrackingScheme
        << ". Available options are 'MULES' and 'isoAdvector'"
        << exit(FatalError);
}

dictionary& transportPropertiesMetal = transportProperties.subDict("metal");
dictionary& transportPropertiesGas = transportProperties.subDict("gas");


const dimensionedScalar Tsolidus1
(
    "Tsolidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesMetal//.lookup("Tsolidus")
);
// Info<<Tsolidus1<<endl;
const dimensionedScalar Tsolidus2
(
    "Tsolidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesGas//.lookup("Tsolidus")
);
// Info<<Tsolidus2<<endl;

const dimensionedScalar Tliquidus1
(
    "Tliquidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesMetal//.lookup("Tliquidus")
);

const dimensionedScalar Tliquidus2
(
    "Tliquidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesGas//.lookup("Tliquidus")
);

const dimensionedScalar LatentHeat1
(
    "LatentHeat",
    dimensionSet(0, 2, -2, 0, 0),
    transportPropertiesMetal//.lookup("LatentHeat")
);

const dimensionedScalar LatentHeat2
(
    "LatentHeat",
    dimensionSet(0, 2, -2, 0, 0),
    transportPropertiesGas//.lookup("LatentHeat")
);

const dimensionedScalar beta1
(
    "beta",
    dimensionSet(0, 0, 0, -1, 0),
    transportPropertiesMetal//.lookup("beta")
);

const dimensionedScalar beta2
(
    "beta",
    dimensionSet(0, 0, 0, -1, 0),
    transportPropertiesGas//.lookup("beta")
);

const dimensionedScalar Marangoni_Constant
(
    "Marangoni_Constant",
    dimensionSet(1, 0, -2, -1, 0),
    transportProperties//.lookup("dsigmadT")
);

// Temperature-dependent surface-tension gradient (optional).
//   MarangoniModel constant;  -> dsigma/dT = Marangoni_Constant (default, as before)
//   MarangoniModel Sahoo;     -> Fe-S model of Sahoo, DebRoy & McNallan (1988):
//
//     dsigma/dT = -A - R*Gamma_s*ln(1 + K*a_s)
//                 - K*a_s/(1 + K*a_s) * Gamma_s*deltaH0/T,
//     K = k1*exp(-deltaH0/(R*T))
//
//   a_s = sulphur activity (taken as the S content in wt%), read as
//   MarangoniSulphur. The other constants default to the Fe-S values of Sahoo
//   et al. and can be overridden in transportProperties:
//     MarangoniA        4.3e-4    [N/m/K]
//     MarangoniGammaS   1.3e-8    [kmol/m^2]
//     MarangoniDeltaH0  -1.66e8   [J/kmol]
//     MarangoniK1       3.18e-3   [-]
//   With S above ~30 ppm dsigma/dT is positive (inward flow) below a
//   cross-over temperature and negative (outward) above it.
const word MarangoniModel
(
    transportProperties.lookupOrDefault<word>("MarangoniModel", "constant")
);

const bool sahooMarangoni = (MarangoniModel == "Sahoo");

if (!sahooMarangoni && MarangoniModel != "constant")
{
    FatalErrorInFunction
        << "Unknown MarangoniModel " << MarangoniModel
        << ". Valid options: constant, Sahoo" << exit(FatalError);
}

const scalar MarangoniSulphur
(
    sahooMarangoni
  ? transportProperties.get<scalar>("MarangoniSulphur")
  : 0.0
);
const scalar MarangoniA
(
    transportProperties.lookupOrDefault<scalar>("MarangoniA", 4.3e-4)
);
const scalar MarangoniGammaS
(
    transportProperties.lookupOrDefault<scalar>("MarangoniGammaS", 1.3e-8)
);
const scalar MarangoniDeltaH0
(
    transportProperties.lookupOrDefault<scalar>("MarangoniDeltaH0", -1.66e8)
);
const scalar MarangoniK1
(
    transportProperties.lookupOrDefault<scalar>("MarangoniK1", 3.18e-3)
);

// Sahoo dsigma/dT at temperature T [N/m/K]
auto sahooDsigmadT = [&](const scalar Tc) -> scalar
{
    const scalar R = 8314.462618;   // J/(kmol K)
    const scalar K =
        MarangoniK1*Foam::exp(-MarangoniDeltaH0/(R*Tc));
    const scalar Ka = K*MarangoniSulphur;
    return
       -MarangoniA
      - R*MarangoniGammaS*Foam::log(1.0 + Ka)
      - Ka/(1.0 + Ka)*MarangoniGammaS*MarangoniDeltaH0/Tc;
};

volScalarField dsigmadT
(
    IOobject
    (
        "dsigmadT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    Marangoni_Constant
);

Info<< "MarangoniModel: " << MarangoniModel;
if (sahooMarangoni)
{
    Info<< "  S = " << MarangoniSulphur << " wt%" << nl
        << "    dsigma/dT at 1750 / 2000 / 2500 / 3000 K = "
        << sahooDsigmadT(1750) << " / " << sahooDsigmadT(2000) << " / "
        << sahooDsigmadT(2500) << " / " << sahooDsigmadT(3000)
        << " N/m/K" << endl;
}
else
{
    Info<< "  dsigma/dT = " << Marangoni_Constant.value() << " N/m/K" << endl;
}

const dimensionedScalar p0
(
    "p0",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    transportProperties//.lookup("p0")
); // atmospheric pressure

const dimensionedScalar Tvap
(
    "Tvap",
    dimensionSet(0, 0, 0, 1, 0),
    transportProperties//.lookup("Tvap")
);

const dimensionedScalar Mm
(
    "Mm",
    dimensionSet(1, 0, 0, 0, -1),
    transportProperties//.lookup("Mm")
);

const dimensionedScalar LatentHeatVap
(
    "LatentHeatVap",
    dimensionSet(0, 2, -2, 0, 0),
    transportProperties//.lookup("LatentHeatVap")
);

const dimensionedScalar R
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
        8.314
);

volScalarField em
(
    IOobject
    (
        "em",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

const dimensionedScalar SB
(
    "SB",
    dimensionSet(1, 0, -3, -4, 0, 0, 0),
    scalar(5.67e-8)
); // Stefan-Boltzmann Constant

[[maybe_unused]]const scalar emL(readScalar(transportProperties.lookup("emL")));
[[maybe_unused]]const scalar emS(readScalar(transportProperties.lookup("emS")));
const dimensionedScalar TRef ("TRef",transportProperties);

// Deactivated polynomials
/*
const Polynomial<8> polykappa_m(transportPropertiesMetal.lookup("poly_kappa"));
const Polynomial<8> polykappa_g(transportPropertiesGas.lookup("poly_kappa"));

const Polynomial<8> polycp_m(transportPropertiesMetal.lookup("poly_cp"));
const Polynomial<8> polycp_g(transportPropertiesGas.lookup("poly_cp"));
*/

// added piecewise linear interpolated tabulated data based on Function1

autoPtr<Function1<scalar>> kappa_m = 
	Function1<scalar>::New("kappaTable", transportPropertiesMetal);

autoPtr<Function1<scalar>> kappa_g = 
	Function1<scalar>::New("kappaTable", transportPropertiesGas);

autoPtr<Function1<scalar>> cp_m = 
	Function1<scalar>::New("cpTable", transportPropertiesMetal);

autoPtr<Function1<scalar>> cp_g = 
	Function1<scalar>::New("cpTable", transportPropertiesGas);



volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//change to no write when working
    ),
    mesh,
    dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0.0)
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//change to no write when working
    ),
    mesh,
    dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0)
);

#include "updateKappaCp.H"


volScalarField TSolidus
(
    IOobject
    (
        "TSolidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*Tsolidus1 + (1.0 - alpha1)*Tsolidus2
);


volScalarField TLiquidus
(
    IOobject
    (
        "TLiquidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*Tliquidus1 + (1.0 - alpha1)*Tliquidus2
);

volScalarField DC
(
    IOobject
    (
        "DC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DC", dimensionSet(1,-3,-1,0,0), 1.0e14)
);

volScalarField metalaverage
(
    IOobject
    (
        "metalaverage",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::average(alpha1)
);

volScalarField epsilon1
(
    IOobject
    (
        "epsilon1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    max(min((T - TSolidus)/(TLiquidus - TSolidus), 1.0), 0.0)
);

volVectorField n_filtered
(
    IOobject
    (
        "n_filtered",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("n_filtered", dimless, vector::zero)//,
    // zeroGradientFvPatchVectorField::typeName
);

volScalarField epsilon1mask
(
    IOobject
    (
        "epsilon1mask",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask",dimless,0.0)
);

volScalarField epsilon1mask2
(
    IOobject
    (
        "epsilon1mask2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask2", dimless, 0.0)
);

volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradT", dimensionSet(0, -1, 0, 1, 0), vector::zero)
);

volScalarField LatentHeat
(
    IOobject
    (
        "LatentHeat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*LatentHeat1 + (1.0 - alpha1)*LatentHeat2
);

volScalarField TRHS
(
    IOobject
    (
        "TRHS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("TRHS", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

volScalarField Tcorr
(
    IOobject
    (
        "Tcorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Tcorr", dimensionSet(0, 0, 0, 1, 0), 0.0)
);

volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*beta1 + (1.0-alpha1)*beta2
);

volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhok", dimless, 0.0)
);

volScalarField Num_divU
(
    IOobject
    (
        "Num_divU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Num_divU", dimensionSet(0,0,-1,0,0), 0.0)
);

volVectorField Marangoni
(
    IOobject
    (
        "Marangoni",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Marangoni", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

const volVectorField zeroJcrossB
(
    IOobject
    (
        "zeroJcrossB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zeroJcrossB", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

const volScalarField zeroJoulesource
(
    IOobject
    (
        "zeroJoulesource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroJoulesource", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

// Recoil Pressure
volScalarField pVap
(
    IOobject
    (
        "pVap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

// Evaporation Cooling
volScalarField Qv
(
    IOobject
    (
        "Qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
);


volScalarField damper
(
    IOobject
    (
        "damper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("damper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField thermalDamper
(
    IOobject
    (
        "thermalDamper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalDamper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField alpha_filtered
(
    IOobject
    (
        "alpha_filtered",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_filtered", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

volScalarField alpha_smoothed
(
    IOobject
    (
        "alpha_smoothed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_smoothed", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

volScalarField e1temp
(
    IOobject
    (
        "e1temp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("e1temp", dimensionSet(0, 0, 0, -0, 0), 0.0)
);
e1temp = fvc::average(epsilon1);



volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoCp", dimensionSet(1, -1, -2, -1, 0), 0.0)
);
rhoCp = (rho*cp);
// volScalarField rhoCp(rho*cp);
rhoCp.oldTime();


// surfaceScalarField rhophicp(fvc::interpolate(cp)*rhoPhi);

surfaceScalarField rhophicp
(
    IOobject
    (
        "rhophicp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhophicp", dimensionSet(1, 2, -3, -1, 0), 0.0)
);
rhophicp = fvc::interpolate(cp)*rhoPhi;



const dimensionedScalar DarcyConstantlarge
(
    "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0), 1.0e6
);

const dimensionedScalar DarcyConstantsmall
(
    "DarcyConstantsmall", dimless, 1.0e-12
);


// volVectorField gradAlpha(fvc::grad(alpha1));

volVectorField gradAlpha
(
    IOobject
    (
        "gradAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradAlpha", dimensionSet(0, -1, 0, -0, 0), vector::zero)
);
gradAlpha = fvc::grad(alpha1);

// volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN));


volVectorField nHatM
(
    IOobject
    (
        "nHatM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("nHatM", dimensionSet(0, 0, 0, -0, 0), vector::zero)
);
nHatM = gradAlpha/(mag(gradAlpha) + deltaN);


// This variable is 1 if (alpha.metal >= 0.5 and epsilon1 >=0.5), 0 otherwise.
// These are the cells that are liquid and metallic at the current time
volScalarField liquidMetalCells
(
    IOobject
    (
        "liquidMetalCells",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("liquidMetalCells", dimensionSet(0, 0, 0, -0, 0), 0.0)
);


typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"




// ***** 
// Electrical resistivity modifications, to be temperature dependent
// *****
autoPtr<Function1<scalar>> electricalResistivity_m =
	Function1<scalar>::New
	(
		"electricalResistivityTable",
		transportPropertiesMetal
	);
	
/*
dimensionedScalar elec_resistivity_metal
(
    "elec_resistivity",
    dimensionSet(1, 3, -3, 0, 0),
    transportPropertiesMetal
    //.lookup("elec_resistivity")
);
*/

// Creating the field electrical_resistivity
volScalarField electrical_resistivity
(
    IOobject
    (
        "electrical_resistivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    //elec_resistivity_metal
	dimensionedScalar
	(
		"electricalResistivity",
		dimensionSet(1, 3, -3, 0, 0, -2, 0),
		scalar(0)
	) 
);
// end of modifications

// Weight the metal resistivity by alpha.metal in mixed cells?
//   no  (default): the laser sees the metal's own rho_e(T) in interface cells
//   yes          : rho_e = alpha.metal*rho_e(T)  (original PR #140 behaviour)
const Switch resistivityAlphaWeighting
(
    transportPropertiesMetal.getOrDefault<Switch>
    (
        "resistivityAlphaWeighting",
        false
    )
);
Info<< "resistivityAlphaWeighting = " << resistivityAlphaWeighting << endl;



autoPtr<mthdModel> mthd(nullptr);

if (mthdModel::electromagneticsEnabled(mesh, "transportProperties"))
{
    mthd.reset(new mthdModel(mesh, "transportProperties"));
}

const dimensionedScalar eleccond1 =
    mthdModel::readElectricalConductivity(transportPropertiesMetal);


// (Original constant electrical_resistivity field removed: it is now
//  declared above and filled from electricalResistivityTable each step.)



// Temperature dependent viscosity implementation

autoPtr<Function1<scalar>> nu_m = 
	Function1<scalar>::New("nuTable", transportPropertiesMetal);

autoPtr<Function1<scalar>> nu_g = 
	Function1<scalar>::New("nuTable", transportPropertiesGas);
	
dimensionedScalar nuTable0
(
	"nuTable0",
	dimViscosity,
	0.0
);


volScalarField nuTableField
(
    IOobject
    (
        "nuTableField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    nuTable0
);	

// end of modification of nu



// Create laser heat source
laserHeatSource laser(mesh);

#include "updateYDim.H"
#include "updateProps.H"

// volScalarField depositionTemp
// (
//     IOobject
//     (
//         "depositionTemp",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("depositionTemp", dimensionSet(1, -1, -2, -1, 0), 0.0),
//     zeroGradientFvPatchScalarField::typeName
// );


// depositionTemp=laser.deposition();
// depositionTemp.correctBoundaryConditions();

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh*rhok
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh*rhok;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"


#include "isoAdvector/createPorosity.H"
isoAdvection advector(alpha1, phi, U);

// Log how much of the interface gradient survives the gradSchemes limiter in
// molten interface cells (Marangoni force uses grad(alpha.metal) and grad(T)).
// Compares the gradients the solver uses with an unlimited Gauss linear reference.
const Switch gradientDiagnostics
(
    transportProperties.getOrDefault<Switch>("gradientDiagnostics", true)
);
Info<< "gradientDiagnostics = " << gradientDiagnostics << endl;
