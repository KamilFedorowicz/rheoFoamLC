#include "LE_1constant.H"
#include "addToRunTimeSelectionTable.H"
 
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(LE_1constant, 0);
    addToRunTimeSelectionTable(constitutiveEq2, LE_1constant, dictionary);
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// declare stress (tau_), director (n_) and the total stress (tauTotal_)
Foam::constitutiveEqs::LE_1constant::LE_1constant
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    constitutiveEq2(name, U, phi),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    n_
    (
        IOobject
        (
            "n" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    tauTotal_
    (
        IOobject
        (
            "tauTotal" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

// alpha_i are Leslie viscosities
// etaS=alpha4/2 is the Newtonian viscosity
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    alpha1_(dict.lookup("alpha1")),
    alpha2_(dict.lookup("alpha2")),
    alpha3_(dict.lookup("alpha3")),
    alpha5_(dict.lookup("alpha5")),
    alpha6_(dict.lookup("alpha6")),
    K_(dict.lookup("K"))


{
 checkForStab(dict);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::LE_1constant::correct()
{
	// Velocity gradient tensor
	volTensorField L = fvc::grad(U());
	// Symmetric velocity gradient tensor
	volTensorField D = ( L+T(L) )/2;
	// Antisymmetric velocity gradient tensor
	volTensorField omega = skew(L);

//director evolution equation
    fvVectorMatrix nEqn
    (
	(alpha3_ - alpha2_)*(fvm::ddt(n_) + fvm::div(phi(), n_) - (n_ & omega) )
	- K_*fvc::laplacian(n_)
	+ ( (n_ *n_) &   K_*fvc::laplacian(n_)    )
	==
	-(alpha6_ - alpha5_)*( (n_ & D) - n_*tr(n_*n_ & D)  )
    );
	nEqn.relax();
	nEqn.solve();
	n_=n_/mag(n_);

// corotational derivative of the director field
	volVectorField N=(  fvc::ddt(n_) + fvc::div(phi(), n_)  - (n_ & omega)  );

// non-Newtonian stress tensor in the Leslie-Ericksen model
	tau_ = (   
	alpha1_ *tr(D & n_ * n_)*(n_*n_) 
	+ alpha2_ * (     n_ *N  )
	+ alpha3_ * (     N *n_  )
	+ alpha5_ * (	(n_*n_) & D )
	+ alpha6_ * (D & (n_*n_)	)
	- K_ * ( 	fvc::grad(n_) & T(fvc::grad(n_))	)
	);

	// total stress tensor
	tauTotal_=tau_+2*etaS_*D;
  	tau_.correctBoundaryConditions();
}
// ************************************************************************* //
