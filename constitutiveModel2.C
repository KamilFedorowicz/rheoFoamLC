#include "constitutiveModel2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(constitutiveModel2, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constitutiveModel2::constitutiveModel2
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "constitutiveProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    eqPtr_(constitutiveEq2::New(word::null, U, phi, subDict("parameters")))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volTensorField> constitutiveModel2::tau() const
{
    return eqPtr_->tau();
}

tmp<volTensorField> constitutiveModel2::tauTotal() const
{
    return eqPtr_->tauTotal();
}

const dimensionedScalar constitutiveModel2::rho() const
{
    return eqPtr_->rho();
}

tmp<fvVectorMatrix> constitutiveModel2::divTau(volVectorField& U) const
{
    return eqPtr_->divTau(U);
}

bool constitutiveModel2::isGNF() const
{
    return eqPtr_->isGNF();
}

void constitutiveModel2::correct()
{
    eqPtr_->correct();
}


bool constitutiveModel2::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
