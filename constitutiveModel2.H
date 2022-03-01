#ifndef constitutiveModel2_H
#define constitutiveModel2_H

#include "constitutiveEq2.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class constitutiveModel2 Declaration
\*---------------------------------------------------------------------------*/

class constitutiveModel2
:
    public IOdictionary
{
private:

    // Private data

        //- Constitutive Eq
        autoPtr<constitutiveEq2> eqPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        constitutiveModel2(const constitutiveModel2&);

        //- Disallow default bitwise assignment
        void operator=(const constitutiveModel2&);


public:

    //- Runtime type information
    TypeName("constitutiveModel2");


    // Constructors

        
        //- Construct from components
        constitutiveModel2
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    // Destructor

        virtual ~constitutiveModel2()
        {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volTensorField> tau() const;
        
        //- Return the total stress tensor (solvent + polymeric)
        virtual tmp<volTensorField> tauTotal() const;
        
        //- Return the density
        virtual const dimensionedScalar rho() const;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;
        
        //- Is GNF or VE
        virtual bool isGNF() const;
       
        //- Correct the viscoelastic stress
        virtual void correct();

        //- Read constitutiveProperties dictionary
        virtual bool read();

         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
