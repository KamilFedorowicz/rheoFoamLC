#ifndef LE_1constant_H
#define LE_1constant_H
#include "constitutiveEq2.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace constitutiveEqs
{ 
/*---------------------------------------------------------------------------*\
                           Class LE_1constant Declaration
\*---------------------------------------------------------------------------*/
class LE_1constant
:
    public constitutiveEq2
{
    // Private data
        //- Transported non-Newtonian nematic stress
        volTensorField tau_;
        //- director field
        volVectorField n_;
        //- Transported total nematic stress
        volTensorField tauTotal_;

        // Model constants
            dimensionedScalar rho_;
            dimensionedScalar etaS_;
            dimensionedScalar alpha1_;
            dimensionedScalar alpha2_;
            dimensionedScalar alpha3_;
            dimensionedScalar alpha5_;
            dimensionedScalar alpha6_;
            dimensionedScalar K_;
 
    // Private Member Functions
        //- Disallow default bitwise copy construct
        LE_1constant(const LE_1constant&);
        //- Disallow default bitwise assignment
        void operator=(const LE_1constant&);
        
protected:

       //- Return the solvent viscosity
       virtual const dimensionedScalar etaS() const
       {
          return etaS_;
       }
       //- Return alpha1
       virtual const dimensionedScalar alpha1() const
       {
          return alpha1_;
       }
       //- Return alpha2
       virtual const dimensionedScalar alpha2() const
       {
          return alpha2_;
       }
       //- Return alpha3
       virtual const dimensionedScalar alpha3() const
       {
          return alpha3_;
       }
       //- Return alpha5
       virtual const dimensionedScalar alpha5() const
       {
          return alpha5_;
       }
       //- Return alpha6
       virtual const dimensionedScalar alpha6() const
       {
          return alpha6_;
       }

public:
    //- Runtime type information
    TypeName("LE_1constant");

    // Constructors
        //- Construct from components
        LE_1constant
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );
    // Destructor
        virtual ~LE_1constant()
        {}
    // Member Functions
       //- Return the viscoelastic stress tensor
        virtual tmp<volTensorField> tau() const
        {
            return tau_;
        }   
        //- Return the density
        virtual const dimensionedScalar rho() const
        {
            return rho_;
        }     
        //- Return true if GNF (non-elastic)
        virtual bool isGNF() const
        {
          return false;
        };

        //- Correct the viscoelastic stress
        virtual void correct();
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace constitutiveEqs 
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// ************************************************************************* //
