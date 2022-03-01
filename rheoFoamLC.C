#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "simpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

#include "adjustCorrPhi.H"
#include "ppUtilInterface.H"
#include "constitutiveModel2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createControls.H" 
    #include "createPPutil.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Inner loop iterations ---
        for (int i=0; i<nInIter; i++)
        {
            Info << "Inner iteration:  " << i << nl << endl; 
            
            if (i==0 || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // --- Pressure-velocity SIMPLEC corrector
            {
               // ---- Solve U and p ----	
               #include "UEqn.H"
               #include "pEqn.H"         
            }
            
            // ---- Solve constitutive equation ----	
            constEq2.correct();

            // --- Passive Scalar transport
            if (sPS)
             {
               #include "CEqn.H"
             }            
        }

        postProc.update();
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
