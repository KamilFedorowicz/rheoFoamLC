# rheoFoamLC
the repository contains an OpenFOAM solver dedicated to modelling the flow of liquid crystals

Solver rheoFoamLC is based on the existing solver rheoFoam.
Compared to rheoFoam, rheoFoamLC allows for non-symmetric stress tensors,
which often occur in liquid crystal flows.


The full solver procedure is contain within the file rheoFoamLC.C. 
The header file createFields.H declares/initialises the relevant fields and variables: 
velocity, pressure, microstructure and stress.
The latter two fields are specified through the constitutiveModel2 object.

The exact form of the constitutive equation is specified in a dedicated file;
look for the LE_model folder in other repository to see an example.


