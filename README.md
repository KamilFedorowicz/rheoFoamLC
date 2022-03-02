# rheoFoamLC
the repository contains an OpenFOAM solver dedicated to modelling the flow of liquid crystals

The full solver procedure is contain within the file rheoFoamLC.C. 
The header file createFields.H declares/initialises the relevant fields and variables: velocity, pressure, microstructure and stress.
The latter two fields are specified through the constitutiveModel2 object; 
in case of the Leslie-Ericksen model the microstructure is described through the director field, 
while in tensorial models, the microstructure is described through the vb{Q}-tensor. 
Finally, the specific form of the constitutive equation is contained in the object LE_1constant, 
where the director evolution equation and stress definition are specified.

