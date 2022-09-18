# Pullout of neo-Hookean fiber embedded in a generalized neo-Hookean matrix
Finite element implementation in FEniCS. This code was used for the simulations in the paper M. Myneni, K. Tůma, P. Kara, K.R. Rajagopal, C.C. Benjamin: "Pullout of neo-Hookean fiber embedded in a generalized neo-Hookean matrix".

The code is based on FEniCS 2019.1.0 (https://fenicsproject.org/) using Python3 environment. It consists of several files:
* axisym_neoHooke_case1.py, axisym_neoHooke_case2.py, axisym_neoHooke_case3.py: scripts corresponding to three different cases
* meshes/mesh1.geo, meshes/mesh2.geo: geometry scripts for description of two cases
* meshes/data/mesh1.msh, meshes/data/mesh2.msh: GMSH files of meshes corresponding to mesh1.geo and mesh2.geo

To run the simulation use <code>python3 axisym_neoHooke_case1.py muf n [kapparatio] [b]</code>, where muf and n are two compulsory parameters and kapparatio and b are optional parameters, muf is the ratio muf:mum, n is the exponent in the power-law model of the matrix, kapparatio is the ratio kappa:mu (default kapparatio = 100), b is the last material parameter of the matrix (default  b = 1).

In the output directory it generates several files:
* u.xdmf, u.h5; norme.xdmf, norme.h5; TmRZ.xdmf, TmRZ.h5; TfRZ.xdmf, TfRZ.h5; TRZ.xdmf, TRZ.h5: HDF5 data for displacement, norm of the Green strain tensor, shear component of the Cauchy stress tensor in the matrix, fiber and its combination readable by Paraview (https://www.paraview.org/).

For support please contact Karel Tůma (ktuma@karlin.mff.cuni.cz).
