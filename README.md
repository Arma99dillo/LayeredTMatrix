# LayeredTMatrix
This repository implements the TDG T-matrix method for scattering by periodic layered structures. This repository contains the scripts used to carry out the numerical tests of the method described in the paper

_Trefftz DG Approximation of the T-Matrix for Scattering by Periodic Layered Structures (2026)_, by Armando Maria Monforte, Andrea Moiola and Simone Zanotto

# Reproducibility Instructions
This directory contains all the code necessary to reproduce the numerical results presented in the paper. All the codes have been tested on MATLAB R2024b release. The MATLAB Partial Differential Equation Toolbox™ is needed for the correct working of the code.

# Contents and Structure

Experiment scripts:
-
The following files contain the main code to run all the experiments in Section 6 of the paper, allowing to reproduce the figures therein:
* The files `pConvergenceTwoFlatLayers.m`, `pConvergenceRelative.m` and `MConvergenceRelative.m` are used to derive the _p_- and _M_-convergence plots in the numerical experiments of the paper;
* The file `RepeatedShapes.m` is used to generate the field plots whene there are repeated layers;
* The file `ComputationalTimesComparison.m` compares the TDG T-matrix method performance with the DtN-TDG method available at https://github.com/Arma99dillo/DtN-TDG;
* The file `BraggReflector.m` generates the plots for the distributed Bragg reflector structure used in the comparison tests with the RCWA method.

Other scripts:
-
* The file `SolveSingleProblem.m` is used to run the TDG T-matrix method and solve the scattering problem on a arbitrary chosen domain, plotting only the numerical solution, without convergence or error plots.

Directories:
-
All the auxiliary files are int the `src` directory.
* All the domain configurations used in the paper are implemented in the `GenerateMesh.m` file;
* The files `MatrixDtNTDG.m` and `rhsDtNTDG.m` implement the TDG linear system as described in the paper, both with and withouth impenetrable obstacles inclusions;
* The file `BuildTMatrix.m` builds the layer T-matrix as described in the paper;
* The files `MultiLayerSolve.m` and `MultiLayerSolveSameType.m` build and solve the coupling system used to compute the solution on the full domain;

The `quadtriangle` directory implements the Duffy quadrature rule on triangular elements and is needed to compute numerical errors, and is taken from https://www.mathworks.com/matlabcentral/fileexchange/72131-quadtriangle.

Change the parameters and new configurations
-
You can easily change the problem parameters such as the wavenumber, the incident angle, the mesh width and the domain configuration in the `SolveSingleProblem.m` file and test the convergence changing the parameters in any of the _p_-convergence scripts. It's easy to add custom configurations: the 'penetrable_poly' and 'dir_poly' domains allow to define any penetrable and impenetrable polygon. Otherwise, other geometries can be added in `GenerateMesh.m`, it is only needed to implement a file to generate the mesh and a file to mark the mesh elements depending on the value of the relative permittivity. I suggest to use `MeshDouble.m` and `EpsValDouble.m` as templates.
