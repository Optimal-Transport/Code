[![DOI](https://zenodo.org/badge/305434774.svg)](https://zenodo.org/badge/latestdoi/305434774)

# Code files for Optimal Transport Taster Project

The repo has the following folders:
## Exact
A set of instances from the DOTmark base solved using linear programming methods (namely simplex and barrier methods) with different costs. The generated OT problems were solved using the state of the art solver Gurobi.

## Python Algorithms Folder
This folder contains these files:

* ```ℓ1-Projection.ipynb``` $\ell_1$ Projection algorithm from Condat https://doi.org/10.1007/s10107-015-0946-6
* ```Primal Dual and Generalised Forward Backward.ipynb``` DOT solver using the Primal Dual and Generalised Forward Backward algorithm.
* ```FISTA_backtracking.ipynb``` DOT solver using the FISTA algorithm with backtracking.
* ```OT - Linear Programming.ipynb``` Linear programming code for minimal example using Scipy and Gurobi.
* ```README.md``` This file.
