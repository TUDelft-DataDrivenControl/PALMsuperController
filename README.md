## PALMsuperController

Developed by Delft University of Technology, 2018


## Summary:
 
This software package offers a communication protocol between the in FORTRAN written high-fidelity wind farm PArallelized Large-eddy simulation Model (PALM) and MATLAB. It allows control developers to evaluate controllers in a high-fidelity wind farm simulation model. Several examples are included.


## Required software:

1) MATLAB (tested with 2016 version).
2) Complete installation of PALM (tested with rev. 2227). 
3) Solver (tested with CPLEX) and Yalmip (both need to be placed in libraries directory of the controller).


## Description:
MPC:
Deterministic model predictive controller that provides wind farm power tracking.

SMPC:
Stochastic model predictive controller that provides wind farm power tracking.

WRC:
Hinf controller that provides wake redirection control.

MAKE_WTM:
Contains a function that generates excitation signals for PALM. 
 
BIN:
Contains MATLAB function that plots PALM flow data (.nc) at hub-height and power signals.


## Publications:
MPC:
S. Boersma, B.M. Doekemeijer, S. Siniscalchi-Minna and J.W. van Wingerden
A constrained wind farm controller providing secondary frequency regulation: an LES study
Renewable Energy, vol. x(y), pp. z-zz, 2018 (in press)

SMPC:
S. Boersma, T. Keviczky and J.W. van Wingerden
Stochastic Model Predictive Control: uncertainty impact on wind farm power tracking
American Control Conference, 2019 (under review)

WRC:
S. Raach, S. Boersma, B.M. Doekemeijer, J.W. van Wingerden and P.W. Cheng
Lidar-based closed-loop wake redirection in high-fidelity simulation
Journal of Physics: Conference Series, 2018


## Notes:

1) The mexcdf toolbox is only used to extract PALM data from an .nc file. This occurs only at the end of a simulation. The toolbox is reduced to its minimum. If problems occur, please download a full version of the mexcdf toolbox and place it in the correct directory.

2) PALM restart data is not provided due to its size. This data should be placed in the RESTART_DATA directory.

3) A solver (like CPLEX) and Yalmip are not provided due to their size. They should be installed and placed in the libraries directory of each controller.
 
4) The folder structure between PALM revisions can slightly differ. Minor modification should be made in the controller .m files when a different PALM folder structure is considered. 
