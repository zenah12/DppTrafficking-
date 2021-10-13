# DppTrafficking-
This repository includes code used for the parameterization procedure in Romanova-Michaelidis et al. The subfolders within contain the experimental data files fitted using the different procedures described in the main article and the supplementary material (Section 2.5.2).

*************
Nanobody_Fits
*************
This folder contains the experimental data for the nanobody experitment for the three conditions discussed in the main text (large discs, L=144um; small discs L=80um; pentagone mutant discs). The code to fit the data to the theoretical curve derived in Section x.y.z in the Supplementary material is provided in the Wolfram Mathematica notebooks contained in the folder (one notebook for each condition provided). The code read in the experimental data, defines the fit function and simultaneously fits the data to the theoretical function using a non linear model. This yields fitted values for A, B and p and their confidences as discussed in Section x.y.z in the Supplementary material. Those fitted values are used to obtain estimates for kr, kN and ko as discussed in the main text and Supplementary material Section x.y.z.

*****************
FRAP numerics
*****************
These scripts are written in C++ and are used to numerically solve Eq. S1-S5 provided in the Supplementary Section 1 attached to Romanova-Michaelides et al, under the conditions of the FRAP experiment. The parameter values used to perform the numerical calculations are chosen in ranges that are defined in the Supplementary Section 2 attached to Romanova-Michaelides et al.

The scripts can be compiled in xCode or similar interface or by using the terminal window.

The output of each run is a file that includes Tth (defined in Global.hpp) rows that summarize specific sets of parameters together with the corresponding quality of fit to the experimental FRAP data as calculated by the R^2.  

For the fitting routine described in the manuscript, we use HPC to run the scripts repeatedly until 3 10^7 parameter sets have been explored. The output files are then used to obtain the parameter sets that provde sattisfactory fits to the data (see below).


*********************
Posteriors something
*********************
