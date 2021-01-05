# DppTrafficking-
This repository includes code used for the parameterization procedure in Romanova-Michaelidis et al.

The scripts are written in C++ and are used to numerically solve Eq. S1-S5 provided in the Supplementary Section 1 attached to Romanova-Michaelides et al (link), under the conditions of the FRAP experiment. The parameter values used to perform the numerical calculations are chosen in ranges that are defined in the Supplementary Section 2 attached to Romanova-Michaelides et al (link).

The scripts can be compiled in xCode or similar interface or by using the terminal window.

The output of each run is a file that includes 300 rows that summarize specific sets of parameters together with the corresponding quality of fit to the experimental FRAP data. 

For the fitting routine described in the manuscript, we use HPC to run the scripts repeatedly until 3 10^7 parameter sets have been explored. The output files are then used to obtain the parameter sets that provde sattisfactory fits to the data (see Supplementary Section 2 attached to Romanova-Michaelides et al, link). 
