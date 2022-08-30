This folder contains the Matlab scripts and functions, as well the necessary data files, to implement the parameter
estimation methods in Gubbins (2022) "Quantifying the relationship between within-host dynamics and transmission
for viral diseases of livestock"

MATLAB REQUIREMENTS
The scripts/functions were run using Matlab version 2020b and require the Statistics and Machine Learning and
Parallel Computing toolboxes. However, they can be easily adapted to run without the Parallel Computing toolbox
by changing any "parfor" loops to "for" loops in the functions with names beginning "ParEst".

PARAMETER ESTIMATION
For each [virus] (i.e. FMDV or SwIV):
ParEst_[virus].m - loads the data, implements the adaptive Metropolis scheme for each model parameterisation
                   and computes the DIC
Lhood_[virus].m - computes the log likelihood and prior for the input parameters for the specified model
                  parameterisation

The data are supplied in Excel format as supporting information with the paper: DataS1.xlsx (FMDV) and DataS2.xlsx
(SwIV) and as Matlab data files (here). These files contain (cf. Data S1 and Data S2):

[virus]_TransmissionExperimentData.mat
tStart - days post challenge of the donor at which each challenge started 
tStop - days post challenge of the donor at which each challenge ended 
D - challenge outcomes (0=transmission did not occur; 1=transmission occurred)
tClin - days post challenge of the donor at which it first showed clinical signs (FMDV only)

[virus]_ViralTitreData.mat
tObs - days post challenge of the donor at which samples were taken
VI_B - log10 viral titre in blood (FMDV)
VI_NF - log10 viral titre in nasal fluid (FMDV)
VI_OPF - log10 viral titre in oesophageal-pharyngeal fluid (FMDV)
VI - log10 viral titre in nasal swabs (SwIV)

COMPUTING SUMMARY TRANSMISSION MEASURES
For each [virus]:
computeSummaryMeasures_[virus].m - loads the MCMC samples and computes the reproduction number and generation time
                                   for each animal using the exact formula and the heuristic approximation; for FMDV
                                   it also computes theta for each animal

