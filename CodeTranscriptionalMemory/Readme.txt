Included are the data and scripts to re-create the computational analysis figures for the manuscript:

"Gene-specific transcriptional memory in mammalian cell lineages"
Nicholas E. Phillips, Aleksandra Mandic, Saeed Omidi, Felix Naef, David M. Suter

All analysis was performed using MATLAB R2017b

The analysis method using Gaussian processes is found in the "Data generator" folder within Figure 4. Note that "Data generator" contains the code for generating MCMC posterior estimates for all genes/conditions considered in the "main.m" file. The MCMC runs are then used to create Figure 4, and additionally the mother-daughter analysis used for Figure 5 and panel f of Figure 6. The file to re-run the analysis is found as main.m in Figure 4/ Data generator.

The analysis using the "main.m" was run on a server, and 10,000 MCMC runs takes approx 10
days. The number of MCMC runs can be adjusted in each of the individual scripts in the "Data generator" folder.

Figures 4,5 and 6 can then be recreated, and the expected output
corresponds to the figures in the paper.

Please address questions to: Nicholas.phillips@epfl.ch