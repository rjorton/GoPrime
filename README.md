# GoPrime
GoPrime is a tool to predict the performance of real time PCR primers and probes. GoPrime takes multiple sequence data and a primer/probe set, identiifes all possible binding sites of indivudal primers and probes, combines these into potential primer-probe-primer sets and calculates a predicted CT value for the best set. GoPrime can handle with ambigious bases (degenerative primers), aligned or unaligned sequnece data, and outputs detailed information on primer/probe binding across the genome, as well as predicted CTs for all possible primer-probe-primer sets. 

GoPrime is idealy suited to evaluate the sensitivity and specificty of genotype specific real time PCR primers and probes for viral data sets.

GoPrime is combined with GoPrimeTree for visualising the CT results on a phylogenetic tree of the data, and GoPrimeAlign for identifying possible genotypic specific primer/probe sites.

