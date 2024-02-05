# Nonpareil-LR
metagenomic coverage estimation for long and accurate reads

## Background
Metagenomic long reads sequencing is now more and more popular. Estimating the coverage of a metagenomic sample provides useful information about the diversity of the sample. Lander-Waterman expectation, or the generalized Stevens’ theorem can be used to do such an estimation. For shotgun metagenomic reads (short reads), nonpareil has implemented this idea. However, for long and accurate long reads, the idea behind nonpareil(e.g., ungapped alignment) fails due to much longer sequneces and lower identities. 
