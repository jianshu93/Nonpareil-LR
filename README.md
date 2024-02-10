
# NonpareilX
metagenomic coverage estimation for long and accurate reads

## Background
Metagenomic long reads sequencing is now more and more popular. Estimating the coverage of a metagenomic sample provides useful information about the diversity of the sample. Lander-Waterman expectation, or the generalized Stevens’ theorem can be used to do such an estimation. For shotgun metagenomic reads (short reads), nonpareil has implemented this idea. However, for long and accurate reads, the idea behind nonpareil(e.g., ungapped alignment) fails due to much longer sequneces and lower identities. And more importantly, sampling 10% of the metagenome to estimate coverage is not ideal solution. All versus all comparisons should be performed, despite being more computationally expensive.

## New algorithm for calculating the number of non-redundant reads
The key step is to obtain the number of non-redundant reads via all versus all comparison of each sequence in the metagenome, which is O(N^2), impractical for real-world metagenomic samples. To solve the problem, we adopt cutting-edge sketching algorithms and graph based nearest neighbor search for finding neighbor sequences of each sequence in metagenome. There are three steps involved:
1. (Order) MinHash & HNSW for extracting the most similar sequences in an approximate manner for distantly related sequences in metagenome to each sequence. Order MinHash is kmer-based sketching method, approximating edit distance quite well. Via HNSW, the running time respect to number of sequences in metagenome is O(N*log(N)), see GSearch paper for details. However, this method is not optimal for sequences with different length.
2. To have even accurate estimation of alignment-based identity for sequneces with different length to find the most similar sequneces, we use co-linear chaining of anchors (can ber minimizers, MUMs or MEMs) with overlaps and gap costs, which approximate edit distance well above 90% identity (semi-global or global alignmend mode).
3. Exact semi-global alignment via vsearch implementation (not unit score scheme as in Edlib) only for sequences with a edit score for above 94% sequence identity, obtained above via chaining.

![Alt!](https://github.com/jianshu93/NonpareilX/blob/master/orderminhash.jpg?raw=true)

## Lander-Waterman expectation
We will then follow the same equation for estimation of coverage after obtaining the number of non-redundant reads in above steps.