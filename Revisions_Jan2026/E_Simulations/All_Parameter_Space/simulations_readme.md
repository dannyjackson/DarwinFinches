# Simulation and analysis workflow

## 1. Parameter space definition
To model the detectability of a selection region, we need to define a few parameters:
1. Frequency of the allele under selection at time point 1
2. Effective population size
3. Potential population decline rate
4. Strength of selection
5. Genomic segment length

None of these parameters have a single value that applies to our system. For some, this is intuitive: the starting allele frequency and the strength of selection both vary across the genome.

### These parameters require some additional attention:
1. **Effective population size**: We've estimated the Ne of these samples, which are ~100-200k. However, our sampling scheme is location specific. The 100-200k Ne defines the standing genetic variation, but the Ne of all contributing individuals to our T1 and T2 samples is much smaller and this determines detectibility.
2. **Population decline rate**: If populations are declining, this can look like selection as well. We do not know if this is happening. Ne estimations at this short of a time frame are very limited.
3. **Genomic segment length**: I defined this as 10k because that was the window size of our selection scans, and due to computational limitations. However, modeling genomic segments of only this size does not allow for the influence of linkage.

### Therefore, I modeled a range of possible parameters.

I defined an array of scenarios spanning five starting allele frequencies (0.05–0.95), four effective population sizes (10,000–150,000), five population decline rates (1.00–0.80), three temporal sampling intervals (5–15 generations), five selection coefficients (s = 0.01-0.5), and a fixed genomic segment length (10 kb). The genomic segment length could be varied in a future iteration.

## 2. Forward-time simulations with selection and demography
100 replicate forward-time simulations were run using SLiM for each combination of the above parameters. Simulations modeled a diploid population experiencing selection on standing variation, optional population decline, and two sampling time points separated by a specified number of generations. I sampled 8 individuals at each time point.

## 3. Tree-sequence recording, temporal sampling, recapitation of genealogies, and neutral mutation overlay
All simulations used tree-sequence recording. At each replicate, individuals were sampled at two time points (T1 and T2), and their pedigree identifiers were recorded to enable downstream mapping of samples to genealogies. For each replicate, the forward-time tree sequence was recapitated using msprime under a Wright–Fisher coalescent model to complete the ancestral history beyond the forward simulation. Neutral mutations were overlaid onto the recapitated genealogies at a fixed mutation rate, generating sequence-level variation consistent with the simulated demographic history.

This is a much more computationally efficient way of running simulations rather than the more conceptually straightforward way of tracking every mutation. I ran a subset of Ne=10k simulations using the straightforward method and confirmed that it produced similar results.

## 4. Computation of summary statistics
Using the temporally sampled individuals, I calculated nucleotide diversity (π) at each time point, divergence between time points (dxy), Hudson’s FST between T1 and T2 samples, Tajima’s D at each time point, and the change in Tajima’s D across time (ΔTajima’s D).

## 5. Visualization and discrimination of temporal genomic patterns.
For each combination of effective population size and sampling interval, replicate outcomes were summarized in faceted scatterplots of Hudson’s FST versus ΔTajima’s D (facets by starting allele frequency and decline rate), with points colored by selection strength. I computed 95% multivariate normal ellipses depicting replicate dispersion, and within-facet AUC values quantifying the separation of selected versus neutral simulations. Plots used standardized axis limits for comparability and were exported as pngs.