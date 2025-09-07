# Project Overview
This repository contains scripts and code used to analyze time-series medium-coverage whole-genome sequence data from three species of Darwin's finch, each sampled before and after the introduction of the avian vampire fly.

### Central Question
How does a novel parasite affect patterns of selection across species of Darwin's finches?
To what extent is evolution parallel across taxa?

These analyses aim to identify signatures of selection that have acted on each species in an urban context. 

### Abbreviations:
   - **CRA / cra**: vegetarian finch (*Platyspiza crassirostris*)
   - **FOR / for**: medium ground finch (*Geospiza fortis*)
   - **PAR / par**: small tree finch (*Camarhynchus parvulus*)

### Pipeline
My analyses use generalizable scripts from the repository [Genomics Main](https://github.com/dannyjackson/Genomics-Main), which is retained as a submodule within this repository in its state at the time of publication of these scripts.

To implement the pipeline, follow the markdown files in sequential order. Each markdown filename follows a structured format:
**[Stage Letter][Step Number]_[Description]**
   - The **letter** represents the pipeline stage.
   - The **number** indicates the order in which the scripts should be excecuted within that stage.
   - The **description** briefly summarizes the purpose of the script.
Many markdown files call various scripts within the Genomics Main repository, all labeled using a similar syntax.

## Pipeline Stages

A. **Preprocessing**  
We follow the guidance for mitigating batch effects in low-coverage genomic sequence data from [Lou and Therkildsen 2022](https://doi.org/10.1111/1755-0998.13559).
   - Trims sequence reads
   - Aligns and filters the genome
   - Clips overlapping read pairs
   - Indel realignment
   - Removes individuals with low depth of coverage (<3x)
   - Generates a list of SNPs (MAF filter only applies to *F<sub>ST</sub>*, not *Tajima's D*, analyses)
   - Methods
      - [Trimmomatic](https://github.com/timflutre/trimmomatic)
      - [bwa](https://github.com/lh3/bwa)
      - [BamUtil clipOverlap](https://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap)
      - [gatk](https://gatk.broadinstitute.org/hc/en-us) (requires version 3.7-0 or earlier)
      - [samtools](https://www.htslib.org/)
      - [angsd](https://www.popgen.dk/angsd/)


B. **Population Structure Analysis**  
   - Identifies species and population clusters
   - Also identifies related individuals for filtering from downstream analysis
   - Methods
      - [PCAngsd](https://github.com/Rosemeis/pcangsd/)
         - *Principal component analysis (PCA)*
         - *Admixture modeling*

C. **Selection Analysis**  
   - Detects signatures of natural selection through:
     - **Comparisons of pre- vs post invasion samples** 
       - Implemented in [angsd](https://www.popgen.dk/angsd/)
       - *F<sub>ST</sub>*
       - *∆Tajima's D*
       - Composite statistic between  *F<sub>ST</sub>* and *∆Tajima's D*.
  
D. **Gene Investigations**  
   - Diving deep into the context of interesting genes identified in Stage C.
     


### Acknowledgements
This project is a central project of my postdoctoral fellowship with **Sabrina McNew** (University of Arizona). It is a collaboration between the two of us, a co-mentored Master's student **Logan Vossler** (University of Arizona), **Leonardo Campagna** (Cornell University), **Sangeet Lamichhaney** (Kent State University), **Jessica Rick** (University of Arizona), and **Birgit** Fessl (4.	Charles Darwin Research Station).


All sequence data will be available on the **Sequence Read Archive** upon publicaton, under BioProject ID PRJNA1320981.

