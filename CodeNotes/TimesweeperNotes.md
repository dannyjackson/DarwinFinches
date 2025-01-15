# Timesweeper notes

https://github.com/SchriderLab/Timesweeper/blob/master/README.md#workflow-overview

1. Create a SLiM script based on the example_demo_model.slim example

cd /xdisk/mcnew/dannyjackson/finches/timesweeper

conda create -n blinx -c conda-forge -c bioconda python slim=3.7 samtools bcftools
conda activate blinx

git clone https://github.com/SchriderLab/Timesweeper.git
cd Timesweeper
pip install .

2. Simulate demographic model with time-series sampling sim_custom if using custom SLiM script
    
    Note: If available, we suggest using a job submission platform such as SLURM to parallelize simulations. This is the most resource and time-intensive part of the module by far.
    
    Optional: Preprocess VCFs simulated without timesweepers simulation modules by merging with process_vcfs

3. Create features for the neural network with condense

4. Train networks with train

5. Run detect on VCF of interest using trained models and input data