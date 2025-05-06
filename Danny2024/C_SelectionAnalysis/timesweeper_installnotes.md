# Timesweeper.md

# All phased vcfs are in /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2 and vcf3


# source ~/.bashrc

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper

interactive -a mcnew -t 1:00:00 -g

module load python/3.8/3.8.12
# python3 -m venv --system-site-packages timesweeper_env

source timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 


# Make yaml file
#General
work dir: /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/test_output
slimfile: constant_population.slim
slim path: slim
experiment name: Reps_Size_20k

mut types: [2]

scenarios: ["neut", "ssv", "sdn"]

win_size: 51

num_sample_points : 2

inds_per_tp : 10  # Diploid inds

physLen : 5000000

sample sizes: [10, 10]

ploidy: 2

reps: 10

/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/timesweeper sim_custom -y example_config.yaml 



# Make Training Data (condense)

timesweeper condense -o test.pkl -m 0.02 -y  example_config.yaml 

# train neural networks
timesweeper train -i test.pkl -y example_config.yaml 
