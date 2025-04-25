# params_fst
source /home/u15/dannyjackson/programs/DarwinFinches/params_base.sh

ANGSD=~/programs/angsd/ # path to directory with angsd executables

CUTOFF=0.01

# define two colors to be used (alternating chromosomes in manhattan plots)
COLOR1="#A6C965"
COLOR2="#203000"

# define the names of the two populations that will be compared
POP=parpre

# source the setup file for fst
source ${SCRIPTDIR}/Genomics-Main/C_SelectionAnalysis/tajima/setup_tajima.sh

