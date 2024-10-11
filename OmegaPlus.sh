# OmegaPlus
git clone https://github.com/alachins/omegaplus
cd omegaplus
ls
./compile_all.sh

./OmegaPlus -name Test \
    -input ./examples/fasta.FA \
    -ld RSQUARE \
    -minwin 10 \
    -maxwin 100\
    -grid 100

cd /xdisk/mcnew/dannyjackson/finches/sweed/OmegaPlus/for/pre

~/programs/omegaplus/OmegaPlus -name for_pre -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre/for_pre.vcf -grid 100 -length 100000 -minwin 100 -maxwin 10000 -seed 1500
