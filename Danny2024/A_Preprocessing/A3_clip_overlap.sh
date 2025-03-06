## Clip overlapping read pairs
#!/bin/bash

#SBATCH --job-name=clipoverlap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.clipping.%j

module load parallel

clipping () {
echo "clipping" "$@" >> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/clippingstats.txt 

/home/u15/dannyjackson/programs/bamUtil-master/bin/bam clipOverlap --in /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam --out /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.sorted.marked.clipped.bam --stats --params


echo "done " $@ >>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/clippingstats.txt 
}


export -f clipping 

parallel -j 12 clipping :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt

sbatch clipping.all.sh 
Submitted batch job 1946694