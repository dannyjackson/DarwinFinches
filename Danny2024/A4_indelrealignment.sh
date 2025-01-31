## Indel realignment

#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j

module load samtools

while read -r finch;
do 
  samtools index /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam

  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
    -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam \
    -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch".intervals

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 


# for all files 
#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j

module load samtools

while read -r finch;
do 
  samtools index /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam

  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
    -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam \
    -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch".intervals

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

#!/bin/bash

#SBATCH --job-name=indelrealignment
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelrealignment.%j

while read -r finch;
do 
  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
  -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  --consensusDeterminationModel USE_READS \
  -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam \
  --targetIntervals /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch".intervals \
  -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

sbatch indelrealignment.sh 
Submitted batch job 1942426

# for all files 
#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j

module load samtools

while read -r finch;
do 
  samtools index /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".all.sorted.marked.clipped.bam

  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
    -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".all.sorted.marked.clipped.bam \
    -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch"_all.intervals

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt 

sbatch ~/programs/slurmscripts/indelrealignment.all.sh

#!/bin/bash

#SBATCH --job-name=indelrealignment
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelrealignment.%j

while read -r finch;
do 
  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
  -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  --consensusDeterminationModel USE_READS \
  -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".all.sorted.marked.clipped.bam \
  --targetIntervals /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch"_.intervals \
  -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch"_all.realigned.bam

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt 

sbatch ~/programs/slurmscripts/indelmap.all.sh 
Submitted batch job 9724713

sbatch indelrealignment.all.sh 
Submitted batch job 9750831