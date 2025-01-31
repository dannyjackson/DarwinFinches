# bam statistics 


## Compute statistics on bam files 
# some alignment stats 
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/

stating () { 
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".sorted.marked.clipped.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".depthstats

samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".sorted.marked.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt  


sbatch stating_dupmarkedbams.sh 
Submitted batch job 1934432
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/

stating () { 

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".sorted.marked.clipped.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt  


sbatch stating_clipoverlap.sh 
Submitted batch job 9591145



cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/
ls *all* | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt

#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/

stating () { 
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.sorted.marked.clipped.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.depthstats

samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt  





sbatch stating_all_dupmarkedbams.sh 
Submitted batch job 9544175


#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/

stating () { 

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt  

sbatch stating_all_dupmarkedbams_fixerror.sh 
Submitted batch job 1934522


# compute average and sd depth per individual # sorted marked bam files
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/


while read -r finch;
do 
  echo $finch
  echo $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

while read -r finch;
do 
  echo "all" $finch
  echo "all" $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".all.depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt 

sbatch averagedepth_sortedbams.sh 
Submitted batch job 9566373



# compute average and sd depth per individual # sorted marked bam files
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/


while read -r finch;
do 
  echo $finch
  echo $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

while read -r finch;
do 
  echo "all" $finch
  echo "all" $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".all.depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt 



while read -r finch;
do 
  samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$finch".sorted.marked.clipped.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".depthstats
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 




# Compute per-sample sequencing depth for bases with mapping quality higher than 20


module load samtools
while read -r finch;
do 
  samtools depth -q 30 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.depthstats
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt

ls *bam | awk 'BEGIN {FS = "."} {print $1}' | wc -l > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_indelrealigned.txt



### test for effects of duplicate sequencing effort










# python script for summarizing depth stats from ANGSD
python
import numpy as np
import pandas as pd 

df = pd.read_csv('stats.NC_044571.depthSample', sep='\t', header=None)
cols = df.shape[1]
multiplier = list(range(1, cols+1))

df_adj = df.mul(multiplier, axis = 1)

counts = df.sum(axis = 1)
totaldepth = df_adj.sum(axis = 1)
avgdepth = totaldepth / counts

df_names = pd.read_csv('/xdisk/mcnew/dannyjackson/finches/reference_lists/sample_species_treatment.txt', sep='\t')

df_final = pd.concat([df_names, avgdepth], axis=1)

df_final = df_final.rename({0: 'avgdepth'}, axis=1)

df_final.groupby(['species', 'treatment']).mean()

df_final.sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'CRA'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'FOR'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'PAR'].sort_values(by=['avgdepth']).round(decimals=2)

print(df_final.groupby(['species', 'treatment']).size())

df_final.groupby(['species', 'treatment']).mean('MEAN_DEPTH')





CRA     post        6.712184
        pre        12.027188
FOR     post        5.676069
        pre         7.770516
PAR     post        6.038819
        pre        10.395772