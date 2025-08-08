# Clip Overlapping Read Pairs (bamUtil `clipOverlap`)

**Goal:** For each *duplicate-marked, coordinate-sorted* BAM (`*.all.sorted.marked.bam`), clip bases from overlapping regions of read pairs to avoid double counting in downstream depth/allele calculations.

**Tool:** [`bamUtil clipOverlap`](https://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap) (invoked via the `bam` binary)

---

## SLURM job: clip all samples in parallel

```bash
#!/bin/bash
#
# Clip overlapping paired-end reads using bamUtil's clipOverlap in parallel.
# Assumes inputs live in .../sortedmarkedbamfiles/ and a list of sample IDs
# exists at .../samplenames_dupmarked.all.txt (one ID per line, no extensions).
#
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

# Run 12 samples in parallel (matches --ntasks=12)
parallel -j 12 clipping :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt
```
---

### Submit

```bash
sbatch clipping.all.sh
```

