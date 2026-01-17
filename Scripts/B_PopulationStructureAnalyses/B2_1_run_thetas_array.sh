#!/bin/bash
#SBATCH --account=mcnew
#SBATCH --job-name=thetas
#SBATCH --partition=standard
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/output.thetas.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem=100gb

# must be run on puma

source /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh

OUTDIR=/xdisk/mcnew/finches/dannyjackson/finches
WORKDIR=${OUTDIR}/analyses/effective_pop

SAF_DIR=${WORKDIR}/safs
SFS_DIR=${WORKDIR}/sfs_files
THETA_DIR=${WORKDIR}/thetas

mkdir -p slurm_output "$SAF_DIR" "$SFS_DIR" "$THETA_DIR"
cd "$WORKDIR"

POPS=(crapre crapost forpre forpost parpre parpost)
POP=${POPS[$SLURM_ARRAY_TASK_ID]}

THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "POP=$POP  Threads=$THREADS"

~/programs/angsd/angsd \
  -bam ${OUTDIR}/referencelists/${POP}bams.txt \
  -rf ${OUTDIR}/referencelists/autosomes.txt \
  -out ${SAF_DIR}/${POP} \
  -dosaf 1 \
  -GL 1 \
  -doGlf 2 \
  -doMajorMinor 1 \
  -doCounts 1 \
  -doDepth 1 \
  -setMinDepthInd ${MINDEPTHIND} \
  -minInd ${MININD} \
  -minQ ${MINQ} \
  -minMapQ ${MINMAPQ} \
  -anc ${REF} \
  -nThreads ${THREADS}

#!/bin/bash
#SBATCH --account=mcnew
#SBATCH --job-name=thetas
#SBATCH --partition=standard
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/output.thetas.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem=50gb


source /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh

OUTDIR=/xdisk/mcnew/finches/dannyjackson/finches
WORKDIR=${OUTDIR}/analyses/effective_pop

SAF_DIR=${WORKDIR}/safs
SFS_DIR=${WORKDIR}/sfs_files
THETA_DIR=${WORKDIR}/thetas

mkdir -p slurm_output "$SAF_DIR" "$SFS_DIR" "$THETA_DIR"
cd "$WORKDIR"

POPS=(crapre crapost forpre forpost parpre parpost)
POP=${POPS[$SLURM_ARRAY_TASK_ID]}

THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "POP=$POP  Threads=$THREADS"

~/programs/angsd/misc/realSFS -P ${THREADS} \
  ${SAF_DIR}/${POP}.saf.idx > ${SFS_DIR}/${POP}.sfs

~/programs/angsd/misc/realSFS saf2theta \
  ${SAF_DIR}/${POP}.saf.idx \
  -outname ${THETA_DIR}/${POP} \
  -sfs ${SFS_DIR}/${POP}.sfs

~/programs/angsd/misc/thetaStat do_stat \
  ${THETA_DIR}/${POP}.thetas.idx \
  > ${THETA_DIR}/${POP}.thetas.idx.stats
