
sbatch --account=mcnew \
        --job-name=TENM3_CRA \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.TENM3_CRA.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh TENM3 cra


# Gene lists
# All three species:
## TENM3

sbatch --account=mcnew \
        --job-name=TENM3_CRA \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.TENM3_CRA.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh TENM3 cra


sbatch --account=mcnew \
        --job-name=TENM3_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.TENM3_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh TENM3 for


sbatch --account=mcnew \
        --job-name=TENM3_PAR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.TENM3_PAR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh TENM3 par

## NUDT4
sbatch --account=mcnew \
        --job-name=NUDT4_CRA \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.NUDT4_CRA.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh NUDT4 cra


sbatch --account=mcnew \
        --job-name=NUDT4_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.NUDT4_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh NUDT4 for


sbatch --account=mcnew \
        --job-name=NUDT4_PAR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.NUDT4_PAR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh NUDT4 par


## BMPER
sbatch --account=mcnew \
        --job-name=BMPER_CRA \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.BMPER_CRA.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh BMPER cra


sbatch --account=mcnew \
        --job-name=BMPER_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.BMPER_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh BMPER for


sbatch --account=mcnew \
        --job-name=BMPER_PAR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.BMPER_PAR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.sh BMPER par


# CRA
## PODXL 
sbatch --account=mcnew \
        --job-name=PODXL_CRA \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.PODXL_CRA.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh PODXL cra

## ANGPT1
sbatch --account=mcnew \
        --job-name=ANGPT1_CRA \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.ANGPT1_CRA.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh ANGPT1 cra

## HPSE
sbatch --account=mcnew \
        --job-name=HPSE_CRA \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.HPSE_CRA.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh HPSE cra

## ITGA2B
sbatch --account=mcnew \
        --job-name=ITGA2B_CRA \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.ITGA2B_CRA.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh ITGA2B cra

# FOR
## PRTFDC1
sbatch --account=mcnew \
        --job-name=PRTFDC1_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.PRTFDC1_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh PRTFDC1 for

## FARS2
sbatch --account=mcnew \
        --job-name=FARS2_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.FARS2_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh FARS2 for

## MCM9
sbatch --account=mcnew \
        --job-name=MCM9_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.MCM9_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh MCM9 for

## CELF1
sbatch --account=mcnew \
        --job-name=CELF1_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.CELF1_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh CELF1 for

## ADAT1
sbatch --account=mcnew \
        --job-name=ADAT1_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.ADAT1_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh ADAT1 for

## CHST4 
sbatch --account=mcnew \
        --job-name=CHST4_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.CHST4_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh CHST4 for

## TERF2IP
sbatch --account=mcnew \
        --job-name=TERF2IP_FOR \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.TERF2IP_FOR.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=0:10:00 \
        --mem=50gb \
        gene_pca_win.0.1.sh TERF2IP for




cp */*/*treatment.pdf treatment_plots/
mv treatment_plots/*TENM* treatment_plots/parallel_genes
mv treatment_plots/*NUDT* treatment_plots/parallel_genes
mv treatment_plots/*BMPER* treatment_plots/parallel_genes
mv treatment_plots/*for* treatment_plots/for_genes
mv treatment_plots/*cra* treatment_plots/cra_genes



# redo with 
for SP in cra for par; do
  case $SP in
    cra)
      GENES=("ANGPT1" "BMPER" "HPSE" "ITGA2B" "NUDT4" "PODXL" "TENM3")
      ;;
    for)
      GENES=("ADAT1" "BMPER" "CELF1" "CHST4" "FARS2" "MCM9" "NUDT4" "PRTFDC1" "TENM3" "TERF2IP")
      ;;
    par)
      GENES=("BMPER" "NUDT4" "TENM3")
      ;;
  esac

  for GENE in "${GENES[@]}"; do
    cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/pca/indv_win/${SP}/${GENE} || exit 1
    # Your commands here, e.g.
    echo "Processing $GENE for species $SP"
    while read -r chrom start end; do
    region="${chrom}:${start}-${end}"
    echo "Processing $region"

     # Run R script for plotting
     Rscript /xdisk/mcnew/finches/dannyjackson/finches/analyses/pca/indv_win/gene_pca_win_polygon.R ${GENE}_${SP}_${start}_${end} ${SP}

     done < "${GENE}.merged.windows.bed"
  done
done


cp */*/*poly_treatment.pdf treatment_plots_poly/
mv treatment_plots_poly/*TENM* treatment_plots_poly/parallel_genes
mv treatment_plots_poly/*NUDT* treatment_plots_poly/parallel_genes
mv treatment_plots_poly/*BMPER* treatment_plots_poly/parallel_genes
mv treatment_plots_poly/*for* treatment_plots_poly/for_genes
mv treatment_plots_poly/*cra* treatment_plots_poly/cra_genes