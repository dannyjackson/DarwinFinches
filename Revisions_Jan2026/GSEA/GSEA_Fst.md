# D1. Gene Set Enrichment Analysis

## Overview

This stage tests whether genes associated with particular biological pathways or Gene Ontology (GO) terms tend to occur in genomic regions showing elevated signatures of selection.

For this first implementation, I'm using the selection statistic is **FST**. For each gene, I am taking the average statistic across all windows that overlap the gene. This means that every gene receives one FST score regardless of its length or the number of windows it overlaps.

The same workflow will later be applied to:

- ΔTajima's D
- the composite statistic
- RAiSD μ

---

# Directory structure

Create a directory for Stage D:

```bash
mkdir -p gene_set_enrichment/{genes,windows,gene_scores,gsea}
```

---

# Step 1. Convert the genome annotation to BED

The gene annotation should be the annotation corresponding to the reference genome used for the population genomic analyses.

The output BED contains:

```text
chrom
start
end
gene_id
gene_name
strand
```

BED uses 0-based, half-open coordinates.

Run:



```bash
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/gene_set_enrichment

GFF="/xdisk/mcnew/finches/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file

python3 Scripts/D1_make_gene_bed.py \
    $GFF \
    genes/genes.bed
```

## `Scripts/D1_make_gene_bed.py`

```python
#!/usr/bin/env python3

import argparse


def parse_attributes(attribute_string):
    """Convert GFF3 attribute field into a dictionary."""
    attributes = {}

    for item in attribute_string.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value

    return attributes


parser = argparse.ArgumentParser(
    description="Extract gene coordinates from a GFF3 file and write BED."
)

parser.add_argument(
    "gff",
    help="Input GFF3 annotation"
)

parser.add_argument(
    "output",
    help="Output BED file"
)

args = parser.parse_args()


with open(args.gff) as infile, open(args.output, "w") as outfile:

    for line in infile:

        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")

        if len(fields) != 9:
            continue

        chrom = fields[0]
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]

        if feature != "gene":
            continue

        attributes = parse_attributes(fields[8])

        gene_id = attributes.get("ID", "NA")
        gene_name = attributes.get(
            "Name",
            attributes.get("gene", gene_id)
        )

        # Remove common GFF3 prefix if present
        if gene_id.startswith("gene:"):
            gene_id = gene_id.replace("gene:", "", 1)

        # GFF3 is 1-based inclusive.
        # BED is 0-based, half-open.
        bed_start = start - 1
        bed_end = end

        outfile.write(
            f"{chrom}\t"
            f"{bed_start}\t"
            f"{bed_end}\t"
            f"{gene_id}\t"
            f"{gene_name}\t"
            f"{strand}\n"
        )
```


---

# Step 2. Convert FST windows to a common format

All selection statistics used in Stage D will ultimately be converted into the following format:

```text
chrom   start   end   score
```

No header is used because this file will be passed directly to `bedtools`.

For the current FST analysis, the Stage C output contains a chromosome, window midpoint, and FST value.

For example:

```text
chromo        position       Nsites       fst
chr1       25000        200          0.014
chr1       37500        217          0.021
...
```

Convert this into BED intervals using the known FST window size.

For 50-kb windows:

```bash
FST_CRA=/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/50000/crapre_crapost.50000.fst.depthmapfiltered

FST_FOR=/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/forpre_forpost/50000/forpre_forpost.50000.fst.depthmapfiltered

FST_PAR=/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/parpre_parpost/50000/parpre_parpost.50000.fst.depthmapfiltered

python Scripts/D1_standardize_windows.py \
    --input $FST_CRA \
    --output windows/cra.fst.windows.bed \
    --chrom-col chromo \
    --position-col position \
    --score-col fst \
    --window-size 50000

python Scripts/D1_standardize_windows.py \
    --input $FST_FOR \
    --output windows/for.fst.windows.bed \
    --chrom-col chromo \
    --position-col position \
    --score-col fst \
    --window-size 50000

python Scripts/D1_standardize_windows.py \
    --input $FST_PAR \
    --output windows/par.fst.windows.bed \
    --chrom-col chromo \
    --position-col position \
    --score-col fst \
    --window-size 50000
```

The script can also accept statistics that already contain explicit window start and end positions. This makes it reusable later for ΔTajima's D, the composite statistic, and RAiSD.


---

# Step 3. Identify every FST window overlapping each gene

Use `bedtools intersect`.

For CRA:

```bash
bedtools intersect \
    -a genes/genes.bed \
    -b windows/cra.fst.windows.bed \
    -wa -wb \
    > gene_scores/cra.fst.gene_window_overlap.tsv

bedtools intersect \
    -a genes/genes.bed \
    -b windows/for.fst.windows.bed \
    -wa -wb \
    > gene_scores/for.fst.gene_window_overlap.tsv

bedtools intersect \
    -a genes/genes.bed \
    -b windows/par.fst.windows.bed \
    -wa -wb \
    > gene_scores/par.fst.gene_window_overlap.tsv
```

The resulting file contains:

```text
gene chromosome
gene start
gene end
gene ID
gene name
strand
window chromosome
window start
window end
FST
```

Example:

```text
chr1  100000  140000  geneA  ABC1  +  chr1  87500   137500  0.04
chr1  100000  140000  geneA  ABC1  +  chr1  100000  150000  0.06
chr1  100000  140000  geneA  ABC1  +  chr1  112500  162500  0.05
```

---

# Step 4. Calculate mean FST for each gene

Run:

```bash
python3 Scripts/D1_mean_gene_scores.py \
    gene_scores/cra.fst.gene_window_overlap.tsv \
    gene_scores/cra.fst.gene_scores.tsv

python3 Scripts/D1_mean_gene_scores.py \
    gene_scores/for.fst.gene_window_overlap.tsv \
    gene_scores/for.fst.gene_scores.tsv

python3 Scripts/D1_mean_gene_scores.py \
    gene_scores/par.fst.gene_window_overlap.tsv \
    gene_scores/par.fst.gene_scores.tsv
```

The resulting file looks like:

```text
gene_id gene_name chrom gene_start gene_end mean_score n_windows
geneA   ABC1      chr1  100000     140000   0.050      3
geneB   DEF2      chr1  180000     195000   0.041      2
geneC   GHI3      chr1  250000     310000   0.037      5
```

Importantly, **every gene appears only once** in the GSEA ranking.

---

# Step 5. Prepare gene sets

`fgsea` needs a mapping between pathway IDs and gene IDs.

For GO enrichment, create:

```text
gene_set_enrichment/gene_sets/go_terms.tsv
```

With columns:

```text
term_id     gene_id     term_name                       ontology
GO:0002376  geneA       immune_system_process           BP
GO:0002376  geneB       immune_system_process           BP
GO:0006955  geneB       immune_response                 BP
```

The **gene IDs must match the IDs in `genes.bed`**.

No downstream script needs to change as long as the file contains:

```text
term_id
gene_id
```

---

# Step 6. Run pre-ranked GSEA

For FST, larger values are the values of interest.

By using option "pos", genes are ranked from:

```text
highest mean FST
        ↓
lowest mean FST
```

## Download GO Term
```
mkdir -p gene_sets
cd gene_sets

wget https://current.geneontology.org/annotations/gaf/CHICK-uniprot.gaf.gz
gunzip -c CHICK-uniprot.gaf.gz > CHICK-uniprot.gaf

grep -v '^!' CHICK-uniprot.gaf \
    | awk -F '\t' 'BEGIN{OFS="\t"} {print $3,$5,$9}' \
    | sort -u \
    > chicken_symbol_GO.tsv

cd ..

tail -n +2 gene_scores/cra.fst.gene_scores.tsv \
    | cut -f2 \
    | sort -u \
    > /tmp/cra_genes.txt

cut -f1 gene_sets/chicken_symbol_GO.tsv \
    | sort -u \
    > /tmp/chicken_go_symbols.txt

echo "Finch genes:"
wc -l /tmp/cra_genes.txt

echo "Chicken GO symbols:"
wc -l /tmp/chicken_go_symbols.txt

echo "Exact symbol matches:"
comm -12 /tmp/cra_genes.txt /tmp/chicken_go_symbols.txt | wc -l

python3 Scripts/D1_map_GO_finch_chicken.py
```
Run the fgsea.

The final argument specifies the type of enrichment:

```text
pos = enrichment toward high statistic values
std = enrichment at either end of a signed statistic
neg = enrichment toward low statistic values
```

Command:

```bash
Rscript Scripts/D1_run_fgsea.R \
    gene_scores/cra.fst.gene_scores.tsv \
    gene_sets/go_terms.tsv \
    gsea/cra.fst \
    pos

Rscript Scripts/D1_run_fgsea.R \
    gene_scores/for.fst.gene_scores.tsv \
    gene_sets/go_terms.tsv \
    gsea/for.fst \
    pos

Rscript Scripts/D1_run_fgsea.R \
    gene_scores/par.fst.gene_scores.tsv \
    gene_sets/go_terms.tsv \
    gsea/par.fst \
    pos
```
