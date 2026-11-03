#!/bin/bash
# Robustness check: rerun the per-species fgsea with GO annotations propagated to ancestor
# terms (is_a + part_of) and compare the rank of the focal terms with the direct-annotation
# results used in the manuscript. Nothing in the existing gsea/ or gene_sets/ directories is
# overwritten; all new files go to gene_sets/propagated/ and gsea_propagated/.
#
# Run from the GSEA project directory, the one that holds gene_scores/, gene_sets/ and gsea/
# (the layout in Revisions_Jan2026/GSEA/GSEA_Fst.md). Put the three scripts from this folder
# somewhere and set CHECKDIR to it.
#
#   export CHECKDIR=/path/to/go_propagation_check
#   bash $CHECKDIR/run_propagation_check.sh
#
# Requires: python3 with pandas and scipy (no other dependencies);
#           R with fgsea, data.table, GO.db (already used by D1_run_fgsea.R).
set -euo pipefail

CHECKDIR="${CHECKDIR:?set CHECKDIR to the folder holding propagate_gaf.py etc.}"
FGSEA_R="${FGSEA_R:-DarwinFinches/Revisions_Jan2026/GSEA/D1_run_fgsea.R}"   # path to the existing script
GO_RELEASE="${GO_RELEASE:-2026-07-26}"   # from the !go-version header line of CHICK-uniprot.gaf
GAF="${GAF:-gene_sets/CHICK-uniprot.gaf}"
OBO="gene_sets/propagated/go-basic.${GO_RELEASE}.obo"

for f in "$GAF" "$FGSEA_R" gene_scores/cra.fst.gene_scores.tsv gsea/cra.fst.gsea.tsv; do
    [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }
done
mkdir -p gene_sets/propagated gsea_propagated

# 1. Ontology from the same GO release as the annotation file
if [ ! -s "$OBO" ]; then
    wget -q -O "$OBO" "https://release.geneontology.org/${GO_RELEASE}/ontology/go-basic.obo" \
        || { echo "could not fetch go-basic.obo for release ${GO_RELEASE}; download it by hand to $OBO" >&2; exit 1; }
fi
echo "ontology: $(grep -m1 '^data-version' "$OBO")"

# 2. Propagate chicken annotations (BP only). Also write the direct table with the same NOT /
#    obsolete filters, so the comparison isolates the effect of propagation.
python3 "$CHECKDIR/propagate_gaf.py" --gaf "$GAF" --obo "$OBO" --aspect P \
    --out gene_sets/propagated/chicken_symbol_GO_BP.propagated.tsv \
    --direct-out gene_sets/propagated/chicken_symbol_GO_BP.direct_refiltered.tsv

# 3. Map to finch genes (same join as D1_map_GO_finch_chicken.py). The join uses the gene list
#    of one species; all three share the same gene annotation, as in the original pipeline.
python3 "$CHECKDIR/map_GO_finch_chicken.py" gene_scores/cra.fst.gene_scores.tsv \
    gene_sets/propagated/chicken_symbol_GO_BP.propagated.tsv gene_sets/propagated/go_terms.propagated.tsv

# 4. fgsea per species with the existing script, unchanged (minSize 10, maxSize 500, pos, eps 0, seed 1)
for SP in cra for par; do
    echo "== fgsea, propagated annotations: $SP"
    Rscript "$FGSEA_R" gene_scores/${SP}.fst.gene_scores.tsv \
        gene_sets/propagated/go_terms.propagated.tsv gsea_propagated/${SP}.fst pos
done

# 5. Compare with the direct-annotation results in gsea/
python3 "$CHECKDIR/compare_propagated.py" --direct gsea --propagated gsea_propagated \
    --species cra:VF for:MGF par:STF --focal GO:0002250 GO:0007596 GO:0030168 GO:0045087 \
    --out gsea_propagated/comparison | tee gsea_propagated/comparison.log

echo
echo "Done. Key numbers: gsea_propagated/comparison.focal_terms.tsv (rank_mean_NES, direct vs propagated)."
