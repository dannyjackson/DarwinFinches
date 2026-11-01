import re
inp = "/home/claude/runs/run_005/input_terms.txt"
out = "/home/claude/runs/run_005/scores.tsv"

# Ordered rules: first match wins. (regex on lowercase term name, score)
rules = [
    # 5: documented mechanisms (antibody response, blood loss/erythropoiesis, tissue damage/repair)
    (r"^adaptive immune response$", 5),
    (r"^b cell (activation|differentiation|receptor signaling pathway)$", 5),
    (r"positive regulation of b cell proliferation", 5),
    (r"^immune response$", 5),
    (r"erythrocyte differentiation", 5),
    (r"^hemopoiesis$", 5),
    (r"^wound healing$|^response to wounding$|positive regulation of wound healing", 5),
    # 4: core components / strong predictions
    (r"^blood coagulation$|^platelet (activation|aggregation)$", 4),
    (r"hematopoietic progenitor", 4),
    (r"^innate immune response$|activation of innate immune response|regulation of innate immune response", 4),
    (r"^inflammatory response$|regulation of inflammatory response", 4),
    (r"^defense response$|^immune system process$|regulation of immune (system process|response)", 4),
    (r"complement activation", 4),
    (r"\bt cell (activation|differentiation|receptor signaling)|regulation of t cell", 4),
    (r"leukocyte cell-cell adhesion", 4),
    (r"intracellular iron ion homeostasis", 4),
    (r"interleukin|tumor necrosis factor|type ii interferon|cytokine(?!sis)", 4),
    (r"nf-kappab", 3),
    (r"toll-like receptor|pattern recognition|nlrp3", 3),
    (r"jak-stat|via stat", 3),
    # 3: plausible downstream/supporting
    (r"^phagocytosis|regulation of phagocytosis", 3),
    (r"chemotaxis", 3),
    (r"^angiogenesis$|regulation of angiogenesis|sprouting angiogenesis|vascular endothelial growth factor|blood vessel development", 3),
    (r"^collagen fibril organization$|^extracellular matrix organization$|basement membrane organization", 3),
    (r"keratinocyte differentiation|^epidermis development$|^skin development$", 3),
    (r"fibroblast proliferation", 3),
    (r"transforming growth factor beta", 3),
    (r"oxidative stress|reactive oxygen species|hydrogen peroxide|oxidant detoxification|cell redox homeostasis|superoxide", 3),
    (r"response to hypoxia|cellular response to hypoxia", 3),
    (r"glucose homeostasis|^glucose metabolic process$|regulation of glucose metabolic process|^energy homeostasis$", 3),
    (r"starvation|response to nutrient levels|cellular response to nutrient levels", 3),
    (r"sensory perception of pain", 3),
    (r"cytokine(?!sis) production", 4),
    (r"defense response to (gram-negative )?bacterium|response to lipopolysaccharide|cellular response to lipopolysaccharide", 2),
    (r"plasma membrane repair", 3),
    (r"^epithelial cell migration|positive regulation of epithelial cell migration", 3),
    (r"^response to stress$", 3),
    (r"regulation of blood pressure", 2),
    # 2: weak/speculative
    (r"virus|viral|interferon", 2),
    (r"symbiont entry", 2),
    (r"cell migration|cell motility", 2),
    (r"endothelial", 2),
    (r"insulin|igf|insulin-like growth factor|tor signaling|torc1", 2),
    (r"response to hormone|hormone-mediated|hormone stimulus|steroid hormone", 2),
    (r"^circadian rhythm$|regulation of circadian rhythm|entrainment of circadian", 2),
    (r"^learning$|^memory$|^cognition$", 2),
    (r"bone mineralization|^ossification$|osteoblast differentiation|^bone development$|osteoclast", 2),
    (r"^skeletal system development$", 2),
    (r"cold-induced thermogenesis", 2),
    (r"lipid storage|triglyceride|fat cell differentiation|fatty acid beta-oxidation|gluconeogenesis|glycogen biosynthetic", 2),
    (r"^apoptotic process$|extrinsic apoptotic", 2),
    (r"cellular senescence|telomere maintenance", 2),
    (r"response to heat|cellular response to heat", 2),
    (r"cellular response to stress", 2),
    (r"growth factor stimulus|regulation of cell growth|cell population proliferation", 2),
    (r"xenobiotic", 2),
    (r"cell-matrix adhesion|cell adhesion mediated by integrin|integrin-mediated", 2),
    (r"regulated exocytosis|protein secretion", 2),
    (r"heparan sulfate|chondroitin sulfate", 2),
    (r"p38mapk|jnk cascade", 2),
    (r"zinc ion|iron-sulfur", 2),
    (r"heart rate|cardiac muscle contraction", 2),
    (r"lymph|leukocyte|myeloid|macrophage|neutrophil|granulocyte", 3),
]

lines = open(inp).read().splitlines()[1:]
rows = []
for ln in lines:
    if not ln.strip():
        continue
    gid, _, name = ln.partition(" ")
    n = name.lower().strip()
    s = 1
    for pat, sc in rules:
        if re.search(pat, n):
            s = sc
            break
    rows.append((gid, s, name))

with open(out, "w") as f:
    f.write("pathway\tscore\n")
    for gid, s, _ in rows:
        f.write(f"{gid}\t{s}\n")

import collections
print(len(rows), collections.Counter(r[1] for r in rows))
for gid, s, name in rows:
    if s >= 2:
        print(s, gid, name)
