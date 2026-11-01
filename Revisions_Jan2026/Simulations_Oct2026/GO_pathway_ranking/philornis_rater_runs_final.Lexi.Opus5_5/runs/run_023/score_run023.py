import re
src = "/home/claude/project/runs/run_023/input_terms.txt"
out = "/home/claude/project/runs/run_023/scores.tsv"

# Exact-name scores (my judgments); checked first.
exact = {
 # 5: documented mechanisms (antibody-mediated resistance; blood loss / wound responses)
 "adaptive immune response": 5, "B cell activation": 5,
 "response to wounding": 5, "wound healing": 5,
 "blood coagulation": 4, "hemopoiesis": 4, "erythrocyte differentiation": 4,
 "B cell receptor signaling pathway": 4, "B cell differentiation": 4,
 "positive regulation of B cell proliferation": 4, "immune system process": 4,
 "defense response": 4, "innate immune response": 4, "immune response": 4,
 "regulation of immune response": 4, "regulation of immune system process": 4,
 "positive regulation of wound healing": 4, "platelet activation": 4,
 "regulation of inflammatory response": 4, "positive regulation of inflammatory response": 4,
 "inflammatory response": 4, "negative regulation of inflammatory response": 4,
 "complement activation": 4, "T cell activation": 4, "positive regulation of T cell activation": 4,
 "intracellular iron ion homeostasis": 4, "hematopoietic progenitor cell differentiation": 4,
 "activation of innate immune response": 4, "positive regulation of innate immune response": 4,
 "negative regulation of innate immune response": 3,
 "T cell receptor signaling pathway": 3, "negative regulation of T cell receptor signaling pathway": 3,
 "regulation of cytokine production": 3, "positive regulation of cytokine production": 3,
 "cytokine-mediated signaling pathway": 3, "response to cytokine": 3,
 "positive regulation of interleukin-6 production": 3, "negative regulation of interleukin-6 production": 3,
 "positive regulation of interleukin-8 production": 3, "positive regulation of interleukin-12 production": 3,
 "positive regulation of tumor necrosis factor production": 3, "negative regulation of tumor necrosis factor production": 3,
 "tumor necrosis factor-mediated signaling pathway": 3, "response to tumor necrosis factor": 3,
 "cellular response to interleukin-1": 3, "cellular response to type II interferon": 3,
 "positive regulation of chemokine production": 3, "leukocyte cell-cell adhesion": 3,
 "toll-like receptor signaling pathway": 3, "toll-like receptor 4 signaling pathway": 3,
 "canonical NF-kappaB signal transduction": 3, "regulation of canonical NF-kappaB signal transduction": 3,
 "positive regulation of canonical NF-kappaB signal transduction": 3,
 "negative regulation of canonical NF-kappaB signal transduction": 3,
 "non-canonical NF-kappaB signal transduction": 3,
 "positive regulation of NLRP3 inflammasome complex assembly": 3,
 "phagocytosis": 3, "positive regulation of phagocytosis": 3,
 "defense response to bacterium": 3, "response to lipopolysaccharide": 3,
 "cellular response to lipopolysaccharide": 3, "cytoplasmic pattern recognition receptor signaling pathway": 2,
 "superoxide anion generation": 3, "cell chemotaxis": 3, "chemotaxis": 3,
 "cell surface receptor signaling pathway via JAK-STAT": 3,
 "positive regulation of receptor signaling pathway via JAK-STAT": 3,
 "negative regulation of receptor signaling pathway via JAK-STAT": 3,
 "establishment of endothelial barrier": 3, "angiogenesis": 3, "sprouting angiogenesis": 3,
 "positive regulation of angiogenesis": 3, "regulation of angiogenesis": 3, "negative regulation of angiogenesis": 3,
 "response to hypoxia": 3, "cellular oxidant detoxification": 3, "cellular response to oxidative stress": 3,
 "cellular response to hydrogen peroxide": 3, "positive regulation of reactive oxygen species metabolic process": 3,
 "collagen fibril organization": 3, "keratinocyte differentiation": 3, "positive regulation of keratinocyte differentiation": 3,
 "epidermis development": 3, "skin development": 3, "plasma membrane repair": 3,
 "positive regulation of fibroblast proliferation": 3, "energy homeostasis": 3,
 "response to starvation": 3, "cellular response to starvation": 3, "response to nutrient levels": 3,
 "cellular response to nutrient levels": 3, "glucose homeostasis": 3,
 "basement membrane organization": 2, "positive regulation of epithelial cell migration": 2,
 "positive regulation of epithelial cell proliferation": 2, "epithelial cell differentiation": 2,
 "cell-matrix adhesion": 2, "positive regulation of cell migration": 2, "cell migration": 2,
 "regulation of cell migration": 2, "positive regulation of endothelial cell migration": 2,
 "positive regulation of blood vessel endothelial cell migration": 2,
 "positive regulation of cell migration involved in sprouting angiogenesis": 3,
 "negative regulation of cell migration involved in sprouting angiogenesis": 3,
 "cell migration involved in sprouting angiogenesis": 3,
 "endothelial cell migration": 2, "vascular endothelial growth factor receptor signaling pathway": 2,
 "vascular endothelial growth factor signaling pathway": 2,
 "cellular response to vascular endothelial growth factor stimulus": 2,
 "positive regulation of vascular endothelial growth factor production": 2,
 "transforming growth factor beta receptor signaling pathway": 2,
 "cellular response to transforming growth factor beta stimulus": 2,
 "negative regulation of transforming growth factor beta receptor signaling pathway": 2,
 "response to xenobiotic stimulus": 2, "cellular response to xenobiotic stimulus": 2,
 "cellular response to stress": 2, "glucose metabolic process": 2, "regulation of glucose metabolic process": 2,
 "lipid storage": 2, "triglyceride homeostasis": 2, "fatty acid beta-oxidation": 2,
 "gluconeogenesis": 2, "glycogen biosynthetic process": 2, "insulin receptor signaling pathway": 2,
 "insulin-like growth factor receptor signaling pathway": 2, "response to insulin": 2,
 "cellular response to insulin stimulus": 2, "lipid homeostasis": 2,
 "regulation of heart rate": 2, "blood vessel development": 2, "cell population proliferation": 1,
 "receptor-mediated endocytosis": 2, "endocytosis": 1, "antigen processing and presentation": 4,
 "regulation of pH": 1, "symbiont entry into host cell": 1, "response to virus": 1,
 "defense response to virus": 1, "antiviral innate immune response": 2, "cellular response to virus": 1,
 "positive regulation of defense response to virus by host": 1,
 "positive regulation of interferon-beta production": 2, "positive regulation of interferon-alpha production": 2,
 "positive regulation of type I interferon production": 2, "negative regulation of type I interferon production": 2,
 "skeletal system development": 2, "bone development": 2, "ossification": 2, "cartilage development": 2,
 "osteoblast differentiation": 2, "regulation of bone mineralization": 2, "positive regulation of bone mineralization": 2,
 "negative regulation of ossification": 2, "chondrocyte differentiation": 2,
 "cellular response to mechanical stimulus": 2, "extrinsic apoptotic signaling pathway": 2,
 "learning": 2, "memory": 2, "neuropeptide signaling pathway": 2, "circadian rhythm": 1,
 "sensory perception of sound": 1, "hormone-mediated signaling pathway": 2,
 "cellular response to hormone stimulus": 2, "nuclear receptor-mediated steroid hormone signaling pathway": 2,
 "cellular response to heat": 1, "positive regulation of cold-induced thermogenesis": 2,
 "negative regulation of cold-induced thermogenesis": 2, "regulation of cytosolic calcium ion concentration": 1,
 "protein secretion": 1, "regulated exocytosis": 2, "exocytosis": 1, "response to cytokine": 3,
 "regulation of cytokine production": 3, "membrane protein ectodomain proteolysis": 2,
 "proteolysis": 1, "cell adhesion mediated by integrin": 2, "integrin-mediated signaling pathway": 2,
 "focal adhesion assembly": 1, "substrate adhesion-dependent cell spreading": 2,
 "lamellipodium assembly": 1, "intermediate filament organization": 2,
 "fat cell differentiation": 1, "positive regulation of fat cell differentiation": 1,
 "stem cell population maintenance": 1, "stem cell differentiation": 2,
 "positive regulation of smooth muscle cell proliferation": 1,
 "positive regulation of vascular associated smooth muscle cell proliferation": 2,
 "negative regulation of endothelial cell apoptotic process": 2, "positive regulation of endothelial cell apoptotic process": 2,
 "negative regulation of endothelial cell proliferation": 2, "response to endoplasmic reticulum stress": 1,
 "cellular response to glucose starvation": 2, "cellular response to amino acid starvation": 2,
 "response to estrogen": 1, "metal ion transport": 1, "heparan sulfate proteoglycan biosynthetic process": 1,
 "negative regulation of ferroptosis": 2, "cellular response to growth factor stimulus": 1,
 "embryonic skeletal system morphogenesis": 1, "embryonic skeletal system development": 1,
 "roof of mouth development": 2,
}


extra = {
 "positive regulation of erythrocyte differentiation": 4, "negative regulation of type I interferon-mediated signaling pathway": 2,
 "platelet aggregation": 4, "extracellular matrix organization": 3, "defense response to Gram-negative bacterium": 3,
 "cellular response to reactive oxygen species": 3, "response to oxidative stress": 3, "cellular response to hypoxia": 3,
 "positive regulation of non-canonical NF-kappaB signal transduction": 3, "negative regulation of non-canonical NF-kappaB signal transduction": 3,
 "phagocytosis, engulfment": 3, "osteoclast differentiation": 2, "positive regulation of osteoclast differentiation": 2,
 "cell redox homeostasis": 2, "sensory perception of pain": 2, "regulation of blood pressure": 2,
 "epithelial to mesenchymal transition": 2, "positive regulation of epithelial to mesenchymal transition": 2,
 "negative regulation of epithelial to mesenchymal transition": 2, "cell surface receptor signaling pathway via STAT": 3,
 "xenobiotic transport": 2, "response to stress": 2, "cognition": 2, "bone mineralization": 2,
 "positive regulation of p38MAPK cascade": 2, "positive regulation of endothelial cell proliferation": 2,
 "regulation of reactive oxygen species metabolic process": 3, "intracellular glucose homeostasis": 2,
 "lipid catabolic process": 2, "triglyceride biosynthetic process": 2, "positive regulation of D-glucose import": 2,
 "bicellular tight junction assembly": 2, "positive chemotaxis": 2, "glycolytic process": 2, "JNK cascade": 2,
}
exact.update(extra)

# Pattern fallbacks for names not in the dict.
patterns = [
 (r"immun|lymphocyte|leukocyte|antigen|immunoglobulin|antibody|B cell|T cell|macrophage|neutrophil|mast cell|granulocyte|natural killer|myeloid", 3),
 (r"interleukin|chemokine|cytokine(?!sis)|interferon(?!.*vir)(?!-mediated)|tumor necrosis", 3),
 (r"inflamm", 4), (r"coagulation|hemostasis|platelet|fibrinolysis", 4),
 (r"erythro|hemoglobin|heme |iron(?!-sulfur)", 3), (r"wound", 4),
 (r"vir(al|us)", 1),
]

rows = []
with open(src) as f:
    next(f)
    for line in f:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        gid, _, name = line.partition(" ")
        s = exact.get(name)
        if s is None:
            s = 1
            for p, v in patterns:
                if re.search(p, name):
                    s = v
                    break
        rows.append((gid, name, s))

with open(out, "w") as o:
    o.write("pathway\tscore\n")
    for gid, name, s in rows:
        o.write(f"{gid}\t{s}\n")
print(len(rows))
import collections
print(collections.Counter(r[2] for r in rows))
for gid, name, s in rows:
    if s >= 2 and name not in exact:
        print("PATTERN", s, name)
