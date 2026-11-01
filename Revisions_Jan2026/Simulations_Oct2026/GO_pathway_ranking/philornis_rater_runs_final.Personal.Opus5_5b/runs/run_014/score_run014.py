import re
S = {}
def s(v, *names):
    for n in names: S[n] = v
s(5, "adaptive immune response","immune response","B cell activation","wound healing","response to wounding")
s(4, "blood coagulation","platelet activation","platelet aggregation","inflammatory response","hemopoiesis",
  "erythrocyte differentiation","positive regulation of erythrocyte differentiation","intracellular iron ion homeostasis",
  "B cell receptor signaling pathway","B cell differentiation","positive regulation of B cell proliferation",
  "positive regulation of wound healing","regulation of inflammatory response","positive regulation of inflammatory response",
  "innate immune response","complement activation","regulation of immune response","immune system process","defense response")
s(3, "regulation of immune system process","hematopoietic progenitor cell differentiation","T cell activation",
  "T cell receptor signaling pathway","T cell differentiation","positive regulation of T cell activation","negative regulation of T cell activation",
  "cytokine-mediated signaling pathway","response to cytokine","cellular response to cytokine stimulus",
  "regulation of cytokine production","positive regulation of cytokine production","negative regulation of cytokine production",
  "negative regulation of inflammatory response","positive regulation of innate immune response","activation of innate immune response",
  "negative regulation of innate immune response","canonical NF-kappaB signal transduction","regulation of canonical NF-kappaB signal transduction",
  "positive regulation of canonical NF-kappaB signal transduction","negative regulation of canonical NF-kappaB signal transduction",
  "positive regulation of interleukin-6 production","negative regulation of interleukin-6 production","positive regulation of interleukin-8 production",
  "positive regulation of interleukin-12 production","positive regulation of tumor necrosis factor production","negative regulation of tumor necrosis factor production",
  "tumor necrosis factor-mediated signaling pathway","response to tumor necrosis factor","cellular response to tumor necrosis factor",
  "cellular response to interleukin-1","interleukin-17-mediated signaling pathway","positive regulation of type II interferon production",
  "cellular response to type II interferon","positive regulation of NLRP3 inflammasome complex assembly","toll-like receptor signaling pathway",
  "toll-like receptor 4 signaling pathway","positive regulation of chemokine production","chemotaxis","cell chemotaxis","leukocyte cell-cell adhesion",
  "phagocytosis","phagocytosis, engulfment","positive regulation of phagocytosis","superoxide anion generation",
  "defense response to bacterium","defense response to Gram-negative bacterium","response to oxidative stress","cellular response to oxidative stress",
  "response to hypoxia","cellular response to hypoxia","regulation of blood pressure","glucose homeostasis","energy homeostasis",
  "response to starvation","cellular response to starvation","response to nutrient levels","nuclear receptor-mediated steroid hormone signaling pathway",
  "skin development","epidermis development","keratinocyte differentiation","positive regulation of keratinocyte differentiation",
  "positive regulation of epithelial cell proliferation","positive regulation of epithelial cell migration","positive regulation of fibroblast proliferation",
  "extracellular matrix organization","collagen fibril organization","angiogenesis","positive regulation of angiogenesis",
  "transforming growth factor beta receptor signaling pathway","establishment of endothelial barrier","sensory perception of pain")
s(2, "regulation of angiogenesis","negative regulation of angiogenesis","sprouting angiogenesis","endothelial cell migration",
  "positive regulation of endothelial cell migration","positive regulation of blood vessel endothelial cell migration","vascular endothelial growth factor signaling pathway",
  "vascular endothelial growth factor receptor signaling pathway","cellular response to vascular endothelial growth factor stimulus","positive regulation of vascular endothelial growth factor production",
  "positive regulation of endothelial cell proliferation","blood vessel development","epidermal growth factor receptor signaling pathway","cellular response to epidermal growth factor stimulus",
  "fibroblast growth factor receptor signaling pathway","cellular response to transforming growth factor beta stimulus","positive regulation of transforming growth factor beta receptor signaling pathway",
  "negative regulation of transforming growth factor beta receptor signaling pathway","SMAD protein signal transduction","epithelial to mesenchymal transition",
  "positive regulation of epithelial to mesenchymal transition","basement membrane organization","bicellular tight junction assembly","plasma membrane repair",
  "epithelial cell differentiation","cell migration","regulation of cell migration","positive regulation of cell migration","positive chemotaxis",
  "heterotypic cell-cell adhesion","cell adhesion mediated by integrin","integrin-mediated signaling pathway","response to lipopolysaccharide","cellular response to lipopolysaccharide",
  "cell surface receptor signaling pathway via JAK-STAT","cell surface receptor signaling pathway via STAT","positive regulation of receptor signaling pathway via JAK-STAT",
  "negative regulation of receptor signaling pathway via JAK-STAT","non-canonical NF-kappaB signal transduction","positive regulation of non-canonical NF-kappaB signal transduction",
  "negative regulation of non-canonical NF-kappaB signal transduction","cytoplasmic pattern recognition receptor signaling pathway","negative regulation of T cell receptor signaling pathway",
  "positive regulation of p38MAPK cascade","response to xenobiotic stimulus","cellular response to xenobiotic stimulus","xenobiotic transport",
  "response to stress","cellular response to stress","response to reactive oxygen species","cellular response to reactive oxygen species","cellular response to hydrogen peroxide",
  "cell redox homeostasis","cellular oxidant detoxification","regulation of reactive oxygen species metabolic process","positive regulation of reactive oxygen species metabolic process",
  "regulation of heart rate","gluconeogenesis","glucose metabolic process","glycogen biosynthetic process","lipid storage","triglyceride homeostasis","lipid catabolic process",
  "fatty acid beta-oxidation","response to insulin","insulin receptor signaling pathway","insulin-like growth factor receptor signaling pathway","positive regulation of cell growth",
  "regulation of cell growth","TOR signaling","cellular response to amino acid starvation","cellular response to glucose starvation","cellular response to nutrient levels",
  "response to hormone","hormone-mediated signaling pathway","cellular response to hormone stimulus","positive regulation of cold-induced thermogenesis","negative regulation of cold-induced thermogenesis",
  "circadian rhythm","learning","memory","cognition","neuropeptide signaling pathway","roof of mouth development","skeletal system development","bone development","ossification",
  "G protein-coupled purinergic nucleotide receptor signaling pathway","response to ATP","regulated exocytosis","positive regulation of type I interferon production",
  "extrinsic apoptotic signaling pathway","angiogenesis involved")
seen=set(); rows=[]
lines=open("/home/claude/runs/run_014/input_terms.txt").read().splitlines()[1:]
for l in lines:
    if not l.strip(): continue
    gid,name=l.split(" ",1)
    sc=S.get(name,1); seen.add(name); rows.append((gid,sc))
with open("/home/claude/runs/run_014/scores.tsv","w") as f:
    f.write("pathway\tscore\n")
    for g,sc in rows: f.write(f"{g}\t{sc}\n")
print(len(rows), [n for n in S if n not in seen])
