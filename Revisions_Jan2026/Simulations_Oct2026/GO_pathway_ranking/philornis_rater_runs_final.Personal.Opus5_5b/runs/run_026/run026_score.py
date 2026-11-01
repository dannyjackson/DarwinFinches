import re
S = {}
def s(score, *names):
    for n in names: S[n] = score

# 5: directly documented mechanisms (blood loss/haemoglobin, wound, antibody, inflammation)
s(5,"blood coagulation","platelet aggregation","platelet activation","wound healing","response to wounding",
  "positive regulation of wound healing","erythrocyte differentiation","positive regulation of erythrocyte differentiation",
  "inflammatory response","adaptive immune response","immune response","B cell activation","B cell receptor signaling pathway",
  "positive regulation of B cell proliferation","B cell differentiation")
# 4
s(4,"innate immune response","immune system process","regulation of immune response","regulation of immune system process",
  "regulation of inflammatory response","positive regulation of inflammatory response","negative regulation of inflammatory response",
  "defense response","hemopoiesis","hematopoietic progenitor cell differentiation","intracellular iron ion homeostasis",
  "glucose homeostasis","intracellular glucose homeostasis","T cell activation","positive regulation of T cell activation","T cell receptor signaling pathway",
  "T cell differentiation","complement activation","cytokine-mediated signaling pathway","response to cytokine","cellular response to cytokine stimulus",
  "positive regulation of cytokine production","regulation of cytokine production","negative regulation of cytokine production",
  "positive regulation of interleukin-6 production","negative regulation of interleukin-6 production","positive regulation of interleukin-8 production",
  "positive regulation of tumor necrosis factor production","negative regulation of tumor necrosis factor production",
  "positive regulation of interleukin-12 production","positive regulation of chemokine production","positive regulation of type II interferon production",
  "tumor necrosis factor-mediated signaling pathway","response to tumor necrosis factor","cellular response to tumor necrosis factor",
  "cellular response to interleukin-1","interleukin-17-mediated signaling pathway","leukocyte cell-cell adhesion",
  "keratinocyte differentiation","positive regulation of keratinocyte differentiation","epidermis development","skin development",
  "response to hypoxia","cellular response to hypoxia","activation of innate immune response","positive regulation of innate immune response",
  "negative regulation of innate immune response","nuclear receptor-mediated steroid hormone signaling pathway","energy homeostasis",
  "positive regulation of NLRP3 inflammasome complex assembly","negative regulation of T cell activation","negative regulation of T cell receptor signaling pathway")
# 3
s(3,"positive regulation of canonical NF-kappaB signal transduction","regulation of canonical NF-kappaB signal transduction",
  "negative regulation of canonical NF-kappaB signal transduction","canonical NF-kappaB signal transduction","non-canonical NF-kappaB signal transduction",
  "positive regulation of non-canonical NF-kappaB signal transduction","negative regulation of non-canonical NF-kappaB signal transduction",
  "cell surface receptor signaling pathway via JAK-STAT","positive regulation of receptor signaling pathway via JAK-STAT",
  "negative regulation of receptor signaling pathway via JAK-STAT","cell surface receptor signaling pathway via STAT",
  "toll-like receptor signaling pathway","toll-like receptor 4 signaling pathway","phagocytosis","phagocytosis, engulfment",
  "positive regulation of phagocytosis","cell chemotaxis","chemotaxis","cellular response to type II interferon","superoxide anion generation",
  "angiogenesis","sprouting angiogenesis","regulation of angiogenesis","positive regulation of angiogenesis","negative regulation of angiogenesis",
  "vascular endothelial growth factor signaling pathway","vascular endothelial growth factor receptor signaling pathway",
  "cellular response to vascular endothelial growth factor stimulus","positive regulation of vascular endothelial growth factor production",
  "transforming growth factor beta receptor signaling pathway","cellular response to transforming growth factor beta stimulus",
  "collagen fibril organization","extracellular matrix organization","basement membrane organization","positive regulation of fibroblast proliferation",
  "establishment of endothelial barrier","positive regulation of epithelial cell proliferation","positive regulation of epithelial cell migration",
  "intermediate filament organization","glucose metabolic process","regulation of glucose metabolic process","gluconeogenesis",
  "response to insulin","cellular response to insulin stimulus","insulin receptor signaling pathway","insulin-like growth factor receptor signaling pathway",
  "response to starvation","cellular response to starvation","response to nutrient levels","cellular response to nutrient levels",
  "regulation of cell growth","lipid storage","fatty acid beta-oxidation","triglyceride homeostasis","response to hormone","cellular response to hormone stimulus",
  "hormone-mediated signaling pathway","sensory perception of pain","response to oxidative stress","cellular response to oxidative stress",
  "cellular oxidant detoxification","cell redox homeostasis","response to stress","cellular response to stress",
  "bone development","ossification","skeletal system development","roof of mouth development","bone mineralization",
  "regulation of bone mineralization","positive regulation of bone mineralization","negative regulation of ossification",
  "osteoclast differentiation","positive regulation of osteoclast differentiation","osteoblast differentiation",
  "regulation of blood pressure","regulation of heart rate","cytoplasmic pattern recognition receptor signaling pathway",
  "positive regulation of protein secretion","leukocyte migration","regulated exocytosis","learning","memory","heme metabolic process",
  "cellular response to lipopolysaccharide","response to lipopolysaccharide","defense response to bacterium","defense response to Gram-negative bacterium",
  "positive regulation of type I interferon production","glycogen biosynthetic process","positive regulation of D-glucose import",
  "TOR signaling","TORC1 signaling","positive regulation of cell growth","negative regulation of cell growth",
  "positive regulation of endothelial cell proliferation","cell-matrix adhesion","heterotypic cell-cell adhesion",
  "membrane protein ectodomain proteolysis","circadian rhythm")
# 2
s(2,"antiviral innate immune response","response to virus","cellular response to virus","defense response to virus",
  "positive regulation of defense response to virus by host","positive regulation of interferon-alpha production","positive regulation of interferon-beta production",
  "negative regulation of type I interferon production","negative regulation of type I interferon-mediated signaling pathway",
  "positive regulation by host of viral transcription","positive regulation of viral genome replication","symbiont entry into host cell",
  "response to xenobiotic stimulus","cellular response to xenobiotic stimulus","xenobiotic transport",
  "positive regulation of p38MAPK cascade","MAPK cascade","regulation of MAPK cascade","positive regulation of MAPK cascade","negative regulation of MAPK cascade",
  "JNK cascade","positive regulation of JNK cascade","negative regulation of JNK cascade","ERK1 and ERK2 cascade","regulation of ERK1 and ERK2 cascade",
  "positive regulation of ERK1 and ERK2 cascade","negative regulation of ERK1 and ERK2 cascade",
  "apoptotic process","regulation of apoptotic process","positive regulation of apoptotic process","negative regulation of apoptotic process",
  "extrinsic apoptotic signaling pathway","positive regulation of extrinsic apoptotic signaling pathway","negative regulation of extrinsic apoptotic signaling pathway",
  "negative regulation of extrinsic apoptotic signaling pathway via death domain receptors",
  "cell migration","regulation of cell migration","positive regulation of cell migration","negative regulation of cell migration","cell motility",
  "endothelial cell migration","positive regulation of endothelial cell migration","positive regulation of blood vessel endothelial cell migration",
  "cell migration involved in sprouting angiogenesis","positive regulation of cell migration involved in sprouting angiogenesis",
  "negative regulation of cell migration involved in sprouting angiogenesis","negative regulation of endothelial cell proliferation",
  "negative regulation of endothelial cell apoptotic process","positive regulation of endothelial cell apoptotic process",
  "blood vessel development","epidermal growth factor receptor signaling pathway","cellular response to epidermal growth factor stimulus",
  "positive regulation of epidermal growth factor receptor signaling pathway","negative regulation of epidermal growth factor receptor signaling pathway",
  "fibroblast growth factor receptor signaling pathway","negative regulation of fibroblast growth factor receptor signaling pathway",
  "positive regulation of transforming growth factor beta receptor signaling pathway","negative regulation of transforming growth factor beta receptor signaling pathway",
  "SMAD protein signal transduction","positive regulation of SMAD protein signal transduction","BMP signaling pathway","cellular response to BMP stimulus",
  "positive regulation of BMP signaling pathway","negative regulation of BMP signaling pathway",
  "epithelial to mesenchymal transition","positive regulation of epithelial to mesenchymal transition","negative regulation of epithelial to mesenchymal transition",
  "negative regulation of epithelial cell proliferation","epithelial cell differentiation","morphogenesis of an epithelium",
  "cell-cell junction assembly","bicellular tight junction assembly","cell adhesion","regulation of cell adhesion","positive regulation of cell adhesion",
  "negative regulation of cell adhesion","cell adhesion mediated by integrin","integrin-mediated signaling pathway","focal adhesion assembly",
  "cartilage development","chondrocyte differentiation","negative regulation of osteoblast differentiation","positive regulation of osteoblast differentiation",
  "embryonic skeletal system development","embryonic skeletal system morphogenesis","negative regulation of chondrocyte differentiation",
  "lipid metabolic process","fatty acid metabolic process","regulation of lipid metabolic process","lipid homeostasis","lipid catabolic process",
  "triglyceride biosynthetic process","glycolytic process","carbohydrate metabolic process","insulin secretion","regulation of insulin secretion",
  "positive regulation of insulin secretion","negative regulation of insulin secretion","positive regulation of insulin secretion involved in cellular response to glucose stimulus",
  "cellular response to glucose stimulus","negative regulation of insulin receptor signaling pathway","cellular response to glucose starvation",
  "cellular response to amino acid starvation","fat cell differentiation","positive regulation of fat cell differentiation","negative regulation of fat cell differentiation",
  "aerobic respiration","tricarboxylic acid cycle","ATP metabolic process","mitochondrial electron transport, NADH to ubiquinone",
  "positive regulation of cold-induced thermogenesis","negative regulation of cold-induced thermogenesis","response to heat","cellular response to heat",
  "positive regulation of TOR signaling","negative regulation of TOR signaling","positive regulation of TORC1 signaling","negative regulation of TORC1 signaling",
  "regulation of reactive oxygen species metabolic process","positive regulation of reactive oxygen species metabolic process","cellular response to reactive oxygen species",
  "cellular response to hydrogen peroxide","negative regulation of ferroptosis","metal ion transport","zinc ion transport","zinc ion transmembrane transport",
  "intracellular zinc ion homeostasis","cellular response to mechanical stimulus","skeletal muscle tissue development","muscle organ development",
  "skeletal muscle contraction","muscle contraction","heart development","cardiac muscle contraction","regulation of cardiac muscle contraction",
  "regulation of heart rate by cardiac conduction","response to estrogen","response to progesterone","androgen receptor signaling pathway",
  "cellular response to estrogen stimulus","cellular response to estradiol stimulus","intracellular receptor signaling pathway",
  "neuropeptide signaling pathway","cognition","nervous system process","entrainment of circadian clock by photoperiod","regulation of circadian rhythm",
  "circadian regulation of gene expression","positive regulation of cell population proliferation","cell population proliferation",
  "regulation of cell population proliferation","negative regulation of cell population proliferation","cellular senescence","negative regulation of cellular senescence",
  "plasma membrane repair","response to ATP","G protein-coupled purinergic nucleotide receptor signaling pathway","lysosome organization",
  "endocytosis","receptor-mediated endocytosis","exocytosis","protein secretion","autophagy","macroautophagy","response to endoplasmic reticulum stress",
  "endoplasmic reticulum unfolded protein response","positive regulation of protein secretion","sensory perception of sound",
  "positive regulation of smooth muscle cell proliferation","positive regulation of vascular associated smooth muscle cell proliferation",
  "positive regulation of exosomal secretion","negative chemotaxis","positive chemotaxis","positive regulation of cell cycle","limb development",
  "positive regulation of substrate adhesion-dependent cell spreading","substrate adhesion-dependent cell spreading","cell-substrate adhesion","positive regulation of cell-matrix adhesion",
  "heparan sulfate proteoglycan catabolic process","heparan sulfate proteoglycan biosynthetic process","proteoglycan biosynthetic process",
  "chondroitin sulfate proteoglycan biosynthetic process","oligosaccharide metabolic process","protein N-linked glycosylation","glycoprotein biosynthetic process",
  "regulation of pH","bicarbonate transport","cholesterol metabolic process","cholesterol homeostasis","cholesterol efflux","positive regulation of cholesterol efflux",
  "sphingolipid biosynthetic process","ceramide biosynthetic process","ceramide metabolic process","fatty acid biosynthetic process","long-chain fatty acid metabolic process",
  "fatty acid transport","lipid transport","amino acid metabolic process","retinol metabolic process","response to retinoic acid","cellular response to retinoic acid",
  "positive regulation of calcium-mediated signaling","calcium-mediated signaling","stress granule assembly","protein refolding","protein folding",
  "regulation of developmental process","developmental process","system development","tissue development","animal organ development",
  "anatomical structure morphogenesis","animal gross anatomical part developmental process","intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress",
  "regulation of signaling","signaling","cell communication","regulation of multicellular organismal process","positive regulation of multicellular organismal process",
  "visual perception","sensory organ development","camera-type eye development","eye development")

name = {}
ids=[]
for line in open('/home/claude/runs/run_026/input_terms.txt').read().splitlines()[1:]:
    if not line.strip(): continue
    gid, nm = line.split(' ',1)
    ids.append((gid,nm))
unknown=[k for k in S if k not in {n for _,n in ids}]
print("unmatched names:",unknown)
with open('/home/claude/runs/run_026/scores.tsv','w') as f:
    f.write("pathway\tscore\n")
    for gid,nm in ids:
        f.write(f"{gid}\t{S.get(nm,1)}\n")
from collections import Counter
print(len(ids), Counter(S.get(n,1) for _,n in ids))
