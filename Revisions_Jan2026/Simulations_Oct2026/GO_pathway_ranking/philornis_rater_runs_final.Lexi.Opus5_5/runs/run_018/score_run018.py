#!/usr/bin/env python3
# Scores for run_018: explicit judgments for terms scored >1; all other terms
# were reviewed and judged to have no specific link (score 1).
import sys

S = {}
def put(score, names):
    for n in names.strip().split("\n"):
        n = n.strip()
        if n:
            S[n] = score

# 5: directly captures documented host mechanisms (IgY antibody/adaptive immunity,
# inflammation/immune response to larval feeding, blood loss/clotting, wound/tissue damage)
put(5, """
immune response
adaptive immune response
B cell activation
B cell receptor signaling pathway
inflammatory response
blood coagulation
response to wounding
wound healing
positive regulation of wound healing
hemopoiesis
erythrocyte differentiation
positive regulation of erythrocyte differentiation
defense response
""")

# 4: core components / strong a priori predictions
put(4, """
B cell differentiation
positive regulation of B cell proliferation
T cell activation
T cell differentiation
T cell receptor signaling pathway
positive regulation of T cell activation
negative regulation of T cell activation
negative regulation of T cell receptor signaling pathway
innate immune response
activation of innate immune response
positive regulation of innate immune response
negative regulation of innate immune response
immune system process
regulation of immune response
regulation of immune system process
complement activation
cytokine-mediated signaling pathway
positive regulation of cytokine production
negative regulation of cytokine production
regulation of cytokine production
positive regulation of inflammatory response
negative regulation of inflammatory response
regulation of inflammatory response
platelet activation
platelet aggregation
hematopoietic progenitor cell differentiation
intracellular iron ion homeostasis
positive regulation of interleukin-6 production
negative regulation of interleukin-6 production
positive regulation of interleukin-8 production
positive regulation of interleukin-12 production
positive regulation of tumor necrosis factor production
negative regulation of tumor necrosis factor production
positive regulation of chemokine production
positive regulation of type II interferon production
interleukin-17-mediated signaling pathway
tumor necrosis factor-mediated signaling pathway
leukocyte cell-cell adhesion
phagocytosis
phagocytosis, engulfment
positive regulation of phagocytosis
canonical NF-kappaB signal transduction
positive regulation of canonical NF-kappaB signal transduction
negative regulation of canonical NF-kappaB signal transduction
regulation of canonical NF-kappaB signal transduction
cell chemotaxis
keratinocyte differentiation
positive regulation of keratinocyte differentiation
epidermis development
skin development
response to cytokine
cellular response to cytokine stimulus
response to tumor necrosis factor
cellular response to tumor necrosis factor
cellular response to interleukin-1
cellular response to type II interferon
positive regulation of NLRP3 inflammasome complex assembly
""")

# 3: plausible downstream / supporting processes
put(3, """
defense response to bacterium
defense response to Gram-negative bacterium
response to lipopolysaccharide
cellular response to lipopolysaccharide
toll-like receptor signaling pathway
toll-like receptor 4 signaling pathway
cytoplasmic pattern recognition receptor signaling pathway
non-canonical NF-kappaB signal transduction
positive regulation of non-canonical NF-kappaB signal transduction
negative regulation of non-canonical NF-kappaB signal transduction
cell surface receptor signaling pathway via JAK-STAT
cell surface receptor signaling pathway via STAT
positive regulation of receptor signaling pathway via JAK-STAT
negative regulation of receptor signaling pathway via JAK-STAT
positive regulation of p38MAPK cascade
angiogenesis
sprouting angiogenesis
positive regulation of angiogenesis
negative regulation of angiogenesis
regulation of angiogenesis
blood vessel development
collagen fibril organization
extracellular matrix organization
basement membrane organization
positive regulation of fibroblast proliferation
positive regulation of epithelial cell migration
positive regulation of epithelial cell proliferation
epithelial cell differentiation
establishment of endothelial barrier
transforming growth factor beta receptor signaling pathway
cellular response to transforming growth factor beta stimulus
positive regulation of transforming growth factor beta receptor signaling pathway
negative regulation of transforming growth factor beta receptor signaling pathway
vascular endothelial growth factor receptor signaling pathway
vascular endothelial growth factor signaling pathway
positive regulation of vascular endothelial growth factor production
cellular response to vascular endothelial growth factor stimulus
epidermal growth factor receptor signaling pathway
cellular response to epidermal growth factor stimulus
plasma membrane repair
cellular response to hypoxia
response to hypoxia
response to oxidative stress
cellular response to oxidative stress
cellular oxidant detoxification
cellular response to reactive oxygen species
cellular response to hydrogen peroxide
superoxide anion generation
positive regulation of reactive oxygen species metabolic process
regulation of reactive oxygen species metabolic process
cell redox homeostasis
response to stress
cellular response to stress
response to hormone
cellular response to hormone stimulus
hormone-mediated signaling pathway
nuclear receptor-mediated steroid hormone signaling pathway
glucose homeostasis
intracellular glucose homeostasis
energy homeostasis
response to starvation
cellular response to starvation
response to nutrient levels
cellular response to nutrient levels
cellular response to glucose starvation
cellular response to amino acid starvation
gluconeogenesis
lipid storage
triglyceride homeostasis
fatty acid beta-oxidation
insulin-like growth factor receptor signaling pathway
bone development
ossification
osteoblast differentiation
osteoclast differentiation
positive regulation of osteoblast differentiation
negative regulation of osteoblast differentiation
positive regulation of osteoclast differentiation
negative regulation of ossification
bone mineralization
positive regulation of bone mineralization
regulation of bone mineralization
regulation of blood pressure
erythrocyte homeostasis
leukocyte migration
chemotaxis
positive chemotaxis
cell migration
positive regulation of cell migration
endothelial cell migration
positive regulation of endothelial cell migration
positive regulation of blood vessel endothelial cell migration
positive regulation of endothelial cell proliferation
sensory perception of pain
heme biosynthetic process
regulated exocytosis
positive regulation of protein secretion
apoptotic process
cellular response to mechanical stimulus
antiviral innate immune response
defense response to virus
response to virus
cellular response to virus
positive regulation of defense response to virus by host
positive regulation of type I interferon production
negative regulation of type I interferon production
negative regulation of type I interferon-mediated signaling pathway
positive regulation of interferon-alpha production
positive regulation of interferon-beta production
response to xenobiotic stimulus
cellular response to xenobiotic stimulus
xenobiotic transport
""")

# 2: weak, indirect or speculative links (growth/development tolerance, beak/skull
# morphogenesis, metabolism supporting compensation, behaviour, general signalling)
put(2, """
skeletal system development
cartilage development
chondrocyte differentiation
negative regulation of chondrocyte differentiation
roof of mouth development
BMP signaling pathway
cellular response to BMP stimulus
positive regulation of BMP signaling pathway
negative regulation of BMP signaling pathway
embryonic skeletal system morphogenesis
embryonic skeletal system development
skeletal muscle tissue development
muscle organ development
myoblast differentiation
positive regulation of myoblast differentiation
regulation of cell growth
positive regulation of cell growth
negative regulation of cell growth
cell population proliferation
positive regulation of cell population proliferation
negative regulation of cell population proliferation
regulation of cell population proliferation
positive regulation of smooth muscle cell proliferation
positive regulation of vascular associated smooth muscle cell proliferation
fat cell differentiation
positive regulation of fat cell differentiation
negative regulation of fat cell differentiation
lipid metabolic process
lipid catabolic process
lipid homeostasis
regulation of lipid metabolic process
fatty acid metabolic process
long-chain fatty acid metabolic process
fatty acid biosynthetic process
triglyceride biosynthetic process
glucose metabolic process
regulation of glucose metabolic process
carbohydrate metabolic process
glycolytic process
glycogen biosynthetic process
positive regulation of D-glucose import
cellular response to glucose stimulus
insulin receptor signaling pathway
negative regulation of insulin receptor signaling pathway
insulin secretion
regulation of insulin secretion
positive regulation of insulin secretion
negative regulation of insulin secretion
positive regulation of insulin secretion involved in cellular response to glucose stimulus
response to insulin
cellular response to insulin stimulus
amino acid metabolic process
TOR signaling
TORC1 signaling
positive regulation of TOR signaling
negative regulation of TOR signaling
positive regulation of TORC1 signaling
negative regulation of TORC1 signaling
aerobic respiration
tricarboxylic acid cycle
ATP metabolic process
mitochondrial electron transport, NADH to ubiquinone
mitochondrial respiratory chain complex I assembly
positive regulation of cold-induced thermogenesis
negative regulation of cold-induced thermogenesis
cholesterol homeostasis
cholesterol metabolic process
cholesterol biosynthetic process
response to heat
cellular response to heat
circadian rhythm
learning
memory
cognition
neuropeptide signaling pathway
regulation of heart rate
response to estrogen
cellular response to estrogen stimulus
cellular response to estradiol stimulus
response to progesterone
androgen receptor signaling pathway
intracellular receptor signaling pathway
retinoic acid receptor signaling pathway
retinol metabolic process
response to retinoic acid
cellular response to retinoic acid
intracellular zinc ion homeostasis
zinc ion transport
zinc ion transmembrane transport
calcium ion homeostasis
cell adhesion mediated by integrin
integrin-mediated signaling pathway
cell-matrix adhesion
positive regulation of cell-matrix adhesion
focal adhesion assembly
cell adhesion
regulation of cell migration
negative regulation of cell migration
cell motility
fibroblast growth factor receptor signaling pathway
negative regulation of fibroblast growth factor receptor signaling pathway
epithelial to mesenchymal transition
positive regulation of epithelial to mesenchymal transition
negative regulation of epithelial to mesenchymal transition
morphogenesis of an epithelium
bicellular tight junction assembly
cell-cell junction assembly
negative regulation of endothelial cell apoptotic process
positive regulation of endothelial cell apoptotic process
negative regulation of endothelial cell proliferation
cell migration involved in sprouting angiogenesis
positive regulation of cell migration involved in sprouting angiogenesis
negative regulation of cell migration involved in sprouting angiogenesis
apoptotic signaling pathway
extrinsic apoptotic signaling pathway
positive regulation of extrinsic apoptotic signaling pathway
negative regulation of extrinsic apoptotic signaling pathway
negative regulation of extrinsic apoptotic signaling pathway via death domain receptors
regulation of apoptotic process
positive regulation of apoptotic process
negative regulation of apoptotic process
negative regulation of ferroptosis
autophagy
macroautophagy
endocytosis
receptor-mediated endocytosis
exocytosis
protein secretion
lysosome organization
lysosomal lumen acidification
proteolysis
membrane protein ectodomain proteolysis
JNK cascade
positive regulation of JNK cascade
negative regulation of JNK cascade
MAPK cascade
positive regulation of MAPK cascade
negative regulation of MAPK cascade
regulation of MAPK cascade
positive regulation of ERK1 and ERK2 cascade
negative regulation of ERK1 and ERK2 cascade
regulation of ERK1 and ERK2 cascade
phosphatidylinositol 3-kinase/protein kinase B signal transduction
positive regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction
negative regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction
calcium-mediated signaling
positive regulation of calcium-mediated signaling
G protein-coupled purinergic nucleotide receptor signaling pathway
response to ATP
cellular response to cAMP
positive regulation of exosomal secretion
endoplasmic reticulum unfolded protein response
response to endoplasmic reticulum stress
cellular senescence
negative regulation of cellular senescence
DNA damage response
regulation of pH
bicarbonate transport
cardiac muscle contraction
heart development
kidney development
lung development
animal organ development
tissue development
developmental process
regulation of multicellular organismal process
positive regulation of multicellular organismal process
anatomical structure morphogenesis
animal organ morphogenesis
symbiont entry into host cell
heparan sulfate proteoglycan biosynthetic process
chondroitin sulfate proteoglycan biosynthetic process
proteoglycan biosynthetic process
intermediate filament organization
oligosaccharide metabolic process
protein N-linked glycosylation
protein O-linked glycosylation
ceramide biosynthetic process
sphingolipid biosynthetic process
cell surface receptor signaling pathway
signaling
cell communication
cell-cell signaling
""")

def main(inp, out):
    with open(inp) as f:
        lines = f.read().splitlines()
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        gid, _, name = line.partition(" ")
        rows.append((gid, S.get(name.strip(), 1)))
    with open(out, "w") as f:
        f.write("pathway\tscore\n")
        for gid, s in rows:
            f.write(f"{gid}\t{s}\n")
    names = set(l.partition(" ")[2].strip() for l in lines[1:])
    missing = [k for k in S if k not in names]
    print("rows", len(rows), "unmatched keys", missing, file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
