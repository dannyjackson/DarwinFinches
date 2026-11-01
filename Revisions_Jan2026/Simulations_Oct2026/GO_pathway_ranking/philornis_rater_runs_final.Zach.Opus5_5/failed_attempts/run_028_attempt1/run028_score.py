import sys
S = {}
def a(score, names):
    for n in names.strip().split("\n"):
        n = n.strip()
        if n: S[n] = score
a(5, """
immune response
adaptive immune response
B cell activation
B cell receptor signaling pathway
B cell differentiation
positive regulation of B cell proliferation
erythrocyte differentiation
positive regulation of erythrocyte differentiation
""")
a(4, """
immune system process
regulation of immune system process
regulation of immune response
innate immune response
positive regulation of innate immune response
activation of innate immune response
defense response
inflammatory response
regulation of inflammatory response
positive regulation of inflammatory response
negative regulation of inflammatory response
complement activation
T cell activation
positive regulation of T cell activation
T cell differentiation
T cell receptor signaling pathway
hemopoiesis
hematopoietic progenitor cell differentiation
blood coagulation
platelet aggregation
platelet activation
wound healing
positive regulation of wound healing
response to wounding
cytokine-mediated signaling pathway
regulation of cytokine production
positive regulation of cytokine production
positive regulation of interleukin-6 production
positive regulation of tumor necrosis factor production
positive regulation of interleukin-8 production
cellular response to interleukin-1
interleukin-17-mediated signaling pathway
positive regulation of chemokine production
intracellular iron ion homeostasis
""")
a(3, """
negative regulation of cytokine production
negative regulation of interleukin-6 production
negative regulation of tumor necrosis factor production
positive regulation of interleukin-12 production
positive regulation of type II interferon production
cellular response to type II interferon
response to cytokine
cellular response to cytokine stimulus
response to tumor necrosis factor
cellular response to tumor necrosis factor
tumor necrosis factor-mediated signaling pathway
negative regulation of T cell activation
negative regulation of T cell receptor signaling pathway
negative regulation of innate immune response
leukocyte cell-cell adhesion
phagocytosis
phagocytosis, engulfment
positive regulation of phagocytosis
chemotaxis
cell chemotaxis
positive chemotaxis
canonical NF-kappaB signal transduction
regulation of canonical NF-kappaB signal transduction
positive regulation of canonical NF-kappaB signal transduction
negative regulation of canonical NF-kappaB signal transduction
non-canonical NF-kappaB signal transduction
positive regulation of non-canonical NF-kappaB signal transduction
negative regulation of non-canonical NF-kappaB signal transduction
toll-like receptor signaling pathway
toll-like receptor 4 signaling pathway
cytoplasmic pattern recognition receptor signaling pathway
positive regulation of NLRP3 inflammasome complex assembly
defense response to bacterium
defense response to Gram-negative bacterium
cell surface receptor signaling pathway via JAK-STAT
positive regulation of receptor signaling pathway via JAK-STAT
negative regulation of receptor signaling pathway via JAK-STAT
cell surface receptor signaling pathway via STAT
superoxide anion generation
response to hypoxia
cellular response to hypoxia
angiogenesis
regulation of angiogenesis
positive regulation of angiogenesis
sprouting angiogenesis
keratinocyte differentiation
positive regulation of keratinocyte differentiation
epidermis development
skin development
collagen fibril organization
extracellular matrix organization
positive regulation of fibroblast proliferation
positive regulation of epithelial cell migration
transforming growth factor beta receptor signaling pathway
cellular response to transforming growth factor beta stimulus
osteoclast differentiation
positive regulation of osteoclast differentiation
glucose homeostasis
regulation of glucose metabolic process
gluconeogenesis
energy homeostasis
response to starvation
cellular response to starvation
response to nutrient levels
cellular response to nutrient levels
nuclear receptor-mediated steroid hormone signaling pathway
response to oxidative stress
cellular response to oxidative stress
cellular oxidant detoxification
sensory perception of pain
response to stress
regulation of reactive oxygen species metabolic process
positive regulation of reactive oxygen species metabolic process
insulin receptor signaling pathway
regulation of blood pressure
""")
a(2, """
cellular response to lipopolysaccharide
response to lipopolysaccharide
positive regulation of type I interferon production
negative regulation of type I interferon production
positive regulation of interferon-beta production
positive regulation of interferon-alpha production
negative regulation of type I interferon-mediated signaling pathway
antiviral innate immune response
positive regulation of defense response to virus by host
defense response to virus
response to virus
cellular response to virus
platelet-derived
regulation of heart rate
regulation of heart rate by cardiac conduction
blood vessel development
cell migration involved in sprouting angiogenesis
positive regulation of cell migration involved in sprouting angiogenesis
negative regulation of cell migration involved in sprouting angiogenesis
negative regulation of angiogenesis
vascular endothelial growth factor signaling pathway
vascular endothelial growth factor receptor signaling pathway
cellular response to vascular endothelial growth factor stimulus
positive regulation of vascular endothelial growth factor production
establishment of endothelial barrier
endothelial cell migration
positive regulation of endothelial cell migration
positive regulation of blood vessel endothelial cell migration
positive regulation of endothelial cell proliferation
basement membrane organization
plasma membrane repair
epithelial to mesenchymal transition
positive regulation of epithelial cell proliferation
TGF-beta
positive regulation of transforming growth factor beta receptor signaling pathway
negative regulation of transforming growth factor beta receptor signaling pathway
SMAD protein signal transduction
cell migration
positive regulation of cell migration
regulation of cell migration
lamellipodium assembly
bone mineralization
regulation of bone mineralization
positive regulation of bone mineralization
ossification
negative regulation of ossification
bone development
osteoblast differentiation
positive regulation of osteoblast differentiation
negative regulation of osteoblast differentiation
skeletal system development
BMP signaling pathway
roof of mouth development
insulin secretion
regulation of insulin secretion
positive regulation of insulin secretion
response to insulin
cellular response to insulin stimulus
glucose metabolic process
intracellular glucose homeostasis
cellular response to glucose starvation
cellular response to amino acid starvation
insulin-like growth factor receptor signaling pathway
lipid storage
triglyceride homeostasis
glycogen biosynthetic process
fatty acid beta-oxidation
positive regulation of cold-induced thermogenesis
negative regulation of cold-induced thermogenesis
TOR signaling
regulation of cell growth
positive regulation of cell growth
response to hormone
hormone-mediated signaling pathway
cellular response to hormone stimulus
neuropeptide signaling pathway
learning
memory
cognition
cellular response to stress
cellular response to hydrogen peroxide
cellular response to reactive oxygen species
cell redox homeostasis
negative regulation of ferroptosis
metal ion transport
response to xenobiotic stimulus
cellular response to xenobiotic stimulus
xenobiotic transport
extrinsic apoptotic signaling pathway
symbiont entry into host cell
cellular response to mechanical stimulus
proteolysis
positive regulation of proteolysis
heterotypic cell-cell adhesion
cell adhesion mediated by integrin
integrin-mediated signaling pathway
cell-matrix adhesion
regulated exocytosis
positive regulation of protein secretion
protein secretion
receptor-mediated endocytosis
response to calcium ion
JNK cascade
p38MAPK
positive regulation of p38MAPK cascade
positive regulation of JNK cascade
fibroblast growth factor receptor signaling pathway
epidermal growth factor receptor signaling pathway
cellular response to epidermal growth factor stimulus
morphogenesis of an epithelium
epithelial cell differentiation
muscle contraction
heart development
cholesterol metabolic process
retinol metabolic process
cell population proliferation
negative regulation of epithelial cell proliferation
glycolytic process
mitochondrial calcium ion homeostasis
response to heat
cellular response to heat
""")
lines = open(sys.argv[1]).read().splitlines()[1:]
out = ["pathway\tscore"]
used=set()
from collections import Counter
c=Counter()
for l in lines:
    if not l.strip(): continue
    gid, name = l.split(" ",1)
    s = S.get(name.strip(), 1)
    if name.strip() in S: used.add(name.strip())
    c[s]+=1
    out.append(f"{gid}\t{s}")
open(sys.argv[2],"w").write("\n".join(out)+"\n")
print(len(out)-1, c, file=sys.stderr)
print("unmatched:", [k for k in S if k not in used], file=sys.stderr)
