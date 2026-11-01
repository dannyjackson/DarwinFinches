import sys
S = {}
def put(score, names):
    for n in names.strip().split("\n"):
        n = n.strip()
        if n: S[n] = score
put(5, """
adaptive immune response
B cell activation
B cell receptor signaling pathway
B cell differentiation
positive regulation of B cell proliferation
immune response
wound healing
response to wounding
positive regulation of wound healing
inflammatory response
erythrocyte differentiation
positive regulation of erythrocyte differentiation
hemopoiesis
blood coagulation
""")
put(4, """
immune system process
innate immune response
defense response
complement activation
regulation of immune response
regulation of immune system process
regulation of inflammatory response
positive regulation of inflammatory response
negative regulation of inflammatory response
T cell activation
positive regulation of T cell activation
T cell receptor signaling pathway
T cell differentiation
cytokine-mediated signaling pathway
platelet activation
platelet aggregation
intracellular iron ion homeostasis
hematopoietic progenitor cell differentiation
activation of innate immune response
positive regulation of innate immune response
plasma membrane repair
skin development
epidermis development
""")
put(3, """
cellular response to cytokine stimulus
response to cytokine
regulation of cytokine production
positive regulation of cytokine production
negative regulation of cytokine production
positive regulation of tumor necrosis factor production
negative regulation of tumor necrosis factor production
tumor necrosis factor-mediated signaling pathway
response to tumor necrosis factor
cellular response to tumor necrosis factor
positive regulation of interleukin-6 production
negative regulation of interleukin-6 production
positive regulation of interleukin-8 production
positive regulation of interleukin-12 production
interleukin-17-mediated signaling pathway
cellular response to interleukin-1
positive regulation of type II interferon production
cellular response to type II interferon
positive regulation of chemokine production
negative regulation of innate immune response
negative regulation of T cell activation
negative regulation of T cell receptor signaling pathway
positive regulation of NLRP3 inflammasome complex assembly
toll-like receptor signaling pathway
toll-like receptor 4 signaling pathway
cytoplasmic pattern recognition receptor signaling pathway
canonical NF-kappaB signal transduction
regulation of canonical NF-kappaB signal transduction
positive regulation of canonical NF-kappaB signal transduction
negative regulation of canonical NF-kappaB signal transduction
non-canonical NF-kappaB signal transduction
positive regulation of non-canonical NF-kappaB signal transduction
negative regulation of non-canonical NF-kappaB signal transduction
cell surface receptor signaling pathway via JAK-STAT
positive regulation of receptor signaling pathway via JAK-STAT
negative regulation of receptor signaling pathway via JAK-STAT
cell surface receptor signaling pathway via STAT
leukocyte cell-cell adhesion
phagocytosis
phagocytosis, engulfment
positive regulation of phagocytosis
superoxide anion generation
defense response to bacterium
keratinocyte differentiation
positive regulation of keratinocyte differentiation
collagen fibril organization
extracellular matrix organization
angiogenesis
sprouting angiogenesis
positive regulation of angiogenesis
regulation of angiogenesis
blood vessel development
response to hypoxia
cellular response to hypoxia
response to oxidative stress
cellular response to oxidative stress
cellular oxidant detoxification
cell redox homeostasis
sensory perception of pain
energy homeostasis
glucose homeostasis
glucose metabolic process
regulation of glucose metabolic process
response to starvation
cellular response to starvation
response to nutrient levels
cellular response to nutrient levels
transforming growth factor beta receptor signaling pathway
cellular response to transforming growth factor beta stimulus
positive regulation of fibroblast proliferation
establishment of endothelial barrier
positive regulation of cholesterol efflux
""")
put(2, """
cellular response to lipopolysaccharide
response to lipopolysaccharide
defense response to Gram-negative bacterium
defense response to virus
response to virus
cellular response to virus
antiviral innate immune response
positive regulation of defense response to virus by host
positive regulation of type I interferon production
negative regulation of type I interferon production
negative regulation of type I interferon-mediated signaling pathway
positive regulation of interferon-beta production
positive regulation of interferon-alpha production
chemotaxis
cell chemotaxis
positive chemotaxis
cell migration
regulation of cell migration
positive regulation of cell migration
endothelial cell migration
positive regulation of endothelial cell migration
positive regulation of endothelial cell proliferation
positive regulation of epithelial cell migration
positive regulation of epithelial cell proliferation
cell population proliferation
positive regulation of cell population proliferation
regulation of cell population proliferation
regulation of cell growth
positive regulation of cell growth
negative regulation of cell growth
insulin-like growth factor receptor signaling pathway
TOR signaling
TORC1 signaling
response to stress
cellular response to stress
regulation of blood pressure
regulation of heart rate
response to heat
cellular response to heat
response to xenobiotic stimulus
cellular response to xenobiotic stimulus
xenobiotic transport
positive regulation of cold-induced thermogenesis
negative regulation of cold-induced thermogenesis
learning
memory
cognition
skeletal system development
bone development
ossification
bone mineralization
positive regulation of bone mineralization
regulation of bone mineralization
roof of mouth development
cartilage development
osteoblast differentiation
osteoclast differentiation
positive regulation of osteoclast differentiation
epithelial cell differentiation
basement membrane organization
cell-matrix adhesion
integrin-mediated signaling pathway
apoptotic process
regulation of apoptotic process
gluconeogenesis
glycogen biosynthetic process
intracellular glucose homeostasis
insulin receptor signaling pathway
cellular response to insulin stimulus
response to insulin
insulin secretion
regulation of insulin secretion
positive regulation of insulin secretion
lipid storage
fatty acid beta-oxidation
triglyceride homeostasis
lipid homeostasis
regulation of reactive oxygen species metabolic process
positive regulation of reactive oxygen species metabolic process
receptor-mediated endocytosis
heparan sulfate proteoglycan biosynthetic process
response to hormone
cellular response to hormone stimulus
hormone-mediated signaling pathway
androgen receptor signaling pathway
nuclear receptor-mediated steroid hormone signaling pathway
response to ATP
G protein-coupled purinergic nucleotide receptor signaling pathway
vascular endothelial growth factor receptor signaling pathway
vascular endothelial growth factor signaling pathway
positive regulation of vascular endothelial growth factor production
cellular response to vascular endothelial growth factor stimulus
negative regulation of angiogenesis
positive regulation of cell migration involved in sprouting angiogenesis
negative regulation of cell migration involved in sprouting angiogenesis
cell migration involved in sprouting angiogenesis
positive regulation of blood vessel endothelial cell migration
epidermal growth factor receptor signaling pathway
cellular response to epidermal growth factor stimulus
fibroblast growth factor receptor signaling pathway
lysosome organization
cellular response to mechanical stimulus
heme
respiration
aerobic respiration
cellular response to reactive oxygen species
cellular response to hydrogen peroxide
negative regulation of ferroptosis
metal ion transport
response to stress
circadian rhythm
""")
lines = open(sys.argv[1]).read().splitlines()[1:]
out = ["pathway\tscore"]
for l in lines:
    if not l.strip(): continue
    gid, _, name = l.partition(" ")
    out.append(f"{gid}\t{S.get(name.strip(),1)}")
open(sys.argv[2],"w").write("\n".join(out)+"\n")
print(len(out)-1)
