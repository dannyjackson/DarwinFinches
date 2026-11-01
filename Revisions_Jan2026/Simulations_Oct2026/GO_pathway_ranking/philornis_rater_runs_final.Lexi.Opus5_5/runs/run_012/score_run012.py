import sys
S = {}
def a(score, names):
    for n in names.strip().split("\n"):
        n = n.strip()
        if n: S[n] = score
a(5, """adaptive immune response
immune response
B cell activation""")
a(4, """B cell receptor signaling pathway
B cell differentiation
positive regulation of B cell proliferation
innate immune response
defense response
immune system process
regulation of immune system process
regulation of immune response
inflammatory response
complement activation
wound healing
positive regulation of wound healing
response to wounding
blood coagulation
erythrocyte differentiation
positive regulation of erythrocyte differentiation
hemopoiesis
intracellular iron ion homeostasis""")
a(3, """T cell activation
T cell differentiation
T cell receptor signaling pathway
positive regulation of T cell activation
negative regulation of T cell activation
negative regulation of T cell receptor signaling pathway
regulation of inflammatory response
positive regulation of inflammatory response
negative regulation of inflammatory response
negative regulation of innate immune response
positive regulation of innate immune response
activation of innate immune response
regulation of cytokine production
positive regulation of cytokine production
negative regulation of cytokine production
positive regulation of tumor necrosis factor production
negative regulation of tumor necrosis factor production
positive regulation of interleukin-6 production
negative regulation of interleukin-6 production
positive regulation of interleukin-8 production
positive regulation of interleukin-12 production
positive regulation of type II interferon production
positive regulation of chemokine production
positive regulation of NLRP3 inflammasome complex assembly
toll-like receptor signaling pathway
canonical NF-kappaB signal transduction
positive regulation of canonical NF-kappaB signal transduction
cytokine-mediated signaling pathway
response to cytokine
cellular response to cytokine stimulus
response to tumor necrosis factor
cellular response to tumor necrosis factor
tumor necrosis factor-mediated signaling pathway
cellular response to interleukin-1
interleukin-17-mediated signaling pathway
cellular response to type II interferon
cell surface receptor signaling pathway via JAK-STAT
cell surface receptor signaling pathway via STAT
leukocyte cell-cell adhesion
phagocytosis
phagocytosis, engulfment
positive regulation of phagocytosis
cell chemotaxis
superoxide anion generation
platelet activation
platelet aggregation
hematopoietic progenitor cell differentiation
response to hypoxia
cellular response to hypoxia
keratinocyte differentiation
positive regulation of keratinocyte differentiation
epidermis development
skin development
collagen fibril organization
extracellular matrix organization
angiogenesis
positive regulation of fibroblast proliferation
response to oxidative stress
cellular response to oxidative stress
sensory perception of pain
nuclear receptor-mediated steroid hormone signaling pathway
insulin-like growth factor receptor signaling pathway
glucose homeostasis
energy homeostasis
response to starvation
response to nutrient levels""")
a(2, """antiviral innate immune response
defense response to virus
response to virus
cellular response to virus
positive regulation of defense response to virus by host
negative regulation of type I interferon production
positive regulation of type I interferon production
positive regulation of interferon-beta production
positive regulation of interferon-alpha production
negative regulation of type I interferon-mediated signaling pathway
cytoplasmic pattern recognition receptor signaling pathway
toll-like receptor 4 signaling pathway
response to lipopolysaccharide
cellular response to lipopolysaccharide
defense response to bacterium
defense response to Gram-negative bacterium
regulation of canonical NF-kappaB signal transduction
negative regulation of canonical NF-kappaB signal transduction
non-canonical NF-kappaB signal transduction
positive regulation of non-canonical NF-kappaB signal transduction
negative regulation of non-canonical NF-kappaB signal transduction
positive regulation of receptor signaling pathway via JAK-STAT
negative regulation of receptor signaling pathway via JAK-STAT
chemotaxis
positive chemotaxis
regulated exocytosis
response to stress
cellular response to stress
regulation of reactive oxygen species metabolic process
positive regulation of reactive oxygen species metabolic process
cellular oxidant detoxification
cellular response to hydrogen peroxide
cellular response to reactive oxygen species
regulation of blood pressure
regulation of heart rate
sprouting angiogenesis
positive regulation of angiogenesis
regulation of angiogenesis
negative regulation of angiogenesis
vascular endothelial growth factor signaling pathway
transforming growth factor beta receptor signaling pathway
cellular response to transforming growth factor beta stimulus
epidermal growth factor receptor signaling pathway
positive regulation of epithelial cell migration
positive regulation of epithelial cell proliferation
establishment of endothelial barrier
bicellular tight junction assembly
basement membrane organization
plasma membrane repair
cell migration
bone development
ossification
roof of mouth development
learning
memory
cognition
sensory perception of sound
circadian rhythm
neuropeptide signaling pathway
response to hormone
cellular response to hormone stimulus
hormone-mediated signaling pathway
intracellular receptor signaling pathway
cellular response to starvation
cellular response to nutrient levels
insulin receptor signaling pathway
response to insulin
cellular response to insulin stimulus
intracellular glucose homeostasis
regulation of glucose metabolic process
glucose metabolic process
gluconeogenesis
glycogen biosynthetic process
lipid storage
triglyceride homeostasis
fatty acid beta-oxidation
positive regulation of cell growth
regulation of cell growth
TOR signaling
lipid homeostasis""")
lines = open("/home/claude/project/runs/run_012/input_terms.txt").read().splitlines()[1:]
used=set(); out=["pathway\tscore"]
for l in lines:
    if not l.strip(): continue
    parts=l.split(" ",1); gid=parts[0]; name=parts[1].strip() if len(parts)>1 else ""
    sc=S.get(name,1)
    if name in S: used.add(name)
    out.append(f"{gid}\t{sc}")
missing=set(S)-used
if missing: print("UNMATCHED:",missing,file=sys.stderr)
open("/home/claude/project/runs/run_012/scores.tsv","w").write("\n".join(out)+"\n")
print(len(out)-1)
