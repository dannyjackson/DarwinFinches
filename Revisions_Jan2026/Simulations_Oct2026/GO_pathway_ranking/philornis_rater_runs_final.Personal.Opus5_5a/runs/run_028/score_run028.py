import re
rules = [(r'cytokinesis',1),(r'antiviral innate immune',2),
# 5: documented mechanisms (blood feeding/anemia, tissue damage, IgY antibody response)
(r'^blood coagulation$',5),(r'^wound healing$',5),(r'positive regulation of wound healing',5),(r'^response to wounding$',5),
(r'^erythrocyte differentiation$',5),(r'positive regulation of erythrocyte differentiation',5),
(r'^adaptive immune response$',5),(r'^B cell activation$',5),(r'^B cell receptor signaling',5),(r'^immune response$',5),
# 4: core components
(r'^platelet (aggregation|activation)$',4),(r'^hemopoiesis$',4),(r'hematopoietic progenitor',4),
(r'^B cell differentiation$',4),(r'positive regulation of B cell proliferation',4),
(r'^inflammatory response$',4),(r'regulation of inflammatory response',4),(r'^defense response$',4),
(r'^innate immune response$',4),(r'innate immune response',3),(r'^complement activation$',4),
(r'^immune system process$',4),(r'regulation of immune (response|system process)',4),
(r'^T cell (activation|differentiation)$',4),(r'^T cell receptor signaling',3),(r'regulation of T cell activation',3),(r'T cell receptor signaling',3),
(r'intracellular iron ion homeostasis',4),
(r'^plasma membrane repair$',3),
# 3: plausible downstream/supporting
(r'virus|viral|interferon-alpha|interferon-beta|type I interferon',2),
(r'Gram-negative|defense response to bacterium',2),(r'lipopolysaccharide|toll-like receptor',2),
(r'cytokine|interleukin|chemokine|tumor necrosis factor|type II interferon',3),
(r'NF-kappaB',3),(r'NLRP3',3),(r'JAK-STAT|via STAT',3),(r'pattern recognition receptor',2),
(r'leukocyte cell-cell adhesion',3),(r'^phagocytosis',3),(r'regulation of phagocytosis',3),
(r'^chemotaxis$|cell chemotaxis|positive chemotaxis',3),(r'superoxide anion generation',3),
(r'keratinocyte differentiation',3),(r'^epidermis development$',3),(r'^skin development$',3),
(r'collagen fibril',3),(r'^extracellular matrix organization$',3),(r'basement membrane organization',2),
(r'positive regulation of fibroblast proliferation',3),
(r'^angiogenesis$|^sprouting angiogenesis$|regulation of angiogenesis|positive regulation of angiogenesis',3),
(r'vascular endothelial growth factor',2),(r'establishment of endothelial barrier',3),
(r'regulation of blood pressure',3),
(r'response to hypoxia|cellular response to hypoxia',3),
(r'oxidative stress|oxidant detoxification|hydrogen peroxide|reactive oxygen species|cell redox homeostasis',3),
(r'response to starvation|cellular response to starvation|response to nutrient levels|cellular response to nutrient levels|energy homeostasis',3),
(r'insulin-like growth factor receptor',3),
(r'nuclear receptor-mediated steroid hormone',3),
(r'sensory perception of pain',3),
(r'epithelial to mesenchymal transition$',2),(r'regulation of epithelial to mesenchymal',2),
(r'positive regulation of epithelial cell (proliferation|migration)',2),
(r'^cell migration$|regulation of cell migration$|positive regulation of cell migration$',2),
(r'^bone mineralization$|regulation of bone mineralization|^ossification$|osteoclast differentiation|osteoblast differentiation|^bone development$',2),
(r'^response to stress$',2),(r'^cellular response to stress$',2),
(r'^response to hormone$|cellular response to hormone stimulus|hormone-mediated signaling',2),
(r'^glucose homeostasis$|regulation of glucose metabolic|insulin receptor signaling pathway$|^response to insulin$|cellular response to insulin',2),
(r'^lipid storage$|triglyceride homeostasis|fatty acid beta-oxidation|lipid catabolic process',2),
(r'^heme|regulation of heart rate$',2),
(r'cold-induced thermogenesis',2),
(r'^learning$|^memory$|^cognition$',2),
(r'^regulation of cell growth$|positive regulation of cell growth',2),
(r'^cell adhesion mediated by integrin|integrin-mediated signaling',2),
(r'focal adhesion assembly|cell-matrix adhesion',2),
(r'transforming growth factor beta receptor signaling pathway$|cellular response to transforming growth factor',2),
(r'epidermal growth factor receptor signaling pathway$|cellular response to epidermal growth factor',2),
(r'fibroblast growth factor receptor signaling pathway$',2),
(r'^apoptotic process$',2),
(r'retinol metabolic',2),(r'circadian rhythm$',1),
(r'regulated exocytosis',2),(r'^response to xenobiotic stimulus$|cellular response to xenobiotic',2),
(r'p38MAPK',2),(r'response to cytokine',3),
(r'^cell population proliferation$',1),
(r'endothelial cell (proliferation|migration)|blood vessel endothelial cell migration|cell migration involved in sprouting',2),
(r'^blood vessel development$',2),
(r'intracellular zinc ion homeostasis|^zinc ion transport',2),
(r'lysosom',1),
(r'mechanical stimulus',2),
(r'^membrane protein ectodomain proteolysis',2),
(r'^proteolysis$',1),
(r'protein secretion$',1),
(r'cellular response to amino acid stimulus|amino acid starvation|glucose starvation',2),
(r'gluconeogenesis|glycogen biosynthetic',2),
(r'positive regulation of cell-matrix adhesion|positive regulation of cell adhesion',1),
(r'^skeletal muscle tissue development$',1),
(r'regulation of reactive|regulation of macroautophagy',1),
]
rows=[]
lines=open('/home/claude/runs/run_028/input_terms.txt').read().splitlines()[1:]
for l in lines:
    if not l.strip(): continue
    gid,name=l.split(' ',1)
    s=1
    for p,v in rules:
        if re.search(p,name): s=v;break
    rows.append((gid,name,s))
with open('/home/claude/runs/run_028/scores.tsv','w') as f:
    f.write('pathway\tscore\n')
    for g,n,s in rows: f.write(f'{g}\t{s}\n')
import sys
if len(sys.argv)>1:
    for g,n,s in sorted(rows,key=lambda r:-r[2]):
        if s>1: print(s,g,n)
print(len(rows))
