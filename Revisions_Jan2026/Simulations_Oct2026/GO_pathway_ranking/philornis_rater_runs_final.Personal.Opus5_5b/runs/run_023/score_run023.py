import re
rows = [l.rstrip('\n') for l in open('/home/claude/runs/run_023/input_terms.txt')][1:]
rows = [r for r in rows if r.strip()]

# ordered rules: first match wins
rules = [
 # 5: directly documented mechanisms
 (5, r'^(erythrocyte differentiation|positive regulation of erythrocyte differentiation|hemopoiesis|B cell activation|B cell receptor signaling pathway|adaptive immune response|inflammatory response|positive regulation of inflammatory response|regulation of inflammatory response|wound healing|response to wounding|blood coagulation|immune response|defense response)$'),
 # 4: core components / strong predictions
 (4, r'^(B cell differentiation|positive regulation of B cell proliferation|T cell activation|positive regulation of T cell activation|T cell receptor signaling pathway|T cell differentiation|innate immune response|activation of innate immune response|positive regulation of innate immune response|immune system process|regulation of immune response|regulation of immune system process|hematopoietic progenitor cell differentiation|platelet activation|platelet aggregation|positive regulation of wound healing|intracellular iron ion homeostasis|regulation of cytokine production|positive regulation of cytokine production|cytokine-mediated signaling pathway|response to cytokine|cellular response to cytokine stimulus|positive regulation of interleukin-6 production|complement activation|leukocyte cell-cell adhesion|positive regulation of chemokine production|chemotaxis|cell chemotaxis|positive regulation of tumor necrosis factor production|tumor necrosis factor-mediated signaling pathway|cellular response to interleukin-1|positive regulation of interleukin-8 production|interleukin-17-mediated signaling pathway|canonical NF-kappaB signal transduction|positive regulation of canonical NF-kappaB signal transduction|phagocytosis|negative regulation of inflammatory response|plasma membrane repair|response to tumor necrosis factor|cellular response to tumor necrosis factor|skin development|epidermis development|keratinocyte differentiation)$'),
 # 3: plausible downstream/supporting
 (3, r'^(regulation of canonical NF-kappaB signal transduction|negative regulation of canonical NF-kappaB signal transduction|non-canonical NF-kappaB signal transduction|positive regulation of non-canonical NF-kappaB signal transduction|negative regulation of non-canonical NF-kappaB signal transduction|negative regulation of cytokine production|negative regulation of interleukin-6 production|negative regulation of tumor necrosis factor production|positive regulation of interleukin-12 production|positive regulation of type II interferon production|cellular response to type II interferon|negative regulation of innate immune response|negative regulation of T cell activation|negative regulation of T cell receptor signaling pathway|toll-like receptor signaling pathway|toll-like receptor 4 signaling pathway|positive regulation of NLRP3 inflammasome complex assembly|positive regulation of phagocytosis|phagocytosis, engulfment|cell surface receptor signaling pathway via JAK-STAT|positive regulation of receptor signaling pathway via JAK-STAT|negative regulation of receptor signaling pathway via JAK-STAT|cell surface receptor signaling pathway via STAT|response to lipopolysaccharide|cellular response to lipopolysaccharide|defense response to bacterium|defense response to Gram-negative bacterium|superoxide anion generation|response to oxidative stress|cellular response to oxidative stress|cellular oxidant detoxification|cellular response to hydrogen peroxide|cellular response to reactive oxygen species|cell redox homeostasis|response to hypoxia|cellular response to hypoxia|angiogenesis|sprouting angiogenesis|positive regulation of angiogenesis|regulation of angiogenesis|collagen fibril organization|extracellular matrix organization|epithelial cell differentiation|positive regulation of epithelial cell proliferation|positive regulation of epithelial cell migration|positive regulation of fibroblast proliferation|cell migration|regulation of cell migration|positive regulation of cell migration|transforming growth factor beta receptor signaling pathway|cellular response to transforming growth factor beta stimulus|establishment of endothelial barrier|regulation of blood pressure|glucose homeostasis|response to starvation|cellular response to starvation|response to nutrient levels|cellular response to nutrient levels|energy homeostasis|bone development|ossification|bone mineralization|regulation of bone mineralization|positive regulation of bone mineralization|cartilage development|chondrocyte differentiation|osteoblast differentiation|skeletal system development|embryonic skeletal system morphogenesis|embryonic skeletal system development|osteoclast differentiation|positive regulation of osteoclast differentiation|roof of mouth development|leukocyte|response to stress|cellular response to stress|response to hormone|hormone-mediated signaling pathway|nuclear receptor-mediated steroid hormone signaling pathway|sensory perception of pain|response to xenobiotic stimulus|cellular response to xenobiotic stimulus|xenobiotic transport|cellular response to epidermal growth factor stimulus|epidermal growth factor receptor signaling pathway|insulin-like growth factor receptor signaling pathway|regulation of cell growth|positive regulation of cell growth|cellular response to growth factor stimulus|cellular response to vascular endothelial growth factor stimulus|vascular endothelial growth factor receptor signaling pathway|vascular endothelial growth factor signaling pathway|positive regulation of vascular endothelial growth factor production|cell-matrix adhesion|regulation of reactive oxygen species metabolic process|positive regulation of reactive oxygen species metabolic process|p38|positive regulation of p38MAPK cascade|regulation of heart rate|response to amino acid|negative regulation of ferroptosis|lipid storage|glycogen biosynthetic process|fatty acid beta-oxidation|gluconeogenesis|protein catabolic process|blood vessel development)$'),
 # 2: weak/indirect
 (2, r'(MAPK|ERK1 and ERK2|JNK cascade|apoptotic process|apoptotic signaling|interferon|virus|viral|pattern recognition|^regulation of cell population proliferation|^positive regulation of cell population proliferation|^cell population proliferation|^cell adhesion$|regulation of cell adhesion|integrin|endothelial cell|smooth muscle cell proliferation|BMP|SMAD|insulin|glucose|^lipid metabolic|fatty acid beta|fatty acid metabolic|triglyceride|lipid homeostasis|lipid storage|lipid catabolic|^macroautophagy$|^autophagy$|unfolded protein|endoplasmic reticulum stress|response to heat|^cellular response to heat|zinc ion homeostasis|^calcium ion homeostasis|focal adhesion|^epithelial|stem cell|morphogenesis of an epithelium|^response to hormone|hormone stimulus|hormone-mediated|retinol|thermogenesis|sensory perception of pain|senescence|^DNA damage response$|fat cell|tissue development|anatomical structure morphogenesis|^animal organ development|^cell motility|phosphatidylinositol 3-kinase|^TORC1 signaling|^TOR signaling|aerobic respiration|glycolytic|ATP metabolic|tricarboxylic|amino acid metabolic|purinergic|response to ATP|developmental process|epithelial to mesenchymal|^proteolysis$|ceramide|ferroptosis|^response to stress|cellular response to stress|growth factor stimulus|carbohydrate metabolic|xenobiotic|starvation|nutrient|cytosolic calcium|wound|cell migration|cell growth|cholesterol homeostasis|glycogen|gluconeogenesis|lipid biosynthetic|^fatty acid|^Wnt signaling|canonical Wnt|^Notch signaling|regulation of Notch|^hippo|cell differentiation$|limb|bone|ossification|osteo|chondro|cartilage|skeletal system|keratin|epiderm|skin|collagen|extracellular matrix|basement membrane|response to estrogen|estradiol|androgen|steroid|progesterone|calcium-mediated signaling|cellular response to calcium|circadian rhythm$|energy|heart rate|blood pressure|hypoxia|iron|heme|reactive oxygen|oxidative|redox|oxidant|hydrogen peroxide|superoxide|lysosom|endocytosis$|receptor-mediated endocytosis|phagocyt|exosomal|protein secretion|exocytosis$|regulated exocytosis|^vesicle-mediated transport$|muscle organ|skeletal muscle tissue|myoblast|angiogen|blood vessel|vascular)'),
]

out = []
for r in rows:
    gid, name = r.split(' ', 1)
    s = 1
    for score, pat in rules:
        if re.search(pat, name):
            s = score
            break
    out.append((gid, name, s))

# manual overrides after review
over = {
 'GO:0160240': 1,
 'GO:0006508': 2, 'GO:0030163': 2,
 'GO:0002376': 4, 'GO:0006955': 5,
 'GO:0006952': 5,
 'GO:0007596': 4,  # host coagulation vs hematophagy: strong prediction, not documented in host
 'GO:0042060': 4, 'GO:0009611': 5,
 'GO:0090303': 4,
 'GO:0046718': 1,
 'GO:0043923': 1, 'GO:0045070': 1,
 'GO:0006954': 5, 'GO:0050729': 5, 'GO:0050727': 4,
 'GO:0050728': 4,
 'GO:0002250': 5, 'GO:0042113': 5, 'GO:0050853': 5, 'GO:0030183': 4, 'GO:0030890': 4,
 'GO:0030218': 5, 'GO:0045648': 5, 'GO:0030097': 5, 'GO:0002244': 4,
 'GO:0006879': 4,
 'GO:0001666': 3, 'GO:0071456': 3,
 'GO:0042593': 3, 'GO:0006006': 2, 'GO:0010906': 2, 'GO:0001678': 2,
 'GO:0050776': 4, 'GO:0002682': 4,
 'GO:0042981': 2, 'GO:0006915': 2,
 'GO:0033554': 2, 'GO:0006950': 3,
 'GO:0009410': 2, 'GO:0071466': 2, 'GO:0042908': 2,
 'GO:0019233': 3,
 'GO:0008217': 3,
 'GO:0001501': 3, 'GO:0060021': 3, 'GO:0048704': 2, 'GO:0048706': 2,
 'GO:0030316': 2, 'GO:0045672': 2,
 'GO:0016477': 3, 'GO:0030334': 2, 'GO:0030335': 2,
 'GO:0009725': 2, 'GO:0009755': 2,
 'GO:0030518': 3,  # steroid hormone (corticosterone) signaling
 'GO:0032868': 2, 'GO:0032869': 2,
 'GO:0097009': 3,
 'GO:0005978': 2, 'GO:0006635': 2, 'GO:0006094': 2, 'GO:0019915': 2,
 'GO:0110076': 2,
 'GO:1900745': 3,
 'GO:0002027': 2,
 'GO:0071230': 2,
 'GO:0048009': 3,
 'GO:0007173': 3, 'GO:0071364': 3,
 'GO:0001568': 3,
 'GO:0007160': 3,
 'GO:0061028': 3,
 'GO:0098869': 3,
 'GO:0050870': 4, 'GO:0042110': 4,
 'GO:0050852': 4, 'GO:0030217': 4,
 'GO:0070527': 4, 'GO:0030168': 4,
 'GO:0007159': 4,
 'GO:0001778': 3,
 'GO:0043588': 4, 'GO:0008544': 4, 'GO:0030216': 4, 'GO:0045618': 3,
 'GO:0006956': 4,
 'GO:0006909': 4, 'GO:0006911': 3, 'GO:0050766': 3,
 'GO:0060326': 4, 'GO:0006935': 3, 'GO:0050918': 3, 'GO:0050919': 2,
 'GO:0032722': 4,
 'GO:0032755': 4, 'GO:0032757': 4, 'GO:0032760': 4,
 'GO:0071347': 4, 'GO:0097400': 4,
 'GO:0033209': 3, 'GO:0034612': 3, 'GO:0071356': 3,
 'GO:0019221': 4, 'GO:0034097': 4, 'GO:0071345': 4,
 'GO:0001817': 4, 'GO:0001819': 4, 'GO:0001818': 3,
 'GO:0007249': 3, 'GO:0043123': 3, 'GO:0043122': 3,
 'GO:0045087': 4, 'GO:0002218': 4, 'GO:0045089': 4, 'GO:0045824': 3,
 'GO:0002224': 3, 'GO:0034142': 3,
 'GO:0051607': 2, 'GO:0009615': 2, 'GO:0098586': 2, 'GO:0140374': 2,
 'GO:0002230': 2, 'GO:0002753': 2,
 'GO:0032728': 2, 'GO:0032727': 2, 'GO:0032481': 2, 'GO:0032480': 2, 'GO:0060339': 2,
 'GO:0032729': 3, 'GO:0071346': 3, 'GO:0032735': 3,
 'GO:0042742': 3, 'GO:0050829': 3, 'GO:0032496': 3, 'GO:0071222': 3,
 'GO:0042554': 3,
 'GO:0001525': 3, 'GO:0002040': 2, 'GO:0045766': 3, 'GO:0045765': 3,
 'GO:0030199': 3, 'GO:0030198': 3,
 'GO:0048146': 3,
 'GO:0007179': 3, 'GO:0071560': 3,
 'GO:0042594': 3, 'GO:0009267': 3, 'GO:0031667': 3, 'GO:0031669': 3,
 'GO:0060348': 3, 'GO:0001503': 3, 'GO:0030282': 2, 'GO:0030500': 2, 'GO:0030501': 2,
 'GO:0051216': 3, 'GO:0002062': 2, 'GO:0001649': 2,
}

for g in ['GO:0001755','GO:2000300','GO:0048488','GO:0043524','GO:0043525','GO:0051402','GO:0060317','GO:0035987','GO:0034198','GO:0016226','GO:0050919','GO:0042771','GO:0008630','GO:0032570','GO:0071230','GO:0019827','GO:0048863']:
    over[g]=1
over.update({'GO:0007204':2,'GO:0006897':2,'GO:0006898':2,'GO:0002230':2,'GO:0045669':2,'GO:0071711':2})

names2 = ['cellular response to mechanical stimulus','cell-cell junction assembly','bicellular tight junction assembly','positive regulation of cell-matrix adhesion','cell-substrate adhesion','substrate adhesion-dependent cell spreading','positive regulation of substrate adhesion-dependent cell spreading','lamellipodium assembly','membrane protein ectodomain proteolysis','intermediate filament organization','sensory perception of sound','regulation of cell-cell adhesion','cell-cell adhesion','negative regulation of epithelial cell proliferation','positive regulation of autophagy','negative regulation of autophagy','regulation of autophagy','regulation of macroautophagy','positive regulation of macroautophagy','negative regulation of transforming growth factor beta receptor signaling pathway','positive regulation of transforming growth factor beta receptor signaling pathway','adherens junction organization','regulation of intracellular pH','protein folding','positive regulation of proteolysis']
for i, (g, n, s) in enumerate(out):
    if n in names2: out[i]=(g,n,2); continue
    if g in over:
        out[i] = (g, n, over[g])

with open('/home/claude/runs/run_023/scores.tsv', 'w') as f:
    f.write('pathway\tscore\n')
    for g, n, s in out:
        f.write(f'{g}\t{s}\n')

import sys
if len(sys.argv) > 1:
    for g, n, s in out:
        if s >= int(sys.argv[1]):
            print(s, g, n)
print(len(out), {k: sum(1 for x in out if x[2] == k) for k in range(1, 6)})
