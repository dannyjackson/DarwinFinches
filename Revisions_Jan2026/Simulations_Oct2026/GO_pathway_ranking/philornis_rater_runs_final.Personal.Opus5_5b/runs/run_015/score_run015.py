import re
terms=[l.rstrip('\n').split(' ',1) for l in open('/home/claude/runs/run_015/input_terms.txt')][1:]
terms=[t for t in terms if t and t[0]]
S5=['GO:0042060','GO:0009611','GO:0002250','GO:0006954','GO:0007596','GO:0030218','GO:0006955','GO:0050853','GO:0042113','GO:0090303','GO:0006952','GO:0045648']
rules=[
(5,r'^(wound healing|response to wounding|inflammatory response|blood coagulation|immune response|adaptive immune response|defense response|B cell receptor signaling|B cell activation|erythrocyte differentiation|positive regulation of erythrocyte differentiation|positive regulation of wound healing)$'),
(4,r'regulation of inflammatory response|positive regulation of inflammatory|negative regulation of inflammatory|B cell (proliferation|differentiation)|complement activation|platelet (aggregation|activation)|hemopoiesis|innate immune response$|regulation of immune response|T cell (activation|receptor signaling|differentiation)|positive regulation of T cell|negative regulation of T cell|interleukin-(6|8|1\b|17|12)|tumor necrosis factor|cytokine production|cytokine-mediated|response to cytokine|cellular response to cytokine|chemokine|^chemotaxis|cell chemotaxis|leukocyte cell-cell adhesion|hematopoietic progenitor|intracellular iron ion homeostasis|keratinocyte differentiation|positive regulation of keratinocyte|epidermis development|skin development|immune system process|regulation of immune system process|sensory perception of pain|phagocytosis|type II interferon|NF-kappaB|toll-like receptor|activation of innate immune|regulation of innate immune|NLRP3|defense response to bacterium|Gram-negative|response to lipopolysaccharide|cellular response to lipopolysaccharide|angiogenesis$|sprouting angiogenesis$|regulation of angiogenesis|positive regulation of angiogenesis|collagen fibril|response to oxidative stress|cellular response to oxidative stress|cellular oxidant detoxification|superoxide anion generation|response to hypoxia|cellular response to hypoxia|interleukin-1$'),
(3,r'JAK-STAT|via STAT|cellular response to interleukin-1|MAPK|p38|JNK|ERK1|response to tumor|interferon|pattern recognition|extracellular matrix organization|fibroblast proliferation|TGF|transforming growth factor|epithelial cell (proliferation|migration)|positive regulation of epithelial cell migration|endothelial|vascular endothelial growth factor|blood vessel development|basement membrane|cell migration$|regulation of cell migration|positive regulation of cell migration|cell motility|reactive oxygen|hydrogen peroxide|redox homeostasis|response to stress$|cellular response to stress|response to hormone|hormone-mediated|steroid hormone|cellular response to hormone|regulation of blood pressure|regulation of heart rate$|glucose homeostasis|energy homeostasis|response to starvation|cellular response to starvation|response to nutrient|cellular response to nutrient|lipid storage|fatty acid beta-oxidation|gluconeogenesis|insulin-like growth factor|insulin receptor signaling$|response to insulin|cellular response to insulin|growth$|regulation of cell growth|positive regulation of cell growth|apoptotic process$|regulation of apoptotic|cell population proliferation$|hemostasis|response to xenobiotic|xenobiotic transport|cellular response to xenobiotic|establishment of endothelial barrier|tight junction|plasma membrane repair|epidermal growth factor|EGF|wound|leukocyte|immune|antiviral|defense response to virus|response to virus|cellular response to virus|lysosom|autophagy|intracellular zinc|zinc ion trans|intracellular glucose|positive regulation of cytokine|negative regulation of cytokine|regulation of cytokine|bone mineralization$|ossification$|bone development|skeletal system development|osteoblast differentiation$|chondrocyte differentiation$|cartilage development|cellular response to mechanical|response to heat|cellular response to heat|ferroptosis|cell adhesion mediated by integrin|integrin-mediated|focal adhesion assembly$|response to wound|smooth muscle cell proliferation|erythrocyte|iron-sulfur|regulation of pH|bicarbonate'),
(2,r'^cell adhesion$|^cell-cell adhesion$|regulation of cell adhesion|regulation of cell-cell adhesion|positive regulation of cell adhesion|negative regulation of cell adhesion|cell-matrix adhesion|^lipid metabolic|lipid catabolic|fatty acid metabolic|fatty acid beta|lipid homeostasis|triglyceride|glucose metabolic|glycolytic|carbohydrate metabolic|glycogen|aerobic respiration|tricarboxylic|ATP metabolic|electron transport|respiratory chain|DNA damage response$|DNA repair$|apoptotic|proliferation$|regulation of cell population proliferation|endoplasmic reticulum stress|unfolded protein|protein folding|refolding|misfolded|senescence|intracellular calcium ion homeostasis|^calcium ion homeostasis|thermogenesis|insulin|estrogen|estradiol|androgen|progesterone|retinoic acid|retinol|cardiac muscle contraction|regulation of heart rate|muscle contraction$|stress granule|\bTOR signaling|TORC1|BMP signaling|cellular response to BMP|canonical Wnt|^Wnt signaling|Notch signaling pathway$|SMAD|fibroblast growth factor|platelet|^proteolysis$|positive regulation of proteolysis|ectodomain|protein secretion|^exocytosis$|regulated exocytosis|learning|memory|cognition|circadian rhythm|neuropeptide|regulation of multicellular|ionizing|UV|amino acid metabolic|amino acid starvation|glucose starvation|glucose stimulus|D-glucose import|epithelial cell differentiation|morphogenesis of an epithelium|tissue development|skeletal system|bone|ossification|osteoblast|osteoclast|chondrocyte|cartilage|adherens junction|cell-cell junction|apical/basal|keratin|intermediate filament|filopodium|lamellipodium|stress fiber|actin filament polymerization|Arp2/3|calcium-mediated signaling|cytosolic calcium|mitochondrial membrane potential|cytochrome c|mitophagy|ferroptosis|xenobiotic|response to calcium|purinergic|response to ATP|heparan sulfate|chondroitin|proteoglycan biosynthetic|cell growth|fat cell differentiation|response to growth factor|epithelial to mesenchymal|ubiquinone|iron-sulfur|stem cell population|hematopo|lung development|blood pressure|myoblast|skeletal muscle tissue|muscle organ|sarcomere|heat|phospholipid scrambling|mitochondrion organization|regulation of mitochondrion|mitochondrial fission|mitochondrial fusion|cholesterol|regulation of lipid metabolic|lipid storage|lipid biosynthetic|fatty acid biosynthetic|fatty acid transport|long-chain fatty acid|receptor-mediated endocytosis|^endocytosis$|ceramide|sphingolipid|limb|Rho protein signal|regulation of pH|intracellular pH|bicarbonate|zinc|PI3|phosphatidylinositol 3-kinase|TGF|activin|hippo|semaphorin|ephrin|embryo development ending|gluconeogenesis|smooth muscle|cellular response to amino acid|L-leucine|response to estrogen|nutrient|vascular|kidney|reactive|p53|DNA damage checkpoint|signal transduction in response to DNA damage|cell migration|chemotaxis|cell motility|autophagy|lysosom|autophagosome|phagocytosis|macroautophagy|epidermal')
]
out=[]
for gid,name in terms:
    s=1
    for sc,pat in rules:
        if re.search(pat,name,re.I): s=sc;break
    out.append([gid,name,s])
ov={
# manual adjustments
'GO:0050776':4,'GO:0045824':4,'GO:0045087':4,'GO:0002224':3,'GO:0034142':3,'GO:0007249':3,'GO:0038061':3,'GO:0043124':3,'GO:0043122':3,'GO:1901224':2,'GO:1901223':2,'GO:0043123':3,
'GO:0098586':2,'GO:0140374':2,'GO:0051607':2,'GO:0009615':2,'GO:0002230':2,'GO:0045070':1,'GO:0043923':1,'GO:0046718':1,
'GO:0032728':2,'GO:0060339':2,'GO:0032480':2,'GO:0032481':2,'GO:0032727':2,'GO:0071346':3,'GO:0032729':3,
'GO:0042742':4,'GO:0050829':3,'GO:0032496':3,'GO:0071222':3,
'GO:0006909':4,'GO:0006911':3,'GO:0050766':3,
'GO:0002218':4,'GO:0045089':4,'GO:0002753':3,'GO:1900227':3,
'GO:0002376':4,'GO:0002682':4,
'GO:0019233':3,'GO:0008544':4,'GO:0043588':4,'GO:0030216':4,'GO:0045618':4,
'GO:0006879':4,'GO:0002244':4,'GO:0030097':4,
'GO:0001525':4,'GO:0002040':3,'GO:0045765':3,'GO:0045766':4,'GO:0016525':3,
'GO:0030199':4,'GO:0030198':3,'GO:0071711':3,
'GO:0006979':4,'GO:0034599':4,'GO:0098869':4,'GO:0042554':4,'GO:0001666':4,'GO:0071456':4,'GO:0070301':3,'GO:0034614':3,'GO:0045454':3,'GO:2000379':3,'GO:2000377':3,
'GO:0070527':4,'GO:0030168':4,'GO:0006956':4,
'GO:0034612':4,'GO:0071356':4,'GO:0033209':4,'GO:0032760':4,'GO:0032720':4,
'GO:0032755':4,'GO:0032715':4,'GO:0032757':4,'GO:0097400':4,'GO:0032735':3,'GO:0071347':4,
'GO:0001817':4,'GO:0001819':4,'GO:0001818':4,'GO:0019221':4,'GO:0034097':4,'GO:0071345':4,'GO:0032722':4,
'GO:0060326':4,'GO:0006935':4,'GO:0050918':3,'GO:0050919':2,'GO:0007159':4,
'GO:0042110':4,'GO:0050870':4,'GO:0050868':4,'GO:0050852':4,'GO:0030217':4,'GO:0050860':4,
'GO:0030890':4,'GO:0030183':4,
'GO:0007259':3,'GO:0046427':3,'GO:0046426':3,'GO:0097696':3,
'GO:0006950':3,'GO:0033554':3,
'GO:0045742':3,'GO:0007173':3,'GO:0071364':3,'GO:0042059':3,
'GO:0010575':3,'GO:0048010':3,'GO:0038084':3,'GO:0035924':3,'GO:0001568':3,'GO:0001938':3,'GO:0001937':2,'GO:0043542':3,'GO:0010595':3,'GO:0043536':3,'GO:0002042':3,'GO:0090050':3,'GO:0090051':3,'GO:0061028':3,'GO:2000352':2,'GO:2000353':2,
'GO:0048146':3,'GO:0050679':3,'GO:0050680':3,'GO:0010634':3,'GO:0001778':3,'GO:0070830':3,
'GO:0071560':3,'GO:0007179':3,'GO:0030511':3,'GO:0030512':3,'GO:0001837':2,'GO:0010718':2,'GO:0010719':2,
'GO:0016477':3,'GO:0030334':3,'GO:0030335':3,'GO:0048870':2,'GO:0030336':2,
'GO:0009725':3,'GO:0009755':3,'GO:0030518':3,'GO:0032870':3,'GO:0030522':2,
'GO:0042593':3,'GO:0097009':3,'GO:0042594':3,'GO:0009267':3,'GO:0031667':3,'GO:0031669':3,'GO:0042149':2,'GO:0034198':2,'GO:0019915':3,'GO:0006635':2,'GO:0006094':3,
'GO:0048009':3,'GO:0008286':2,'GO:0032868':2,'GO:0032869':2,'GO:0046627':2,
'GO:0001558':3,'GO:0030307':3,'GO:0030308':2,
'GO:0006915':3,'GO:0042981':3,'GO:0043065':2,'GO:0043066':2,'GO:0008283':2,'GO:0042127':2,
'GO:0009410':2,'GO:0071466':2,'GO:0042908':2,
'GO:0008217':3,'GO:0002027':2,'GO:0086091':1,
'GO:0010506':2,'GO:0006914':2,'GO:0016236':2,'GO:0016239':2,'GO:0010508':2,'GO:0010507':2,'GO:0016241':2,
'GO:0006882':3,'GO:0071577':2,'GO:0006829':2,
'GO:0001678':2,'GO:0030282':2,'GO:0001503':2,'GO:0060348':2,'GO:0001501':2,'GO:0001649':2,'GO:0002062':2,'GO:0051216':2,'GO:0030500':2,'GO:0030501':2,
'GO:0071260':3,'GO:0009408':2,'GO:0034605':2,'GO:0110076':2,
'GO:0033627':3,'GO:0007229':3,'GO:0048041':2,
'GO:0048661':2,'GO:1904707':2,
'GO:0016226':2,'GO:0006885':2,'GO:0051453':1,'GO:0015701':1,
'GO:0007042':2,'GO:0032418':2,'GO:0007040':2,'GO:0008333':2,'GO:0007041':2,'GO:0006622':2,'GO:0061909':2,'GO:0097401':1,
'GO:0006509':2,'GO:0009306':2,'GO:0050714':2,
'GO:0007623':2,'GO:0042752':2,'GO:0032922':1,'GO:0043153':1,
'GO:0007612':2,'GO:0007613':2,'GO:0050890':2,
'GO:0007218':2,'GO:0120162':2,'GO:0120163':2,
'GO:0006954':5,'GO:0050853':4,'GO:0042113':5,'GO:0007596':5,'GO:0042060':5,'GO:0009611':5,'GO:0002250':5,'GO:0006955':5,'GO:0006952':4,'GO:0030218':5,'GO:0045648':4,'GO:0090303':4,'GO:0050727':4,'GO:0050729':4,'GO:0050728':4,
'GO:0034976':2,'GO:0030968':2,'GO:0070059':2,
'GO:0090398':2,'GO:2000773':2,
'GO:0045444':2,'GO:0045600':2,'GO:0045599':2,
'GO:0048731':1,'GO:0032502':1,'GO:0048513':1,'GO:0009888':1,'GO:0048468':1,'GO:0030154':1,'GO:0009653':1,'GO:0048856':1,
'GO:0160108':1,'GO:0050793':1,'GO:0051094':1,'GO:0045595':1,'GO:0045596':1,'GO:0045597':1,'GO:0009887':1,'GO:0048568':1,
'GO:0051239':2,'GO:0051240':2,
'GO:0007154':1,'GO:0023052':1,'GO:0023051':1,'GO:0007166':1,'GO:0009966':1,'GO:0009967':1,'GO:0009968':1,'GO:0035556':1,'GO:1902531':1,'GO:1902533':1,'GO:0007267':1,'GO:0007165':1,
'GO:0033198':2,'GO:0035589':2,
'GO:0060021':2,'GO:0048704':2,'GO:0048706':2,
'GO:0071363':2,
}
for r in out:
    if r[0] in ov: r[2]=ov[r[0]]
with open('/home/claude/runs/run_015/scores.tsv','w') as f:
    f.write('pathway\tscore\n')
    for g,n,s in out: f.write(f'{g}\t{s}\n')
import sys
if len(sys.argv)>1:
    for g,n,s in out:
        if s>=int(sys.argv[1]): print(s,g,n)
print(len(out))
