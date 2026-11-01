import re
src = "/home/claude/runs/run_002/input_terms.txt"
out = "/home/claude/runs/run_002/scores.tsv"

# Explicit term-level judgments (by exact term name). Anything not listed gets
# a pattern-based score below, which was reviewed term-by-term.
S5 = {
 "immune response", "adaptive immune response", "B cell activation", "B cell receptor signaling pathway",
 "B cell differentiation", "positive regulation of B cell proliferation", "inflammatory response",
 "wound healing", "response to wounding", "hemopoiesis", "erythrocyte differentiation",
 "positive regulation of erythrocyte differentiation", "blood coagulation", "defense response",
 "innate immune response", "immune system process",
}
S4 = {
 "regulation of immune system process", "regulation of immune response", "cytokine-mediated signaling pathway",
 "positive regulation of cytokine production", "regulation of cytokine production", "negative regulation of cytokine production",
 "response to cytokine", "cellular response to cytokine stimulus", "positive regulation of inflammatory response",
 "negative regulation of inflammatory response", "regulation of inflammatory response", "complement activation",
 "T cell activation", "T cell receptor signaling pathway", "T cell differentiation", "positive regulation of T cell activation",
 "negative regulation of T cell activation", "negative regulation of T cell receptor signaling pathway",
 "positive regulation of interleukin-6 production", "negative regulation of interleukin-6 production",
 "positive regulation of tumor necrosis factor production", "negative regulation of tumor necrosis factor production",
 "tumor necrosis factor-mediated signaling pathway", "response to tumor necrosis factor", "cellular response to tumor necrosis factor",
 "positive regulation of interleukin-8 production", "positive regulation of interleukin-12 production",
 "positive regulation of type II interferon production", "cellular response to type II interferon",
 "positive regulation of chemokine production", "interleukin-17-mediated signaling pathway", "cellular response to interleukin-1",
 "canonical NF-kappaB signal transduction", "positive regulation of canonical NF-kappaB signal transduction",
 "regulation of canonical NF-kappaB signal transduction", "negative regulation of canonical NF-kappaB signal transduction",
 "activation of innate immune response", "positive regulation of innate immune response", "negative regulation of innate immune response",
 "positive regulation of NLRP3 inflammasome complex assembly", "leukocyte cell-cell adhesion", "platelet activation",
 "platelet aggregation", "positive regulation of wound healing", "phagocytosis", "phagocytosis, engulfment",
 "positive regulation of phagocytosis", "intracellular iron ion homeostasis", "hematopoietic progenitor cell differentiation",
 "response to hypoxia", "cellular response to hypoxia", "skin development", "epidermis development",
 "keratinocyte differentiation", "positive regulation of keratinocyte differentiation", "cell chemotaxis",
 "superoxide anion generation", "response to stress",
}
S3 = {
 "toll-like receptor signaling pathway", "toll-like receptor 4 signaling pathway", "cytoplasmic pattern recognition receptor signaling pathway",
 "non-canonical NF-kappaB signal transduction", "positive regulation of non-canonical NF-kappaB signal transduction",
 "negative regulation of non-canonical NF-kappaB signal transduction", "cell surface receptor signaling pathway via JAK-STAT",
 "positive regulation of receptor signaling pathway via JAK-STAT", "negative regulation of receptor signaling pathway via JAK-STAT",
 "cell surface receptor signaling pathway via STAT", "chemotaxis", "positive chemotaxis", "response to lipopolysaccharide",
 "cellular response to lipopolysaccharide", "defense response to bacterium", "defense response to Gram-negative bacterium",
 "angiogenesis", "sprouting angiogenesis", "positive regulation of angiogenesis", "regulation of angiogenesis",
 "negative regulation of angiogenesis", "collagen fibril organization", "extracellular matrix organization",
 "response to oxidative stress", "cellular response to oxidative stress", "cellular oxidant detoxification",
 "cellular response to reactive oxygen species", "cellular response to hydrogen peroxide",
 "positive regulation of reactive oxygen species metabolic process", "regulation of reactive oxygen species metabolic process",
 "response to starvation", "cellular response to starvation", "response to nutrient levels", "cellular response to nutrient levels",
 "glucose homeostasis", "glucose metabolic process", "energy homeostasis", "intracellular glucose homeostasis",
 "cellular response to glucose starvation", "cellular response to amino acid starvation", "gluconeogenesis",
 "positive regulation of epithelial cell proliferation", "positive regulation of epithelial cell migration",
 "positive regulation of fibroblast proliferation", "epithelial to mesenchymal transition", "plasma membrane repair",
 "transforming growth factor beta receptor signaling pathway", "cellular response to transforming growth factor beta stimulus",
 "positive regulation of transforming growth factor beta receptor signaling pathway",
 "negative regulation of transforming growth factor beta receptor signaling pathway",
 "regulation of blood pressure", "regulation of heart rate", "sensory perception of pain",
 "positive regulation of type I interferon production", "negative regulation of type I interferon production",
 "cell migration", "leukocyte", "establishment of endothelial barrier", "basement membrane organization",
 "cellular response to mechanical stimulus", "iron-sulfur cluster assembly", "positive regulation of protein secretion",
 "regulated exocytosis", "protein secretion", "intermediate filament organization", "response to hormone",
 "cellular response to hormone stimulus", "hormone-mediated signaling pathway", "nuclear-transcribed mRNA catabolic process" if False else "_",
 "cellular response to vascular endothelial growth factor stimulus", "vascular endothelial growth factor signaling pathway",
 "vascular endothelial growth factor receptor signaling pathway", "positive regulation of vascular endothelial growth factor production",
 "cellular response to epidermal growth factor stimulus", "epidermal growth factor receptor signaling pathway",
 "regulation of cell migration", "positive regulation of cell migration", "endothelial cell migration",
 "positive regulation of endothelial cell migration", "positive regulation of blood vessel endothelial cell migration",
 "apoptotic process", "lipid storage", "triglyceride homeostasis", "response to insulin", "insulin receptor signaling pathway",
 "cellular response to insulin stimulus", "insulin secretion", "regulation of insulin secretion", "positive regulation of insulin secretion",
 "negative regulation of insulin secretion", "positive regulation of D-glucose import", "regulation of glucose metabolic process",
 "cellular response to glucose stimulus", "positive regulation of insulin secretion involved in cellular response to glucose stimulus",
 "response to ATP", "G protein-coupled purinergic nucleotide receptor signaling pathway", "osteoclast differentiation",
 "positive regulation of osteoclast differentiation", "fibroblast growth factor receptor signaling pathway",
 "positive regulation of p38MAPK cascade", "nucleotide" if False else "_",
 "response to xenobiotic stimulus", "cellular response to xenobiotic stimulus", "proteolysis", "membrane protein ectodomain proteolysis",
 "fatty acid beta-oxidation", "lipid catabolic process", "blood vessel development", "cellular response to heat", "response to heat",
 "cholesterol homeostasis", "circadian rhythm", "learning", "memory", "cognition", "neuropeptide signaling pathway",
 "positive regulation of cell adhesion", "cell adhesion mediated by integrin", "integrin-mediated signaling pathway",
 "heme" if False else "_", "stem cell differentiation", "protein N-linked glycosylation" if False else "_",
 "focal adhesion assembly", "positive regulation of focal adhesion assembly", "skeletal system development", "bone development",
 "ossification", "bone mineralization", "positive regulation of bone mineralization", "regulation of bone mineralization",
 "osteoblast differentiation",
 "glycogen biosynthetic process", "response to endoplasmic reticulum stress", "endoplasmic reticulum unfolded protein response",
 "cellular response to stress", "positive regulation of cell population proliferation", "cell population proliferation",
 "regulation of cell population proliferation", "negative regulation of ferroptosis", "positive regulation of exosomal secretion",
 "receptor-mediated endocytosis", "antiviral innate immune response",
}
S1_patterns = [
 r"virus|viral", r"meio", r"sperm", r"oogenesis", r"germ cell", r"gonad", r"rRNA", r"tRNA", r"snRNA", r"snRNP", r"piRNA", r"miRNA",
 r"splic", r"RNA polymerase", r"transcription", r"chromatin", r"chromosome", r"sister chromatid", r"kinetochore", r"spindle",
 r"centrosome|centriole", r"cilium|ciliary|axonem|dynein|cilia", r"synap", r"neurotransmitter", r"axon", r"dendrite", r"neuron",
 r"neural|neurogenesis", r"brain|cerebral|retina|lens|eye|photoreceptor|inner ear|sound|visual", r"heart|cardiac|ventricular|valve|outflow",
 r"kidney|ureteric|lung|roof of mouth|limb|digit|somito|pattern|left/right|dorsal/ventral|endoderm",
 r"DNA (replication|repair|recombination|biosynth|metabolic)|double-strand|mismatch|excision|cross-link|replication fork|telomere|DNA-templated DNA",
 r"ribosom", r"mRNA", r"nucle(ar|olus|osome|ocytoplasmic)", r"Golgi", r"peroxisom", r"ubiquit", r"sumoyl", r"proteasom",
 r"mitotic|cell cycle|G1/S|G2/M|cytokinesis|midbody|cell division|nuclear division",
]

def score(term):
    if term in S5: return 5
    if term in S4: return 4
    if term in S3: return 3
    t = term
    for p in S1_patterns:
        if re.search(p, t):
            return 1
    # generic default for remaining terms: weak/indirect (2) vs none (1)
    weak2 = [r"apopto", r"autophag", r"oxidat|redox|reactive oxygen", r"mitochond.*(respir|electron)|aerobic respiration|tricarboxylic|glycoly|ATP metabolic",
             r"lipid|fatty acid|triglyceride|cholesterol", r"glucose|carbohydrate|glycogen", r"amino acid|leucine|starvation|TOR",
             r"cell adhesion|cell-cell adhesion|tight junction|adherens|cadherin", r"actin|stress fiber|filopod|lamellipod|cell migration|motility",
             r"endothelial|vascular|smooth muscle|angiogen", r"growth factor|BMP|SMAD|activin|Wnt|Notch|hippo|smoothened",
             r"MAPK|ERK|JNK|p38|Ras|Rho|Rap|PI3|phosphatidylinositol 3-kinase|calcium|NF-kappa|STAT",
             r"hormone|estrogen|estradiol|progesterone|androgen|steroid|retino",
             r"fibroblast|epithelial|ossification|bone|osteo|cartilage|chondro|skeletal|collagen|keratin|extracellular matrix",
             r"muscle", r"exocytosis|secretion|endocyt|lysosom|vesicle", r"proliferation", r"senescence", r"thermogenesis",
             r"circadian", r"hypoxia", r"xenobiotic", r"heat|UV|ionizing", r"iron|zinc", r"cytokine|interferon|interleukin|tumor necrosis",
             r"pH|bicarbonate|proton", r"protein folding|refolding|misfolded|ERAD", r"receptor", r"stress",
             r"response to|cellular response", r"fat cell|adipo", r"insulin", r"pain"]
    for p in weak2:
        if re.search(p, t):
            return 2
    return 1

rows = []
with open(src) as f:
    next(f)
    for line in f:
        line = line.rstrip("\n")
        if not line.strip(): continue
        gid, name = line.split(" ", 1)
        rows.append((gid, name, score(name)))

# Manual overrides after line-by-line review
over = {
 "GO:0006954": 5, "GO:0042060": 5, "GO:0009611": 5, "GO:0030097": 5, "GO:0030218": 5, "GO:0045648": 5,
 "GO:0007596": 5, "GO:0002250": 5, "GO:0042113": 5, "GO:0050853": 5, "GO:0006955": 5,
 "GO:0034614": 3, "GO:0070301": 3, "GO:0006979": 3, "GO:0034599": 3, "GO:0098869": 3,
 "GO:0001666": 4, "GO:0071456": 4, "GO:0006879": 4,
 "GO:0043923": 1, "GO:0045070": 1, "GO:0046718": 1, "GO:0098586": 1, "GO:0009615": 1, "GO:0051607": 1,
 "GO:0002230": 1, "GO:0140374": 2, "GO:1904646": 1, "GO:0160240": 1,
 "GO:0042554": 3, "GO:0050766": 4, "GO:0007159": 4,
 "GO:0007160": 2, "GO:0031589": 2,
 "GO:0001503": 2, "GO:0060348": 2, "GO:0001501": 2, "GO:0030282": 2, "GO:0030501": 2, "GO:0030500": 2, "GO:0001649": 2,
 "GO:0007612": 2, "GO:0007613": 2, "GO:0050890": 2, "GO:0007623": 2, "GO:0042752": 2, "GO:0032922": 1, "GO:0043153": 1,
 "GO:0034605": 2, "GO:0009408": 2, "GO:0042632": 2, "GO:0006508": 2, "GO:0006509": 3,
 "GO:0019233": 3, "GO:0008217": 3, "GO:0002027": 2, "GO:0086091": 1, "GO:0071260": 2,
 "GO:0008283": 2, "GO:0042127": 2, "GO:0008284": 2, "GO:0033554": 3, "GO:0006950": 4,
 "GO:1903543": 2, "GO:0006898": 2, "GO:0045055": 2, "GO:0009306": 2, "GO:0050714": 3,
 "GO:0110076": 2, "GO:0016226": 2, "GO:0030316": 2, "GO:0045672": 2, "GO:0019915": 2, "GO:0070328": 2,
 "GO:0045109": 3, "GO:0048863": 2, "GO:0006915": 2, "GO:0005978": 2,
 "GO:0030968": 2, "GO:0034976": 2, "GO:0033198": 2, "GO:0035589": 2,
 "GO:0009410": 2, "GO:0071466": 2, "GO:0016042": 2, "GO:0006635": 2,
 "GO:0001568": 3, "GO:0060326": 4, "GO:0006935": 3, "GO:0050918": 3, "GO:0050919": 2,
 "GO:0009725": 3, "GO:0032870": 3, "GO:0009755": 3, "GO:0007218": 2,
 "GO:0030855": 2, "GO:0008544": 4, "GO:0043588": 4,
 "GO:0002244": 4, "GO:0019827": 2,
 "GO:0007179": 3, "GO:0071560": 3, "GO:0030511": 3, "GO:0030512": 3,
 "GO:0030334": 2, "GO:0030335": 2, "GO:0016477": 2, "GO:0048870": 2, "GO:0043542": 2, "GO:0010595": 2, "GO:0043536": 2,
 "GO:0010634": 3, "GO:0050679": 3, "GO:0048146": 3, "GO:0001837": 2, "GO:0001778": 3,
 "GO:0045785": 2, "GO:0033627": 2, "GO:0007229": 2, "GO:0048041": 2, "GO:0051894": 2,
 "GO:0034097": 4, "GO:0071345": 4, "GO:0019221": 4,
 "GO:0045087": 5, "GO:0006952": 5, "GO:0002376": 5,
 "GO:0042742": 3, "GO:0050829": 3, "GO:0032496": 3, "GO:0071222": 3,
 "GO:0042593": 3, "GO:0006006": 3, "GO:0097009": 3, "GO:0001678": 3, "GO:0042149": 3, "GO:0034198": 2, "GO:0006094": 3,
 "GO:0042594": 3, "GO:0009267": 3, "GO:0031667": 3, "GO:0031669": 3,
 "GO:0032868": 2, "GO:0008286": 2, "GO:0032869": 2, "GO:0030073": 2, "GO:0050796": 2, "GO:0032024": 2, "GO:0046676": 2,
 "GO:0046326": 2, "GO:0010906": 2, "GO:0071333": 2, "GO:0035774": 2,
 "GO:0002224": 3, "GO:0034142": 3, "GO:0002753": 3, "GO:0038061": 3, "GO:1901224": 3, "GO:1901223": 3,
 "GO:0007259": 3, "GO:0046427": 3, "GO:0046426": 3, "GO:0097696": 3,
 "GO:0001525": 3, "GO:0002040": 3, "GO:0045766": 3, "GO:0045765": 3, "GO:0016525": 3,
 "GO:0030199": 3, "GO:0030198": 3, "GO:0071711": 2, "GO:0061028": 2,
 "GO:2000379": 3, "GO:2000377": 3,
 "GO:0035924": 2, "GO:0038084": 2, "GO:0048010": 2, "GO:0010575": 2, "GO:0071364": 2, "GO:0007173": 2, "GO:0008543": 2,
 "GO:1900745": 3, "GO:0032481": 2, "GO:0032480": 2, "GO:0032727": 1, "GO:0032728": 1, "GO:0060339": 1,
 "GO:0030168": 4, "GO:0070527": 4, "GO:0090303": 4,
 "GO:0002682": 4, "GO:0050776": 4, "GO:0001817": 4, "GO:0001819": 4, "GO:0001818": 4,
 "GO:0050727": 4, "GO:0050728": 4, "GO:0050729": 4, "GO:0006956": 4,
 "GO:0042110": 4, "GO:0050852": 4, "GO:0030217": 4, "GO:0050870": 4, "GO:0050868": 4, "GO:0050860": 4,
 "GO:0030183": 5, "GO:0030890": 5,
 "GO:0032755": 4, "GO:0032715": 4, "GO:0032760": 4, "GO:0032720": 4, "GO:0033209": 4, "GO:0034612": 4, "GO:0071356": 4,
 "GO:0032757": 4, "GO:0032735": 4, "GO:0032729": 4, "GO:0071346": 4, "GO:0032722": 4, "GO:0097400": 4, "GO:0071347": 4,
 "GO:0007249": 4, "GO:0043123": 4, "GO:0043122": 4, "GO:0043124": 4,
 "GO:0002218": 4, "GO:0045089": 4, "GO:0045824": 4, "GO:1900227": 4,
 "GO:0006909": 4, "GO:0006911": 4, "GO:0016239": 2,
 "GO:0045618": 4, "GO:0030216": 4, "GO:0019216": 2, "GO:0006629": 2,
 "GO:0006885": 2, "GO:0051453": 2,
 "GO:0001558": 2, "GO:0030307": 2, "GO:0030308": 2, "GO:0030518": 2, "GO:0006509": 2,
 "GO:0045109": 2, "GO:0050714": 2, "GO:0006950": 3,
}
for i, (g, n, s) in enumerate(rows):
    if g in over:
        rows[i] = (g, n, over[g])

with open(out, "w") as f:
    f.write("pathway\tscore\n")
    for g, n, s in rows:
        f.write(f"{g}\t{s}\n")
print(len(rows))
import collections
print(collections.Counter(s for _,_,s in rows))
for k in (5,4,3,2):
    print("==", k, [n for _,n,s in rows if s==k])
