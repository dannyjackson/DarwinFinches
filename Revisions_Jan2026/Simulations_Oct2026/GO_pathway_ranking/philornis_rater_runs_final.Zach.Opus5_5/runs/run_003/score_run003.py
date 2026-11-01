import re
rows = []
with open('/home/claude/runs/run_003/input_terms.txt') as f:
    next(f)
    for line in f:
        line = line.rstrip('\n')
        if not line.strip():
            continue
        gid, _, name = line.partition(' ')
        rows.append((gid, name))

# (regex, score) rules; the maximum matching score is used
rules = [
    # 5: direct documented mechanisms (antibody/adaptive immunity, inflammation, wounding/blood loss, erythropoiesis, stress hormone)
    (r'^adaptive immune response$|^immune response$|^inflammatory response$|^response to wounding$|^wound healing$|^blood coagulation$|^hemopoiesis$|^erythrocyte differentiation$|^B cell activation$|^defense response$|^immune system process$', 5),
    (r'positive regulation of erythrocyte differentiation|positive regulation of wound healing|B cell receptor signaling|positive regulation of B cell proliferation|B cell differentiation|^regulation of immune response$|^regulation of inflammatory response$|positive regulation of inflammatory response|negative regulation of inflammatory response|complement activation|platelet aggregation|platelet activation|innate immune response$|regulation of immune system process|intracellular iron ion homeostasis', 4),
    (r'hematopoietic progenitor|T cell (activation|differentiation|receptor)|regulation of T cell activation|interleukin|tumor necrosis factor|cytokine|chemokine|NF-kappaB|toll-like receptor|pattern recognition|inflammasome|phagocytosis|leukocyte|chemotaxis|innate immune|epidermis development|skin development|keratinocyte|collagen fibril|extracellular matrix organization|angiogenesis|response to hypoxia|cellular response to hypoxia|oxidative stress|reactive oxygen species|oxidant detoxification|superoxide|hydrogen peroxide|response to stress$|type II interferon|JAK-STAT|via STAT|interleukin-17|lipopolysaccharide|plasma membrane repair|regulation of blood pressure|hormone-mediated signaling|response to hormone|cellular response to hormone|nuclear receptor-mediated steroid hormone|response to starvation|cellular response to starvation|energy homeostasis|glucose homeostasis|blood vessel development|endothelial barrier|response to nutrient levels|cellular response to nutrient levels', 3),
    (r'apoptotic|apoptosis|autophagy|TGF|transforming growth factor|cell migration|endothelial|wound|response to virus|defense response to (virus|bacterium|Gram)|antiviral|interferon|MAPK|ERK1|JNK|p38|insulin|gluconeogenesis|glycogen|lipid storage|triglyceride|fatty acid beta-oxidation|fatty acid metabolic|lipid catabolic|glucose metabolic|carbohydrate metabolic|growth|ossification|bone|osteo|chondrocyte|cartilage|skeletal system|epithelial cell proliferation|fibroblast|cell population proliferation|epithelial to mesenchymal|cell-matrix adhesion|integrin|focal adhesion|cellular response to stress|ER stress|endoplasmic reticulum stress|unfolded protein|heat|cell redox|ferroptosis|xenobiotic|VEGF|vascular endothelial growth factor|heart rate|cardiac muscle contraction|response to calcium|calcium-mediated|exocytosis|regulated exocytosis|protein secretion|receptor-mediated endocytosis|response to estrogen|estradiol|progesterone|androgen|retinol|retinoic|circadian|sensory perception of pain|learning|memory|^cognition$|neuropeptide|lysosom|proteolysis|cell adhesion|cell-cell adhesion|tight junction|thermogenesis|TOR signaling|TORC1|aerobic respiration|tricarboxylic|glycolytic|mitochondrial electron transport|ATP metabolic|amino acid starvation|glucose starvation|intracellular pH|regulation of pH|zinc ion|response to ATP|purinergic|mechanical stimulus|cellular senescence|epidermal growth factor|FGF|fibroblast growth factor|BMP|smooth muscle|myoblast|muscle', 2),
]

manual = {
    'GO:0060317': 1, 'GO:0060045': 1, 'GO:2000300': 1, 'GO:0060307': 1, 'GO:0086005': 1, 'GO:0001755': 1, 'GO:0008045': 1, 'GO:0048704': 1, 'GO:0048706': 1, 'GO:0099645': 1, 'GO:0099149': 1, 'GO:1904862': 1, 'GO:0060079': 1, 'GO:0099072': 1, 'GO:0006376': 1,
    'GO:0030097': 5, 'GO:0030218': 5, 'GO:0042060': 5, 'GO:0009611': 5, 'GO:0007596': 5,
    'GO:0002250': 5, 'GO:0006955': 5, 'GO:0006954': 5, 'GO:0042113': 5,
    'GO:0001666': 4, 'GO:0071456': 4,  # hypoxia from anemia
    'GO:0006879': 4,                   # iron homeostasis (blood loss)
    'GO:0090303': 4, 'GO:0045648': 4, 'GO:0050853': 4, 'GO:0030890': 4, 'GO:0030183': 4,
    'GO:0006956': 4, 'GO:0070527': 4, 'GO:0030168': 4,
    'GO:0006950': 3, 'GO:0006952': 5, 'GO:0002376': 5,
    'GO:0045087': 4, 'GO:0050776': 4, 'GO:0050727': 4, 'GO:0050729': 4, 'GO:0050728': 4,
    'GO:0002682': 4,
    'GO:0002244': 3, 'GO:0001525': 3, 'GO:0045766': 3, 'GO:0001568': 3,
    'GO:0008544': 4, 'GO:0043588': 4, 'GO:0030216': 3, 'GO:0045618': 3,
    'GO:0030198': 3, 'GO:0030199': 3,
    'GO:0008217': 3,
    'GO:0009725': 3, 'GO:0032870': 3, 'GO:0030518': 3, 'GO:0009755': 3,
    'GO:0042593': 3, 'GO:0097009': 3, 'GO:0042594': 3, 'GO:0009267': 3,
    'GO:0031667': 3, 'GO:0031669': 3,
    'GO:0006979': 3, 'GO:0034599': 3, 'GO:0098869': 3, 'GO:0045454': 3,
    'GO:0042554': 3, 'GO:0034614': 3, 'GO:0070301': 3, 'GO:2000377': 2, 'GO:2000379': 2,
    'GO:0001778': 3, 'GO:0061028': 3,
    'GO:0019233': 3,  # pain perception: response to larval feeding, possible behavioral defense
    'GO:0060326': 3, 'GO:0006935': 2, 'GO:0050918': 2, 'GO:0050919': 2,
    'GO:0006909': 3, 'GO:0006911': 3, 'GO:0050766': 3,
    'GO:0007159': 3, 'GO:0034113': 2,
    'GO:0042110': 3, 'GO:0030217': 3, 'GO:0050852': 3, 'GO:0050870': 3, 'GO:0050868': 3, 'GO:0050860': 3,
    # antiviral/antibacterial: weak
    'GO:0009615': 2, 'GO:0098586': 2, 'GO:0051607': 2, 'GO:0140374': 2, 'GO:0002230': 2,
    'GO:0042742': 2, 'GO:0050829': 2, 'GO:0043923': 1, 'GO:0045070': 1, 'GO:0046718': 1,
    'GO:0032728': 2, 'GO:0032727': 2, 'GO:0032480': 2, 'GO:0032481': 2, 'GO:0060339': 2,
    'GO:0071346': 3, 'GO:0032729': 3,
    'GO:0034142': 2, 'GO:0071222': 2, 'GO:0032496': 2, 'GO:0002224': 3, 'GO:0002753': 2,
    'GO:1900227': 3,
    # growth / bone: beak deformation, nestling growth
    'GO:0001503': 2, 'GO:0060348': 2, 'GO:0001501': 2, 'GO:0030282': 2,
    # misc
    'GO:0160240': 1, 'GO:0022008': 1, 'GO:0007399': 1,
    'GO:0043627': 2, 'GO:0071391': 2, 'GO:0071392': 2, 'GO:0032570': 2, 'GO:0030521': 2,
    'GO:0007605': 1, 'GO:0007601': 1,
    'GO:0035924': 2, 'GO:0038084': 2, 'GO:0048010': 2, 'GO:0010575': 2,
    'GO:0090398': 1, 'GO:2000773': 1,
    'GO:0034097': 3, 'GO:0071345': 3, 'GO:0019221': 3,
    'GO:0034612': 3, 'GO:0071356': 3, 'GO:0033209': 3,
    'GO:0071347': 3, 'GO:0097400': 3,
    'GO:0045824': 3, 'GO:0045089': 3, 'GO:0002218': 3,
    'GO:0001817': 3, 'GO:0001818': 3, 'GO:0001819': 3,
    'GO:0007249': 3, 'GO:0043122': 3, 'GO:0043123': 3, 'GO:0043124': 3,
    'GO:0038061': 2, 'GO:1901224': 2, 'GO:1901223': 2,
    'GO:0007259': 3, 'GO:0046426': 2, 'GO:0046427': 2, 'GO:0097696': 2,
    'GO:0009408': 2, 'GO:0034605': 2,
    'GO:0001937': 2, 'GO:0001938': 2,
    'GO:0007612': 2, 'GO:0007613': 2, 'GO:0050890': 2,
    'GO:0002027': 2, 'GO:0086091': 2,
    'GO:0120162': 2, 'GO:0120163': 2,
    'GO:0042908': 1, 'GO:0009410': 2, 'GO:0071466': 2,
}

def score(gid, name):
    if gid in manual:
        return manual[gid]
    s = 1
    for pat, sc in rules:
        if re.search(pat, name):
            s = max(s, sc)
    return s

out = [(g, score(g, n), n) for g, n in rows]
with open('/home/claude/runs/run_003/scores.tsv', 'w') as f:
    f.write('pathway\tscore\n')
    for g, s, _ in out:
        f.write(f'{g}\t{s}\n')
import sys
if len(sys.argv) > 1:
    for g, s, n in sorted(out, key=lambda x: -x[1]):
        if s >= int(sys.argv[1]):
            print(s, g, n)
print(len(out), len(set(g for g, _ in rows)))
