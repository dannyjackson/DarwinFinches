import re
inp = "/home/claude/runs/run_019/input_terms.txt"
out = "/home/claude/runs/run_019/scores.tsv"

# Ordered rules: first match wins. (regex on lowercase term name, score)
rules = [
 # 5: directly documented mechanisms (IgY antibody/adaptive & immune response, inflammation, blood loss/anemia, wounds, coagulation)
 (r"^immune response$|^adaptive immune response$|^b cell activation$|b cell receptor signaling|b cell differentiation|positive regulation of b cell proliferation", 5),
 (r"^inflammatory response$|^wound healing$|^response to wounding$|^positive regulation of wound healing$|^blood coagulation$|^erythrocyte differentiation$|positive regulation of erythrocyte differentiation|^hemopoiesis$|^defense response$|^immune system process$", 5),
 # 4: core components
 (r"regulation of inflammatory response|regulation of immune response|regulation of immune system process|innate immune response|activation of innate immune response|t cell activation|t cell receptor signaling|t cell differentiation|platelet activation|platelet aggregation|cytokine production|cytokine-mediated|response to cytokine|cellular response to cytokine|interleukin|tumor necrosis factor|complement activation|leukocyte cell-cell adhesion|chemotaxis|cell chemotaxis|phagocytosis|hematopoietic progenitor|intracellular iron ion homeostasis|skin development|epidermis development|keratinocyte differentiation|nf-kappab|chemokine|type ii interferon|plasma membrane repair|angiogenesis$|^sprouting angiogenesis|nlrp3", 4),
 # 3: plausible downstream/supporting
 (r"jak-stat|via stat|toll-like|pattern recognition|lipopolysaccharide|defense response to|response to virus|antiviral|interferon|immune|superoxide|oxidative stress|reactive oxygen|oxidant detox|hydrogen peroxide|cell redox|collagen fibril|extracellular matrix organization|basement membrane|epithelial cell proliferation|epithelial cell migration|fibroblast proliferation|establishment of endothelial barrier|vascular endothelial growth factor|blood vessel|angiogenesis|endothelial cell (migration|proliferation)|blood vessel endothelial|regulation of blood pressure|response to hypoxia|cellular response to hypoxia|glucose homeostasis|energy homeostasis|response to starvation|response to nutrient|cellular response to nutrient|cellular response to starvation|glucose starvation|amino acid starvation|insulin-like growth factor|growth hormone|^response to stress$|cellular response to stress|apoptotic process$|^apoptotic signaling|extrinsic apoptotic|ossification|bone mineralization|osteoblast differentiation|osteoclast differentiation|bone development|cartilage|chondrocyte|tgf|transforming growth factor|smad|p38mapk|jnk cascade|heme|iron|response to hormone|hormone-mediated|cellular response to hormone|nuclear receptor-mediated steroid|sensory perception of pain|intermediate filament organization|epithelial to mesenchymal transition$|negative regulation of epithelial to mesenchymal|positive regulation of epithelial to mesenchymal|cell migration$|regulation of cell migration|positive regulation of cell migration|integrin|cell-matrix adhesion|focal adhesion|wound|xenobiotic|circadian rhythm$|learning|memory|cognition|fat cell differentiation|lipid storage|triglyceride|fatty acid beta-oxidation|gluconeogenesis|glycogen biosynthetic|cellular response to insulin|response to insulin|insulin receptor|insulin secretion|heat|muscle organ development|skeletal muscle tissue development|positive regulation of cell growth|regulation of cell growth|negative regulation of cell growth|lysosom|autophagy|macroautophagy|calcium-mediated signaling|cell surface receptor signaling pathway via stat|epithelium|epithelial cell differentiation|tissue development|cellular senescence|regulation of reactive", 3),
]

# explicit overrides (score) by GO ID for nuance
over = {
 "GO:0032496": 3, "GO:0071222": 3,
 "GO:0046718": 1, "GO:0045070": 1, "GO:0043923": 1,  # viral host processes
 "GO:0002230": 2, "GO:0140374": 2, "GO:0051607": 2, "GO:0009615": 2, "GO:0098586": 2,
 "GO:0032727": 2, "GO:0032728": 2, "GO:0032481": 2, "GO:0032480": 2, "GO:0060339": 2,
 "GO:0050829": 3, "GO:0042742": 3,
 "GO:0001503": 2, "GO:0030282": 2, "GO:0030500": 2, "GO:0030501": 2, "GO:0060348": 2,
 "GO:0045669": 2, "GO:0045668": 2, "GO:0001649": 2, "GO:0030279": 2,
 "GO:0045672": 2, "GO:0030316": 2,  # beak remodeling speculative
 "GO:0051216": 2, "GO:0002062": 2, "GO:0032331": 2,
 "GO:0009408": 2, "GO:0034605": 2, "GO:0034644": 1, "GO:0009411": 1,
 "GO:0007612": 2, "GO:0007613": 2, "GO:0050890": 2,
 "GO:0045599": 2, "GO:0045600": 2, "GO:0045444": 2,
 "GO:0006096": 2, "GO:0005978": 2, "GO:0006094": 3,
 "GO:0030073": 2, "GO:0050796": 2, "GO:0032024": 2, "GO:0046676": 2, "GO:0035774": 2,
 "GO:0008286": 3, "GO:0046627": 2,
 "GO:1904646": 1,
 "GO:0090398": 2, "GO:2000773": 2,
 "GO:0006915": 3, "GO:0097190": 2, "GO:0097191": 3,
 "GO:0007519": 3, "GO:0007517": 3,
 "GO:0045109": 3,
 "GO:0001837": 3, "GO:0060317": 1,
 "GO:0003341": 1,
 "GO:0002009": 2, "GO:0009888": 2,
 "GO:0007040": 2, "GO:0007042": 2, "GO:0032418": 2, "GO:0006622": 2, "GO:0008333": 2, "GO:0007041": 2,
 "GO:0006914": 2, "GO:0016236": 2, "GO:0010506": 2, "GO:0010508": 2, "GO:0010507": 2, "GO:0016239": 2, "GO:0016241": 2,
 "GO:0000045": 2, "GO:0097352": 2, "GO:0061909": 2,
 "GO:0019722": 2, "GO:0050850": 2,
 "GO:0030518": 3, "GO:0009725": 3, "GO:0009755": 3, "GO:0032870": 3,
 "GO:0019233": 3, "GO:0007623": 2,
 "GO:0006954": 5, "GO:0042060": 5,
 "GO:0016525": 3, "GO:0045765": 3, "GO:0045766": 3, "GO:0001525": 3, "GO:0002040": 3,
 "GO:0090050": 2, "GO:0090051": 2, "GO:0002042": 2,
 "GO:0001568": 3,
 "GO:0006979": 3, "GO:0034599": 3, "GO:0098869": 3,
 "GO:0042554": 3, "GO:0045454": 2, "GO:2000377": 2, "GO:2000379": 2,
 "GO:0001666": 3, "GO:0071456": 3,
 "GO:0006879": 4, "GO:0016226": 2,
 "GO:0031669": 3, "GO:0031667": 3, "GO:0042594": 3, "GO:0009267": 3, "GO:0042149": 2, "GO:0034198": 2,
 "GO:0097009": 3, "GO:0042593": 3, "GO:0001678": 2,
 "GO:0048009": 3,
 "GO:0001501": 2, "GO:0048706": 2, "GO:0048704": 2,
 "GO:0060021": 2,  # roof of mouth (beak/palate) weak
 "GO:0001778": 3,
 "GO:0043588": 4, "GO:0008544": 4, "GO:0030216": 3, "GO:0045618": 3,
 "GO:0030198": 3, "GO:0030199": 3, "GO:0071711": 2,
 "GO:0050679": 3, "GO:0010634": 3, "GO:0050680": 2, "GO:0030855": 2,
 "GO:0048146": 3,
 "GO:0007229": 2, "GO:0033627": 2, "GO:0007160": 2, "GO:0031589": 2, "GO:0001954": 2,
 "GO:0048041": 2, "GO:0051894": 2, "GO:0051895": 2,
 "GO:0016477": 2, "GO:0030334": 2, "GO:0030335": 2, "GO:0030336": 2, "GO:0048870": 2,
 "GO:0006935": 3, "GO:0060326": 3, "GO:0050918": 3, "GO:0050919": 2,
 "GO:0006909": 3, "GO:0006911": 3, "GO:0050766": 3,
 "GO:0002224": 3, "GO:0034142": 3, "GO:0002753": 3,
 "GO:0030168": 4, "GO:0070527": 4,
 "GO:0007596": 5,
 "GO:0030097": 4, "GO:0030218": 5, "GO:0045648": 5, "GO:0002244": 3,
 "GO:0008217": 3,
 "GO:0061028": 3,
 "GO:0006956": 4,
 "GO:0002250": 5, "GO:0042113": 5, "GO:0050853": 4, "GO:0030183": 4, "GO:0030890": 4,
 "GO:0006955": 5, "GO:0002376": 4, "GO:0006952": 4, "GO:0002682": 4, "GO:0050776": 4,
 "GO:0045087": 4, "GO:0002218": 4, "GO:0045089": 4, "GO:0045824": 4,
 "GO:0050729": 4, "GO:0050728": 4, "GO:0050727": 4,
 "GO:0042110": 4, "GO:0050852": 4, "GO:0050870": 4, "GO:0050868": 3, "GO:0030217": 4, "GO:0050860": 3,
 "GO:0007159": 3,
 "GO:0001817": 4, "GO:0001818": 4, "GO:0001819": 4,
 "GO:0019221": 4, "GO:0034097": 4, "GO:0071345": 4,
 "GO:0032755": 4, "GO:0032715": 3, "GO:0032757": 3, "GO:0032735": 3, "GO:0071347": 4,
 "GO:0097400": 3, "GO:0032729": 3, "GO:0071346": 3,
 "GO:0032760": 4, "GO:0032720": 3, "GO:0033209": 3, "GO:0034612": 3, "GO:0071356": 3,
 "GO:0032722": 3,
 "GO:0043122": 3, "GO:0043123": 3, "GO:0043124": 3, "GO:0007249": 3,
 "GO:0038061": 2, "GO:1901223": 2, "GO:1901224": 2,
 "GO:1900227": 3,
 "GO:0007259": 3, "GO:0097696": 3, "GO:0046426": 2, "GO:0046427": 3,
 "GO:0009611": 5, "GO:0090303": 4,
 "GO:0006950": 2, "GO:0033554": 2,
 "GO:0007179": 3, "GO:0071560": 3, "GO:0030511": 3, "GO:0030512": 3, "GO:0060395": 2, "GO:0060391": 2,
 "GO:1900745": 2, "GO:0007254": 2, "GO:0046330": 2, "GO:0046329": 2,
 "GO:0042981": 2, "GO:0043066": 2, "GO:0043065": 2, "GO:2001237": 2, "GO:2001238": 2,
 "GO:1902042": 2,
 "GO:0001819": 4,
 "GO:0009410": 2, "GO:0071466": 2, "GO:0042908": 2,
 "GO:0006629": 2, "GO:0019915": 2, "GO:0070328": 2, "GO:0019432": 2, "GO:0006635": 2,
 "GO:0006631": 2, "GO:0016042": 2,
 "GO:0032869": 2, "GO:0032868": 2,
 "GO:0030308": 3, "GO:0030307": 3, "GO:0001558": 3,
 "GO:0042026": 1,
 "GO:0031175": 1,
 "GO:0006882": 2, "GO:0071577": 2,
 "GO:0001657": 1,
 "GO:0035589": 2, "GO:0033198": 2,
 "GO:0001937": 2, "GO:0001938": 2, "GO:0043542": 2, "GO:0043536": 2, "GO:0010595": 2,
 "GO:2000352": 2, "GO:2000353": 2,
 "GO:0010575": 3, "GO:0048010": 2, "GO:0038084": 2, "GO:0035924": 2,
 "GO:0048661": 1, "GO:1904707": 1, "GO:0060045": 1,
 "GO:0120162": 2, "GO:0120163": 2,
 "GO:0007173": 2, "GO:0071364": 2, "GO:0045742": 2, "GO:0042059": 2,
 "GO:0008543": 2,
 "GO:0001843": 1,
 "GO:0030166": 1, "GO:0050650": 1, "GO:0015012": 1, "GO:0030200": 1,
 "GO:1903543": 1,
 "GO:0006954": 5,
 "GO:0002027": 2, "GO:0086091": 1, "GO:0060048": 1, "GO:0055117": 1,
 "GO:0042542": 3,
 "GO:0006508": 2, "GO:0045862": 2,
 "GO:0009311": 1,
 "GO:0006336": 1,
 "GO:0006520": 2, "GO:0003333": 1, "GO:0006865": 1,
 "GO:0006006": 2, "GO:0010906": 2, "GO:0071333": 2, "GO:0046326": 2,
 "GO:0006099": 1, "GO:0009060": 2, "GO:0046034": 2,
 "GO:0007218": 2,
 "GO:0042752": 1, "GO:0032922": 1, "GO:0043153": 1,
 "GO:0045785": 2, "GO:0030155": 2, "GO:0022407": 2, "GO:0007162": 2,
 "GO:0034113": 2,
 "GO:0006897": 1, "GO:0016192": 1, "GO:0006887": 2, "GO:0045055": 2, "GO:0009306": 2, "GO:0050714": 2,
 "GO:0051604": 1, "GO:0016485": 1,
 "GO:0030010": 1,
 "GO:0030036": 1,
 "GO:0007155": 2,
 "GO:0016055": 1, "GO:0060070": 1,
 "GO:0048731": 1,
 "GO:0051239": 1, "GO:0051240": 1,
 "GO:0031647": 1,
 "GO:0042542": 3,
 "GO:0042277": 1,
 "GO:0006801": 3,
 "GO:0042573": 1,
 "GO:0042572": 2,  # retinol, vitamin A immunity/epithelia weak
 "GO:0048384": 2, "GO:0032526": 2, "GO:0071300": 2,
 "GO:0070301": 2, "GO:0034614": 3, "GO:0302": 0,
 "GO:0001666": 3,
 "GO:0110076": 2,
 "GO:0090201": 1,
 "GO:0030521": 2,  # androgen, immunosuppression / parental behavior weak
 "GO:0032570": 1, "GO:0043627": 1, "GO:0071391": 1, "GO:0071392": 1,
 "GO:0007605": 1,
 "GO:0007601": 1, "GO:0019233": 3,
 "GO:0050877": 1,
 "GO:0060021": 2,
 "GO:0006813": 1,
 "GO:0001755": 1, "GO:0006376": 1, "GO:0043524": 1, "GO:0043525": 1, "GO:0051402": 1, "GO:0034727": 1,
 "GO:0070830": 2, "GO:0007043": 2, "GO:0030307": 2, "GO:0030308": 2, "GO:0001558": 2,
}

rows = []
with open(inp) as f:
    lines = f.read().splitlines()
for ln in lines[1:]:
    if not ln.strip():
        continue
    gid, name = ln.split(" ", 1)
    n = name.lower()
    s = None
    if gid in over and over[gid] > 0:
        s = over[gid]
    else:
        for rx, sc in rules:
            if re.search(rx, n):
                s = sc
                break
    if s is None:
        s = 1
    rows.append((gid, name, s))

with open(out, "w") as f:
    f.write("pathway\tscore\n")
    for g, _, s in rows:
        f.write(f"{g}\t{s}\n")

import sys
if len(sys.argv) > 1:
    for g, n, s in rows:
        if s >= 2:
            print(s, g, n)
print(len(rows))
