import re
inp = "/home/claude/runs/run_024/input_terms.txt"
out = "/home/claude/runs/run_024/scores.tsv"

# Explicit scores by GO ID (judgment), checked first
S = {}
def setv(score, ids):
    for i in ids.split():
        S[i] = score

# 5: directly documented mechanisms (antibody/IgY response, blood loss & anaemia/erythropoiesis,
# coagulation vs. haematophagy, wound/tissue damage, inflammation, immune response)
setv(5, """GO:0007596 GO:0042060 GO:0009611 GO:0002250 GO:0042113 GO:0030183 GO:0050853 GO:0030218
GO:0045648 GO:0006954 GO:0006955 GO:0006952 GO:0030097 GO:0090303 GO:0070527 GO:0030168""")
# 4: core components / strong predictions
setv(4, """GO:0006956 GO:0030890 GO:0050776 GO:0002682 GO:0045087 GO:0002376 GO:0050729 GO:0050727
GO:0050728 GO:0042110 GO:0050870 GO:0030217 GO:0050852 GO:0034612 GO:0071356 GO:0033209 GO:0071347
GO:0019221 GO:0034097 GO:0071345 GO:0001817 GO:0001819 GO:0001818 GO:0032755 GO:0032715 GO:0032760
GO:0032720 GO:0097400 GO:0032757 GO:0032722 GO:0006935 GO:0060326 GO:0002244 GO:0006879 GO:0042742
GO:0006909 GO:0006911 GO:0050766 GO:0045089 GO:0002218 GO:0045824 GO:0007159 GO:0008544 GO:0043588
GO:0030216 GO:0045618 GO:0001778 GO:0032735 GO:0032729 GO:0071346 GO:0034142 GO:0002224 GO:0032496
GO:0071222 GO:0002040 GO:0001525 GO:0050868 GO:0050860 GO:0007249 GO:0043122 GO:0043123 GO:0043124
GO:0050918 GO:0001503 GO:0060348 GO:0030199 GO:0071346 GO:0042554 GO:0019233 GO:0001666 GO:0071456
GO:0050829 GO:0030198 GO:0008217""")
# 3: plausible downstream/supporting
setv(3, """GO:0007259 GO:0046427 GO:0046426 GO:0097696 GO:0038061 GO:1901224 GO:1901223 GO:1900227
GO:0002753 GO:0006979 GO:0034599 GO:0098869 GO:0070301 GO:0034614 GO:2000377 GO:2000379 GO:0045454
GO:0006882 GO:0045766 GO:0016525 GO:0045765 GO:0038084 GO:0048010 GO:0035924 GO:0010575 GO:0001568
GO:0071560 GO:0007179 GO:0030511 GO:0030512 GO:0048146 GO:0050679 GO:0010634 GO:0030335 GO:0016477
GO:0048870 GO:0030334 GO:0000902 GO:0001501 GO:0030282 GO:0030500 GO:0030501 GO:0001649 GO:0030316
GO:0045672 GO:0051216 GO:0002062 GO:0042593 GO:0097009 GO:0019915 GO:0009725 GO:0032870 GO:0009755
GO:0030518 GO:0006629 GO:0009267 GO:0042594 GO:0031667 GO:0031669 GO:0008286 GO:0032868 GO:0032869
GO:0048009 GO:0006950 GO:0033554 GO:0001937 GO:0001938 GO:0043542 GO:0010595 GO:0002042 GO:0043536
GO:0090050 GO:0090051 GO:0061028 GO:0071711 GO:0007160 GO:0033627 GO:0007229 GO:0098609 GO:0007155
GO:0045785 GO:0030155 GO:0022407 GO:0007162 GO:0001954 GO:0034113 GO:0045597 GO:0002009 GO:0030855
GO:0070830 GO:0007043 GO:0034332 GO:0006915 GO:0042981 GO:0043066 GO:0043065 GO:0097191 GO:2001237
GO:2001238 GO:1902042 GO:0009615 GO:0051607 GO:0140374 GO:0032481 GO:0032480 GO:0032728 GO:0032727
GO:0060339 GO:0098586 GO:0002230 GO:0070588 GO:0001756 GO:0007218 GO:0007613 GO:0007612 GO:0050890
GO:0016042 GO:0006635 GO:0005975 GO:0006006 GO:0006094 GO:0001676 GO:0070328 GO:0055088 GO:0045444
GO:0009408 GO:0034605 GO:0035556 GO:0009410 GO:0071466 GO:0042908 GO:0035589 GO:0033198 GO:0048661
GO:1904707 GO:0010718 GO:0001837 GO:0071260 GO:0019216 GO:0006508 GO:0045862 GO:0009060 GO:0046034
GO:0006099 GO:0006096 GO:0007623 GO:0051209 GO:0019722 GO:0007204 GO:0006954 GO:0070374 GO:0000165
GO:0043410 GO:0007254 GO:1900745 GO:0046330 GO:0008284 GO:0008283 GO:0071363 GO:0007173 GO:0071364
GO:0008543 GO:0001508 GO:0042552 GO:0007605 GO:0050877 GO:0002027 GO:0006936 GO:0003009 GO:0007517
GO:0007519 GO:0045214 GO:0043627 GO:0032570 GO:0071392 GO:0071391 GO:0030521 GO:0110076 GO:0006897
GO:0006898 GO:0045055 GO:0009306 GO:0050714 GO:0030073 GO:0042572 GO:0032526 GO:0071300 GO:0048384
GO:0045669 GO:0045668 GO:0030279 GO:0016236 GO:0006914 GO:0007041 GO:0007040 GO:0051897 GO:0043491
GO:0031929 GO:0038202 GO:0045089 GO:0030166 GO:0050650 GO:0015012 GO:0030200 GO:0006487 GO:0009101
GO:0006493 GO:0046718 GO:0070059 GO:0034976 GO:0030968 GO:0036503 GO:0006457 GO:0006954""")
# 2: weak/indirect
setv(2, """GO:0000122 GO:0045893 GO:0045892""")  # placeholders overwritten below by default 2 logic


# final adjustments
for g in "GO:0001756 GO:0042552".split(): S[g]=1
for g in """GO:0007605 GO:0001508 GO:0006936 GO:0003009 GO:0007517 GO:0007519 GO:0045214 GO:0043627 GO:0032570
GO:0071392 GO:0071391 GO:0030521 GO:0007613 GO:0050890 GO:0007218 GO:0046718 GO:0070059 GO:0042572 GO:0032526
GO:0071300 GO:0048384 GO:0007169 GO:0019827 GO:0007267 GO:0009408 GO:0034605 GO:0007623 GO:0048661 GO:1904707
GO:0010718 GO:0001837 GO:0045444 GO:0030166 GO:0050650 GO:0015012 GO:0030200 GO:0006493 GO:0006487 GO:0009101""".split(): S[g]=-2

# keyword-based fallback for terms not explicitly scored
two_kw = ["signal transduction", "cascade", "tor signaling", "torc1", "insulin", "glucose", "lipid", "fatty acid",
          "cholesterol", "triglyceride", "iron", "zinc", "calcium", "oxidative", "redox", "hypoxia",
          "apoptotic", "autophag", "lysosom", "endocyt", "exocyt", "secretion", "ubiquitin",
          "proteasom", "protein catabolic", "mitochondri", "respiratory", "atp ", "heat", "stress",
          "migration", "adhesion", "angiogen", "endothelial", "vascular", "blood", "epithel",
          "keratinocyte", "skin", "bone", "ossification", "osteo", "cartilage", "chondro", "skeletal",
          "muscle", "heart rate", "growth", "proliferation", "starvation", "nutrient", "amino acid",
          "hormone", "estrogen", "estradiol", "androgen", "steroid", "retino", "wnt", "bmp", "transforming growth",
          "smad", "notch", "hippo", "fibroblast", "epidermal growth", "ras protein", "rho protein",
          "integrin", "actin", "chemotaxis", "membrane potential", "potassium", "sodium", "chloride",
          "regulation of ph", "intracellular ph", "circadian", "learning", "memory", "neuropeptide",
          "pain", "phagocyt", "virus", "viral", "interferon", "symbiont", "xenobiotic", "ferroptosis",
          "senescence", "extracellular matrix", "collagen", "basement", "glycosylation", "proteoglycan",
          "thermogenesis", "energy", "phosphatidylinositol 3-kinase", "camp", "purinergic", "fat cell",
          "cytokine", "nf-kappab", "jak-stat", "stat", "wound", "mechanical stimulus", "protein folding",
          "unfolded protein", "er stress", "endoplasmic reticulum stress", "erad", "refolding", "cell motility",
          "cell shape", "stress fiber", "focal adhesion", "lamellipodium", "filopodium", "tight junction",
          "cell-cell junction", "adherens", "metal ion", "monoatomic cation", "proteolysis", "cell population proliferation",
          "cell growth", "glycogen", "gluconeogenesis", "glycolytic", "carbohydrate", "tricarboxylic",
          "aerobic respiration", "electron transport", "response to calcium", "response to nutrient", "ceramide",
          "sphingolipid", "phospholipid metabolic", "receptor signaling pathway via"]

scores = []
with open(inp) as f:
    lines = f.read().splitlines()
for ln in lines[1:]:
    if not ln.strip():
        continue
    gid, _, name = ln.partition(" ")
    if gid in S and S[gid] == -2:
        s = 2
    elif gid in S and S[gid] != 2:
        s = S[gid]
    else:
        nm = name.lower()
        s = 1
        for k in two_kw:
            if k.lower() in nm:
                s = 2
                break
        # explicit weak-only categories pulled back to 1 despite keyword hits
        lows = ["synap", "axon", "dendrite", "neuron projection", "spermat", "sperm", "meiotic",
                "oogenesis", "cilium", "ciliary", "dynein", "neural tube", "lens", "retina",
                "ureteric", "kidney", "cardiac", "heart", "valve", "outflow", "limb", "digit",
                "roof of mouth", "somito", "left/right", "pattern", "neural crest", "lung",
                "inner ear", "eye", "gonad", "germ cell", "rRNA", "tRNA", "snRNA", "spliceosom",
                "splicing", "chromatid", "kinetochore", "spindle", "centro", "centriole",
                "nuclear membrane", "DNA replication", "golgi", "nucleus", "nucleolus",
                "endoderm", "photoreceptor", "neurotransmitter", "presynap", "postsynap",
                "GABA", "glutamat", "acetylcholine", "cholinergic", "amyloid", "piRNA", "miRNA",
                "peroxisom", "clathrin", "vesicle budding", "intraciliary", "mitochondrial translation",
                "ribosom", "magnesium", "sulfate", "bicarbonate", "neuromuscular", "myelin",
                "intra-golgi", "nuclear export", "nuclear import", "protein import into nuc",
                "mrna export", "phospholipid translocation", "scrambling", "neurogenesis",
                "neuron differentiation", "neural precursor", "brain", "cerebral", "central nervous",
                "nervous system development", "motor neuron", "ionizing radiation", "UV", "DNA damage",
                "embryo", "embryonic"]
        if s == 2 and gid not in S:
            for k in lows:
                if k.lower() in nm:
                    s = 1
                    break
    scores.append((gid, s))

with open(out, "w") as f:
    f.write("pathway\tscore\n")
    for g, s in scores:
        f.write(f"{g}\t{s}\n")
print(len(scores))
from collections import Counter
print(Counter(s for _, s in scores))
