#!/usr/bin/env python3
"""Map GO IDs in input_terms.txt to a priori relevance scores (1-5) for a host
response to Philornis downsi. Every score below was assigned by hand after a
literature review (see review.md). Any GO ID not listed explicitly receives 1
(no specific link)."""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.path.join(HERE, "input_terms.txt")
OUT = os.path.join(HERE, "scores.tsv")

S = {}
def s(score, *ids):
    for i in ids:
        S[i] = score

# ---- 5: directly captures a documented host mechanism ----
s(5, "GO:0002250",  # adaptive immune response (P. downsi-binding antibodies documented)
     "GO:0030218", "GO:0045648",  # erythrocyte differentiation / pos reg (blood loss, anaemia, Hb reduction)
     "GO:0009611", "GO:0042060", "GO:0090303",  # response to wounding, wound healing (larval lesions, naris damage)
     "GO:0006954")  # inflammatory response at larval feeding sites

# ---- 4: core component or strong a priori prediction ----
s(4, "GO:0050728", "GO:0050727", "GO:0050729",  # regulation of inflammatory response
     "GO:0050776", "GO:0002682", "GO:0002376", "GO:0006955", "GO:0006952",  # immune response / immune system / defense response
     "GO:0045087", "GO:0045089", "GO:0002218",  # innate immune response and its activation
     "GO:0002224",  # toll-like receptor signaling (innate detection, secondary bacterial infection of wounds)
     "GO:0001817", "GO:0001819", "GO:0019221",  # cytokine production / cytokine-mediated signaling
     "GO:0050853", "GO:0042113", "GO:0030183", "GO:0030890",  # B cell receptor signaling, B cell activation/differentiation/proliferation (antibody response)
     "GO:0006956",  # complement activation
     "GO:0001666", "GO:0071456",  # response to hypoxia (anaemia -> tissue hypoxia, EPO)
     "GO:0002244", "GO:0030097",  # hematopoietic progenitor differentiation, hemopoiesis (replace lost blood)
     "GO:0006879",  # iron ion homeostasis (haemoglobin resynthesis after blood loss)
     "GO:0007596",  # blood coagulation (haemostasis at feeding wounds)
     "GO:0097009",  # energy homeostasis (tolerance via energy compensation)
     "GO:0008544", "GO:0043588", "GO:0030216")  # epidermis / skin development, keratinocyte differentiation (skin barrier, repair)

# ---- 3: plausible downstream or supporting process ----
s(3, "GO:0033209", "GO:0032760", "GO:0032720", "GO:0034612", "GO:0071356",  # TNF signaling/production/response
     "GO:0097696", "GO:0007259", "GO:0046427",  # STAT / JAK-STAT
     "GO:0033554", "GO:0006950",  # response to stress
     "GO:0042554", "GO:2000379", "GO:2000377", "GO:0034614", "GO:0070301", "GO:0034599", "GO:0006979", "GO:0098869",  # ROS / oxidative stress
     "GO:1900745",  # p38 MAPK
     "GO:0030511", "GO:0071560", "GO:0007179",  # TGF-beta (wound healing, immune regulation)
     "GO:0042594", "GO:0009267", "GO:0031669", "GO:0031667",  # starvation / nutrient levels (energy compensation)
     "GO:0001818", "GO:0032715", "GO:0032755", "GO:0032757", "GO:0032735", "GO:0032729", "GO:0032722",  # cytokine/chemokine production regulation
     "GO:0071347", "GO:0071346", "GO:0034097", "GO:0071345",  # response to IL-1, IFN-gamma, cytokines
     "GO:1900227", "GO:0002753",  # inflammasome, cytoplasmic PRR signaling
     "GO:0097400",  # IL-17 signaling (skin/mucosal antibacterial)
     "GO:0050852", "GO:0030217", "GO:0042110", "GO:0050870", "GO:0050868",  # T cell signaling/differentiation/activation
     "GO:0045824", "GO:0038061", "GO:0043122", "GO:0043123", "GO:0043124", "GO:0007249",  # innate regulation, NF-kB
     "GO:0034142", "GO:0071222", "GO:0032496",  # TLR4, LPS
     "GO:0050829", "GO:0042742",  # defense response to bacteria (secondary wound infection)
     "GO:0006935", "GO:0060326", "GO:0007159",  # chemotaxis, leukocyte adhesion
     "GO:0006909", "GO:0050766", "GO:0006911",  # phagocytosis
     "GO:0030855", "GO:0050679", "GO:0010634", "GO:0048146",  # epithelial differentiation/proliferation/migration, fibroblast proliferation (re-epithelialization)
     "GO:0030199", "GO:0071711", "GO:0030198", "GO:0007229",  # collagen, basement membrane, ECM, integrin (dermal repair)
     "GO:0045766", "GO:0010575", "GO:0038084", "GO:0001525",  # angiogenesis / VEGF (granulation tissue)
     "GO:0006006", "GO:0042593", "GO:0006094", "GO:0006635", "GO:0019915",  # glucose metabolism/homeostasis, gluconeogenesis, beta-oxidation, lipid storage (energy budget)
     "GO:0060021",  # roof of mouth development (nasal cavity / palate damaged by first instars)
     "GO:0070527", "GO:0030168",  # platelet aggregation/activation (haemostasis)
     "GO:0045109",  # intermediate filament organization (keratin; larvae consume beak keratin)
     "GO:0019233",  # sensory perception of pain (detection of feeding larvae -> preening/removal)
     "GO:0048009",  # IGF receptor signaling (nestling growth, body-size-mediated tolerance)
     "GO:0045618")  # pos reg keratinocyte differentiation

# ---- 2: weak, indirect or speculative link ----
s(2, "GO:0090051", "GO:0043542", "GO:0010595", "GO:0090050", "GO:0061028", "GO:0002040", "GO:0016525",
     "GO:0045765", "GO:0002042", "GO:0001938", "GO:0043536", "GO:0001568", "GO:0035924", "GO:0048010",
     "GO:0046627", "GO:0008286", "GO:0032868", "GO:0032869", "GO:0032024", "GO:0035774", "GO:0050796",
     "GO:0046676", "GO:0030073", "GO:0071333", "GO:0042149", "GO:0010906", "GO:0001678", "GO:0046326",
     "GO:0050650", "GO:0030200", "GO:0015012", "GO:0030166",
     "GO:0030513", "GO:0030509", "GO:0030514", "GO:0071773", "GO:0040037", "GO:0008543", "GO:0032924",
     "GO:0006865", "GO:0015804", "GO:0003333", "GO:0034198", "GO:0071230", "GO:0006520", "GO:0006541",
     "GO:0030148", "GO:0006672", "GO:0046513",
     "GO:0007043", "GO:0070830", "GO:0034332", "GO:0098609", "GO:0016339", "GO:0031589", "GO:0034113",
     "GO:0007155", "GO:0001954", "GO:0007160", "GO:0045785", "GO:0033627", "GO:0044331",
     "GO:0006509", "GO:0030522", "GO:0030518", "GO:0009725", "GO:0009755", "GO:0032870",
     "GO:0006829", "GO:0006882", "GO:0071577", "GO:0030001",
     "GO:0070534", "GO:0017121",
     "GO:0010506", "GO:0016239", "GO:0006914", "GO:0016241", "GO:0016236", "GO:0010508",
     "GO:0071364", "GO:0045742", "GO:0007173", "GO:0071363", "GO:0007169", "GO:0007178",
     "GO:0042632", "GO:0034063",
     "GO:0005978", "GO:0001837", "GO:0010718",
     "GO:0035589", "GO:0033198",
     "GO:0030501", "GO:0030282", "GO:0030500", "GO:0060348", "GO:0001501", "GO:0001503", "GO:0045669", "GO:0001649",
     "GO:0051216", "GO:0002062",
     "GO:2001237", "GO:2001238", "GO:1902042", "GO:0043065", "GO:0042981", "GO:0097190", "GO:0043066", "GO:0006915", "GO:0097191",
     "GO:0001558", "GO:0030308", "GO:0030307", "GO:0008284",
     "GO:0043408", "GO:0043410", "GO:0043409", "GO:0000165", "GO:0070374", "GO:0070372", "GO:0070373",
     "GO:0046329", "GO:0007254", "GO:0046330",
     "GO:0008610", "GO:0006633", "GO:0016042", "GO:0019432", "GO:0006631", "GO:0001676", "GO:0015908",
     "GO:0055088", "GO:0006869", "GO:0070328", "GO:0019216", "GO:0006629",
     "GO:0048839", "GO:0007605",
     "GO:0120163", "GO:0120162",
     "GO:0009888", "GO:0002009",
     "GO:0060391", "GO:0060395", "GO:0030512",
     "GO:0001755",
     "GO:0006099", "GO:0046034", "GO:0006096", "GO:0005975", "GO:0006120", "GO:0009060",
     "GO:0048870", "GO:0030334", "GO:0016477", "GO:0030335",
     "GO:0009408", "GO:0034605",
     "GO:0032008", "GO:0031929", "GO:0032007", "GO:1904262", "GO:1904263", "GO:0038202",
     "GO:0045444", "GO:0045600", "GO:0045599",
     "GO:0032922", "GO:0007623", "GO:0042752",
     "GO:0110076",
     "GO:1901223", "GO:1901224",
     "GO:0006887", "GO:0045055", "GO:0009306", "GO:0050714",
     "GO:0008217", "GO:0002027",
     "GO:0043491", "GO:0051897",
     "GO:0006493",
     "GO:0007218",
     "GO:0001778",
     "GO:0042572",
     "GO:0009410", "GO:0071466",
     "GO:0050860", "GO:0046426",
     "GO:0032728",
     "GO:0034976", "GO:0030968",
     "GO:0045454",
     "GO:0071260",
     "GO:0050918",
     "GO:0050680")

def main():
    rows = []
    with open(INP) as fh:
        lines = fh.read().splitlines()
    for line in lines[1:]:
        if not line.strip():
            continue
        go = line.split(" ", 1)[0]
        rows.append((go, S.get(go, 1)))
    with open(OUT, "w") as fh:
        fh.write("pathway\tscore\n")
        for go, sc in rows:
            fh.write(f"{go}\t{sc}\n")
    # sanity checks
    ids = [r[0] for r in rows]
    assert len(ids) == len(set(ids)), "duplicate IDs in input"
    missing = [k for k in S if k not in set(ids)]
    if missing:
        print("WARNING: scored IDs not in input:", missing, file=sys.stderr)
    from collections import Counter
    print(Counter(sc for _, sc in rows), file=sys.stderr)
    print(len(rows))

if __name__ == "__main__":
    main()
