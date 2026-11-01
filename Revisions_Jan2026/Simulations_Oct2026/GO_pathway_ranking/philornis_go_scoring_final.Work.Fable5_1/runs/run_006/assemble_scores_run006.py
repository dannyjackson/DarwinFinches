import sys
d="/sessions/sleepy-happy-fermat/mnt/outputs/runs/run_006/"
scores={}
for i in range(1,6):
    for line in open(d+f"scores_part{i}_run006.txt"):
        line=line.strip()
        if not line: continue
        g,s=line.split()
        assert g not in scores, g
        scores[g]=int(s)
ids=[l.split()[0] for l in open(d+"input_terms.txt").read().splitlines()[1:] if l.strip()]
missing=[g for g in ids if g not in scores]
extra=[g for g in scores if g not in set(ids)]
print("input",len(ids),"scored",len(scores),"missing",missing,"extra",extra)
assert not missing and not extra and len(ids)==len(set(ids))
with open(d+"scores.tsv","w") as f:
    f.write("pathway\tscore\n")
    for g in ids: f.write(f"{g}\t{scores[g]}\n")
from collections import Counter
print(Counter(scores[g] for g in ids))
print("rows",len(ids))
