import sys
s = {}
def g(score, *ids):
    for i in ids:
        s["GO:%07d" % i] = score

# 5: directly documented (IgY antibody response; blood loss/hemoglobin recovery; open wounds)
g(5, 2250, 6955, 42113, 30183, 30218, 45648, 42060, 9611, 90303)
# 4: core components / strong a priori predictions
g(4, 45087, 6952, 2376, 2682, 50776, 2218, 45089, 45824 if False else 0)
g(4, 6954, 50727, 50729, 50728, 6956, 50853, 30890,
  42110, 30217, 50852, 50870,
  30097, 2244, 6879, 1666, 71456,
  7596, 30168, 70527)
g(3, 45824, 50868, 50860)
# 3: plausible downstream / supporting
g(3, 32496, 71222, 2224, 34142, 33209, 71356, 34612,
  1817, 1819, 1818, 32760, 32720, 32755, 32715, 32757, 32735, 32722, 32729,
  7249, 43122, 43123, 43124, 38061, 1901223, 1901224,
  19221, 71345, 34097, 7259, 46427, 97696, 71346, 71347, 97400, 1900227,
  6909, 6911, 50766, 42554, 7159, 60326,
  30199, 30198, 48146, 30216, 45618, 8544, 43588,
  97009, 42593, 31667, 48009, 8217, 19233)
# 2: weak / indirect / speculative
g(2, 55117, 71260, 50680, 30512, 38202, 31929, 45785,
  51607, 98586, 9615, 2230, 140374, 32728, 32481, 32480, 32727, 60339,
  2753, 42742, 50829, 46426,
  7613, 7612, 50890, 6882, 16226, 19216, 19915, 120162, 120163,
  50918, 6935, 30335, 1900745, 7254, 46330,
  60395, 60391, 7179, 71560, 30511,
  10595, 1938, 2042, 90050, 1525, 2040, 45765, 45766, 16525, 1568,
  10575, 38084, 48010, 35924, 61028,
  71364, 7173, 8543, 50679, 10634, 71711, 1778, 45109,
  2000379, 2000377, 6979, 34599, 98869, 34614, 70301,
  6950, 33198, 7623, 42752,
  1503, 30282, 30500, 30501, 30279, 60348, 1649, 45669, 45668, 30316, 45672, 1501, 30509,
  8286, 32868, 32869, 6006, 1678, 6094, 9725, 9755, 32870, 30518,
  42594, 31669, 7218, 2027)

out = []
with open("/home/claude/runs/run_016/input_terms.txt") as f:
    next(f)
    for line in f:
        line = line.strip()
        if not line:
            continue
        gid = line.split(" ", 1)[0]
        out.append((gid, s.get(gid, 1)))
with open("/home/claude/runs/run_016/scores.tsv", "w") as f:
    f.write("pathway\tscore\n")
    for gid, sc in out:
        f.write("%s\t%d\n" % (gid, sc))
print(len(out), len(set(o[0] for o in out)))
from collections import Counter
print(Counter(o[1] for o in out))
missing = [k for k in s if k not in set(o[0] for o in out) and k != "GO:0000000"]
print("unmatched ids:", missing)
