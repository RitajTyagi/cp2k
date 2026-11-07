#!/usr/bin/env python3
"""Compare the Gamma-folded G0W0 levels of a Gamma-only n x n x n diamond supercell with the
full-k PR 6447 levels at the primitive k-points that fold onto Gamma (path points).
usage: compare_folded.py <our bandstructure file> <PR bandstructure file>"""
import sys, re
def read(f):
    kp, cur = {}, None
    for l in open(f, encoding='utf-8'):
        m = re.match(r'\s+kpoint:\s+(\d+)\s+coordinate:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', l)
        if m:
            cur = tuple(round(float(x), 4) for x in m.groups()[1:]); kp[cur] = []; continue
        m = re.match(r'\s+(\d+)\s+\((occ|vir)\)\s+\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', l)
        if m and cur is not None:
            kp[cur].append((m.group(2), float(m.group(3)), float(m.group(4)), float(m.group(5)), float(m.group(7))))
    return kp
ours = read(sys.argv[1]); pr = read(sys.argv[2])
o = ours[(0.0, 0.0, 0.0)]
print("ours (Gamma of the supercell): %d levels" % len(o))
ov = [x for x in o if x[0] == 'occ']; oc = [x for x in o if x[0] == 'vir']
print("  DFT VBM/CBM  %8.3f %8.3f  gap %6.3f" % (max(v[1] for v in ov), min(v[1] for v in oc), min(v[1] for v in oc) - max(v[1] for v in ov)))
print("  G0W0 VBM/CBM %8.3f %8.3f  gap %6.3f" % (max(v[4] for v in ov), min(v[4] for v in oc), min(v[4] for v in oc) - max(v[4] for v in ov)))
print("PR path points:")
for k, lev in pr.items():
    v = [x for x in lev if x[0] == 'occ']; c = [x for x in lev if x[0] == 'vir']
    print("  k=%s  DFT VBM/CBM %8.3f %8.3f  G0W0 VBM/CBM %8.3f %8.3f" % (k, max(x[1] for x in v), min(x[1] for x in c), max(x[4] for x in v), min(x[4] for x in c)))
# state-by-state: match our folded DFT levels to PR DFT levels at each point by energy
print("\nState matching by DFT energy (|dE_DFT| < 0.02 eV): DFT, Sigma_c ours/PR, Sigma_x ours/PR, G0W0 ours/PR")
used = set()
for k, lev in pr.items():
    for x in lev:
        best = None
        for i, y in enumerate(o):
            if i in used or y[0] != x[0]: continue
            if abs(y[1] - x[1]) < 0.02 and (best is None or abs(y[1]-x[1]) < abs(o[best][1]-x[1])): best = i
        if best is None: continue
        used.add(best); y = o[best]
        print("  k=%-20s %s %8.3f | %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f  d=%+.3f" % (k, x[0], x[1], y[2], x[2], y[3], x[3], y[4], x[4], y[4]-x[4]))
