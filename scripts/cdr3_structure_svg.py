import json
import numpy as np

T = json.load(open(r"C:\Coding\as-experiment\data\panel\cdr3_structures.json"))
USE = ["7N2S", "7N2N", "7N2R", "7N2O"]          # 7N2Q excluded: second copy in the ASU
pep = np.array(T["7N2S"]["pep_trace"])
loops = {k: np.array(T[k]["cdr3"]) for k in USE}

x = pep[-1] - pep[0]
x /= np.linalg.norm(x)
up = loops["7N2S"].mean(0) - pep.mean(0)
up -= x * (up @ x)
z = up / np.linalg.norm(up)
y = np.cross(z, x)
origin = pep.mean(0)
F = np.vstack([x, y, z])


def proj(P, a, b):
    Q = (P - origin) @ F.T
    return Q[:, [a, b]]


def path(pts, sx, sy, cx, cy):
    d = f"M {cx + pts[0][0]*sx:.1f} {cy - pts[0][1]*sy:.1f}"
    for p in pts[1:]:
        d += f" L {cx + p[0]*sx:.1f} {cy - p[1]*sy:.1f}"
    return d


for name, (a, b) in (("SIDE (along groove x, up z)", (0, 2)), ("TOP (x, y)", (0, 1))):
    print("===", name)
    allpts = np.vstack([proj(pep, a, b)] + [proj(v, a, b) for v in loops.values()])
    lo, hi = allpts.min(0), allpts.max(0)
    print(f"   range x {lo[0]:.1f}..{hi[0]:.1f}   y {lo[1]:.1f}..{hi[1]:.1f}")
    S = 7.0
    cx, cy = -lo[0] * S + 40, hi[1] * S + 30
    print(f"   viewBox 0 0 {int((hi[0]-lo[0])*S+80)} {int((hi[1]-lo[1])*S+60)}")
    print(f'   PEP  <path d="{path(proj(pep, a, b), S, S, cx, cy)}"/>')
    for k in USE:
        print(f'   {k} ({T[k]["tcr"]}/{T[k]["peptide"]})  '
              f'<path d="{path(proj(loops[k], a, b), S, S, cx, cy)}"/>')
    # residue-98 position (index 7 of CASSVGL[Y/F]STDTQYF)
    for k in USE:
        p = proj(loops[k][7:8], a, b)[0]
        print(f'   {k} res98 cx={cx + p[0]*S:.1f} cy={cy - p[1]*S:.1f}')
    for i, p in enumerate(proj(pep, a, b)):
        print(f'   P{i+1} cx={cx + p[0]*S:.1f} cy={cy - p[1]*S:.1f}')
