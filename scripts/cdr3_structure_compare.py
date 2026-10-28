"""Compare the CDR3-beta loop conformations of the AS motif family, from crystal structures.

Yang et al. 2022 (Nature 612:771) deposited seven TCR-peptide-HLA-B*27:05
complexes. Three of them are family members bound to the *same* peptide
(PRPF3), which makes the comparison clean: same groove, same antigen, different
CDR3-beta sequence.

    7N2S  AS3.1  CASSVGLYSTDTQYF  + PRPF3     <- this project's clone
    7N2N  AS4.2  CASSVGLFSTDTQYF  + PRPF3     <- Hanson's headline clone
    7N2R  AS4.3  CASSVATYSTDTQYF  + PRPF3     <- our screen's candidate
    7N2O  AS4.2  + YeiH (bacterial)
    7N2Q  AS4.3  + YeiH
    7N2P  AS4.3  + RNASEH2B
    8CX4  AS8.4  + YeiH

Every structure is superposed on the HLA-B27 alpha1/alpha2 platform, not on the
TCR, so the loops are compared where it matters: their position over the
peptide-binding groove. Superposing on the loops themselves would flatter the
comparison by construction.

Outputs a JSON of traces (for plotting) and prints pairwise RMSDs.

    python scripts/cdr3_structure_compare.py --pdb-dir <dir of .pdb> --out <json>
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}

# CDR3-beta of each TCR, and the HLA-B27 heavy-chain signature used to find it.
CDR3B = {
    "AS3.1": "CASSVGLYSTDTQYF",
    "AS4.2": "CASSVGLFSTDTQYF",
    "AS4.3": "CASSVATYSTDTQYF",
}
ENTRIES = {
    "7N2S": ("AS3.1", "PRPF3"),
    "7N2N": ("AS4.2", "PRPF3"),
    "7N2R": ("AS4.3", "PRPF3"),
    "7N2O": ("AS4.2", "YeiH"),
    "7N2Q": ("AS4.3", "YeiH"),
    "7N2P": ("AS4.3", "RNASEH2B"),
    "8CX4": ("AS8.4", "YeiH"),
}
HLA_SIG = "GSHSMRYFHTSVSRPGRGEPRFITVGYVDDT"   # HLA class I alpha1 N-terminus


def read_pdb(path: str):
    """{chain: [(resseq, resname, {atom: xyz})]} for the first model only."""
    chains: dict[str, dict[int, tuple[str, dict]]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM"):
                continue
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            resname = line[17:20].strip()
            if resname not in AA3:
                continue
            ch = line[21]
            resseq = int(line[22:26])
            atom = line[12:16].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            res = chains.setdefault(ch, {}).setdefault(resseq, (resname, {}))
            res[1][atom] = xyz
    return {c: dict(sorted(r.items())) for c, r in chains.items()}


def seq_of(chain: dict) -> tuple[str, list[int]]:
    keys = sorted(chain)
    return "".join(AA3[chain[k][0]] for k in keys), keys


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation+translation mapping P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def ca(chain: dict, keys) -> np.ndarray:
    return np.array([chain[k][1]["CA"] for k in keys if "CA" in chain[k][1]])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--out", default="data/panel/cdr3_structures.json")
    args = ap.parse_args()

    parsed = {}
    for pdb, (tcr, pep) in ENTRIES.items():
        path = os.path.join(args.pdb_dir, f"{pdb}.pdb")
        if not os.path.exists(path):
            print(f"  missing {path}")
            continue
        chains = read_pdb(path)
        hla_ch = cdr3_ch = pep_ch = None
        cdr3_keys = []
        for c, res in chains.items():
            s, keys = seq_of(res)
            if HLA_SIG[:20] in s and hla_ch is None:
                hla_ch = (c, keys[:180])
            if tcr in CDR3B:
                i = s.find(CDR3B[tcr])
                if i >= 0:
                    cdr3_ch, cdr3_keys = c, keys[i:i + len(CDR3B[tcr])]
            if 8 <= len(s) <= 11 and pep_ch is None and c != (hla_ch or ("", ))[0]:
                pep_ch = (c, keys, s)
        if hla_ch is None:
            print(f"  {pdb}: HLA chain not found")
            continue
        parsed[pdb] = dict(tcr=tcr, peptide=pep, chains=chains, hla=hla_ch,
                           cdr3=(cdr3_ch, cdr3_keys), pep_chain=pep_ch)
        seq = seq_of(chains[cdr3_ch])[0][
            seq_of(chains[cdr3_ch])[1].index(cdr3_keys[0]):
        ][:len(CDR3B.get(tcr, ""))] if cdr3_ch else "-"
        print(f"{pdb}  {tcr:<6} {pep:<9} HLA chain {hla_ch[0]}  "
              f"CDR3b chain {cdr3_ch or '-'} ({len(cdr3_keys)} res)  "
              f"peptide {pep_ch[2] if pep_ch else '-'}")

    ref = "7N2S" if "7N2S" in parsed else next(iter(parsed))
    rc, rk = parsed[ref]["hla"]
    ref_ca = ca(parsed[ref]["chains"][rc], rk)

    traces = {}
    for pdb, d in parsed.items():
        c, k = d["hla"]
        P = ca(d["chains"][c], k)
        n = min(len(P), len(ref_ca))
        R, t = kabsch(P[:n], ref_ca[:n])
        rms = np.sqrt(((P[:n] @ R.T + t - ref_ca[:n]) ** 2).sum(1).mean())
        entry = dict(tcr=d["tcr"], peptide=d["peptide"], hla_rmsd=round(float(rms), 2))
        if d["cdr3"][0]:
            L = ca(d["chains"][d["cdr3"][0]], d["cdr3"][1])
            entry["cdr3"] = (L @ R.T + t).round(2).tolist()
            entry["cdr3_seq"] = CDR3B[d["tcr"]]
        if d["pep_chain"]:
            pc, pk, ps = d["pep_chain"]
            entry["pep_trace"] = (ca(d["chains"][pc], pk) @ R.T + t).round(2).tolist()
            entry["pep_seq"] = ps
        traces[pdb] = entry
        print(f"  {pdb}: HLA superposition RMSD {rms:.2f} A")

    print("\nPAIRWISE CDR3-beta RMSD (A), after superposing on the HLA platform")
    ids = [p for p in traces if "cdr3" in traces[p]]
    print("        " + "".join(f"{p:>9}" for p in ids))
    for a in ids:
        row = f"{a:<8}"
        for b in ids:
            A_, B_ = np.array(traces[a]["cdr3"]), np.array(traces[b]["cdr3"])
            n = min(len(A_), len(B_))
            row += f"{np.sqrt(((A_[:n]-B_[:n])**2).sum(1).mean()):>9.2f}"
        print(row + f"   {traces[a]['tcr']} + {traces[a]['peptide']}")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(traces, fh)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
