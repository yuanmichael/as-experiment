"""Rebuild the 167-sequence AS burden panel, the artefact every downstream score needs.

The panel was originally computed inline and never persisted, which made the
result unreproducible outside the session that produced it. This script is the
definition of record. It writes `data/panel/as_burden_panel_167.tsv`.

The selection rule, in full:

    keep a CDR3 if it appears in >= MIN_AS_CARRIERS of the 25 Komech AS patients
    (presence = >=2 reads after thinning every subject to THIN_DEPTH junction
    reads) AND in <= MAX_EMERSON_PCT of Emerson 2017's 786 healthy repertoires.

MIN_AS_CARRIERS is 5, not 4. The original inline code applied `>= 4`, but to a
candidate pool already restricted to clonotypes in >= 5 AS patients, so the 4
never bound and the effective rule was 5. Relaxing it to a true 4 yields 1,014
sequences rather than 167 -- a different panel, which has not been evaluated.

Both arms matter and neither substitutes for the other. The AS arm alone selects
public clonotypes that most of humanity carries -- of 21,206 sequences in >=5 AS
patients, 90.6% sit in more than a tenth of Emerson. The Emerson arm alone
selects deposit-specific artefacts: the CRVSGL[WGCR]AKNIQYF family is absent
from all 786 external repertoires yet carried by 30-50% of Komech's own healthy
controls. Four of those artefacts survive into the panel by design -- they are
disease-indifferent, so they dilute the score rather than inflate it, and
removing them post hoc would be fitting the panel to its own evaluation.

Thinning is not optional. Presence is depth-sensitive, Komech's libraries span
17-fold in depth, and an unthinned panel measures sequencing depth.

Scores from this panel: AS 38 of 167 sequences (median, range 10-56), healthy
6-9 (range 0-16). At the 0.10 threshold (17 sequences) sensitivity is 84% and
specificity 100% against 51 healthy controls from two cohorts.

    python scripts/build_panel.py --tables <dir of Komech *_cdr3.tsv.gz> \
        --emerson data/airr/emerson-2017-natgen
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clono_sweep import carrier_matrix, load_tables, thin  # noqa: E402

THIN_DEPTH = 800_000        # near the shallowest Komech library, not far below it
MIN_AS_CARRIERS = 5
MAX_EMERSON_PCT = 0.05
SEED = 0
TARGET = "CASSVGLYSTDTQYF"  # the known clone, expected in the panel as a sanity check


def emerson_prevalence(root: str, candidates: set[str]) -> tuple[dict[str, int], int]:
    """How many of Emerson's 786 healthy repertoires carry each candidate."""
    meta = pd.read_csv(f"{root}/metadata.tsv", sep="\t")
    counts = dict.fromkeys(candidates, 0)
    n = 0
    for sid in meta["subject_id"]:
        f = f"{root}/reps/{sid}.parquet"
        if not os.path.exists(f):
            continue
        t = pd.read_parquet(f, columns=["junction_aa"])
        n += 1
        for c in candidates.intersection(t["junction_aa"].values):
            counts[c] += 1
    return counts, n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True,
                    help="directory of Komech *_cdr3.tsv.gz plus *_cdr3.json")
    ap.add_argument("--emerson", default="data/airr/emerson-2017-natgen")
    ap.add_argument("--out", default="data/panel/as_burden_panel_167.tsv")
    ap.add_argument("--case", default="as")
    ap.add_argument("--min-as-carriers", type=int, default=MIN_AS_CARRIERS)
    ap.add_argument("--max-emerson-pct", type=float, default=MAX_EMERSON_PCT)
    args = ap.parse_args()

    rng = np.random.default_rng(SEED)
    df = load_tables(args.tables)
    thinned = thin(df, THIN_DEPTH, rng)
    mat, subjects, cdr3s = carrier_matrix(thinned, "presence")
    grp = (thinned.drop_duplicates("subject").set_index("subject")["group"]
           .reindex(subjects).to_numpy())
    is_case = grp == args.case
    print(f"{int(is_case.sum())} cases, {int((~is_case).sum())} controls, "
          f"{mat.shape[1]:,} clonotypes after thinning to {THIN_DEPTH:,}")

    as_carriers = np.asarray(is_case.astype(np.int32) @ mat).ravel()
    hc_carriers = np.asarray((~is_case).astype(np.int32) @ mat).ravel()

    # Screen only what could possibly qualify, so Emerson is read once.
    shortlist = as_carriers >= args.min_as_carriers
    cands = set(cdr3s[shortlist].tolist())
    print(f"{len(cands):,} clonotypes in >= {args.min_as_carriers} cases; "
          f"screening against Emerson")
    prev, n_emerson = emerson_prevalence(args.emerson, cands)
    print(f"screened {n_emerson} Emerson healthy repertoires")

    ev = np.array([prev.get(c, 0) for c in cdr3s])
    keep = shortlist & (ev <= args.max_emerson_pct * n_emerson)
    panel = pd.DataFrame({
        "cdr3": cdr3s[keep],
        "as_carriers": as_carriers[keep],
        "komech_healthy_carriers": hc_carriers[keep],
        "emerson": ev[keep],
        "emerson_pct": ev[keep] / n_emerson,
    }).sort_values(["as_carriers", "emerson"], ascending=[False, True])

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    panel.to_csv(args.out, sep="\t", index=False, float_format="%.5f")
    print(f"\npanel: {len(panel)} sequences -> {args.out}")
    if TARGET in set(panel["cdr3"]):
        print(f"  sanity check: {TARGET} is in the panel")
    else:
        print(f"  WARNING: {TARGET} is NOT in the panel -- selection has drifted")


if __name__ == "__main__":
    main()
