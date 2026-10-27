"""Test the paper's functionally-near clonotypes against our own AS cohorts.

Where the candidates come from. `prp_divergence.py` computes, from Wang/Yeh/
Garcia's cached peptide recognition profiles, the functional distance (motif
Jensen-Shannon, their definition) between our clonotype `CASSVGLYSTDTQYF` and
the other 15 betas in their HLA-B*27:05 panel. Edit distance predicts that
distance not at all (Spearman +0.113, p = 0.10, flat from 0 to 7 edits), and the
target's nearest functional relatives sit **three** edits away -- invisible to
the radius-1 search that killed both of our earlier widening attempts.

Three of those relatives are absent from our own 167-sequence burden panel:

    CASSLEVYSTDTQYF  JS 0.1750, 3 edits, 10/786 Emerson healthy
    CASSLGLYSTMEQYF  JS 0.1873, 3 edits, <7/786 (absent from the public set)
    CASSVATYSTDTQYF  JS 0.2375, 2 edits,  7/786, and present in GSE288581 AS

This scans all 16 panel sequences -- not just the three -- against the Komech
deposits, because one pass costs the same and 16 sequences with known functional
distances turn a lookup into a testable question: does functional distance to
the target predict AS carriage?

Two things this gets right that the obvious version does not.

**A uniform detection floor.** The cached CDR3 tables kept every sequence seen
twice or more, *plus* singletons within 2 edits of `CASSVGLYSTDTQYF`. So the
target and its 2-edit neighbours are detectable at one read while everything
else needs two. Comparing them as deposited would hand the target a sensitivity
advantage over the very candidates being tested. `MIN_READS = 2` is therefore
applied to all 16, the target included, and a subject whose floor sits above the
carriage threshold is **unevaluable rather than negative** -- the rule that
already set aside 7 of 32 AS subjects in the blood result.

**Numerator and denominator in one code path.** Both come from the same run's
table: reads for the sequence, `junction_reads` for the depth. No motif counted
on two strands against a denominator counted on one, which is the bug that once
inflated every frequency about 34-fold.

Carriage is aggregated per *subject*, never per run, and only unsorted whole
PBMC or blood is used -- a CD4-sorted library cannot contain a CD8 clone, so
counting one as a negative manufactures a false negative.

    modal run scripts/prp_candidates_scan.py::scan
"""
import modal

app = modal.App("as-prp-scan")
vol = modal.Volume.from_name("as-data", create_if_missing=True)

TARGET = "CASSVGLYSTDTQYF"

# The paper's 16 distinct panel betas with their motif JS to the target,
# computed by prp_divergence.py from the cached profiles. Written here rather
# than recomputed so this scan is reproducible without the 32 MB of profiles.
PANEL_JS = {
    "CASSVGLYSTDTQYF": 0.1065,
    "CASSLEVYSTDTQYF": 0.1750,
    "CASSLGLYSTMEQYF": 0.1873,
    "CASSPGTFSTDTQYF": 0.2052,
    "CASSPATYSTDTQYF": 0.2064,
    "CASSLATYSTDTQYF": 0.2207,
    "CASSVGLFSTDTQYF": 0.2238,
    "CASSVAYYSTDTQYF": 0.2321,
    "CASSVGLYSTGELFF": 0.2341,
    "CASSDGYSTDTQYF": 0.2343,
    "CASSPWGTYSTDTQYF": 0.2368,
    "CASSVATYSTDTQYF": 0.2375,
    "CASSLAYYSTDTQYF": 0.2479,
    "CASSVGTYSTDTQYF": 0.2485,
    "CASSGGLYSTDTQYF": 0.2600,
    "CASSPAYYSTDTQYF": 0.2622,
}

MIN_READS = 2          # uniform floor; see the module docstring
PROJECTS = ("PRJNA393498", "PRJEB51315")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("pandas==2.2.3", "pyarrow==17.0.0")
)


@app.function(image=image, volumes={"/data": vol}, cpu=4.0, memory=8192,
              timeout=3600)
def scan():
    """Pull the 16 panel sequences out of every cached CDR3 table."""
    import glob
    import gzip
    import json
    import os

    rows = []
    for proj in PROJECTS:
        metas = sorted(glob.glob(f"/data/tcr/{proj}/*_cdr3.json"))
        print(f"{proj}: {len(metas)} run tables", flush=True)
        for mp in metas:
            with open(mp) as fh:
                meta = json.load(fh)
            if meta.get("status") != "ok":
                continue
            tp = mp.replace("_cdr3.json", "_cdr3.tsv.gz")
            if not os.path.exists(tp):
                continue
            found = {}
            with gzip.open(tp, "rt") as fh:
                next(fh, None)
                for line in fh:
                    c, _, n = line.partition("\t")
                    if c in PANEL_JS:
                        found[c] = int(n)
            rows.append({
                "project": proj,
                "run": meta["run"],
                "disease": meta.get("disease", ""),
                # 2018 uses `subject`/`cell_type`; 2022 uses `patient`/
                # `tissue`/`subset`. Normalise so one aggregation serves both.
                "subject": meta.get("subject") or meta.get("patient", ""),
                "cell_type": meta.get("cell_type", ""),
                "tissue": meta.get("tissue", ""),
                "subset": meta.get("subset", ""),
                "junction_reads": int(meta.get("junction_reads", 0)),
                "n_reads": int(meta.get("n_reads", 0)),
                **{f"q_{c}": found.get(c, 0) for c in PANEL_JS},
            })
    out = "/data/tcr/prp_panel_scan.json"
    with open(out, "w") as fh:
        json.dump(rows, fh)
    vol.commit()
    print(f"wrote {out}: {len(rows)} runs")
    return rows
