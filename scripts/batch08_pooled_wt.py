#!/usr/bin/env python3
"""Pool bb54-g2's WT arm across batch-07 and batch-08, criterion verbatim.

n=3 is not enough for a WT margin - batch-06's -59.2 %/-73.1 % was a three-draw
excursion that pooled (9 G12D vs 12 WT) to -33.4 %/-52.4 %.  Batch-08 re-ran the
same sequence in the same cell (`nopot` == r3np) as the paired control for the
de-wetted twin, which adds 3 G12D and 3 WT draws.  This pools every r3np draw of
that one sequence across both batches and applies
`batch06_panel.wt_verdict` - the batch-06 criterion, unchanged and not retuned.
No new API calls; this only re-reads checkpoints.
"""
from __future__ import annotations

import json, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_panel as P  # noqa: E402

PARENT = "MMLPDDAIAELSGTPEARELARRLGGDERVAINYRAGRRGRVGVEA"
B07_G = ["r3np", "r3np_rep3", "r3np_rep4"]
B07_W = ["r3np", "r3np_rep2", "r3np_rep3", "r3np_rep4"]


def rows(path: Path):
    return [json.loads(l) for l in path.read_text().splitlines() if l.strip()]


def main() -> int:
    b07 = rows(REPO / "results/batch-07/panel_checkpoint.jsonl")
    b08 = rows(REPO / "results/batch-08/panel_checkpoint.jsonl")
    def sel(src, arm, cells):
        return [r["row"] for r in src if r["sequence"] == PARENT
                and r["arm"] == arm and r["condition"] in cells]
    g = sel(b07, "g12d", B07_G) + sel(b08, "g12d", ["nopot"])
    w = sel(b07, "wt", B07_W) + sel(b08, "wt", ["nopot"])
    gates = sum(1 for r in g if str(r.get("asp6_gate_pass")).lower() == "true")
    out = {"design": "bb54-g2 (B6-c6-46_w2_bb54 donor-anchored-g2)",
           "cell": "r3np (batch-08 label `nopot`)",
           "pooled_across": "batch-07 r3np/_rep2/_rep3/_rep4 + batch-08 nopot",
           "n_g12d": len(g), "n_wt": len(w),
           "g12d_asp6_gate_passes": gates,
           "batch08_only_n_g12d": 3, "batch08_only_n_wt": 3,
           "wt_criterion": P.WT_CRITERION,
           "pooled_wt_verdict": P.wt_verdict(g, w)}
    (REPO / "results/batch-08/pooled_wt.json").write_text(
        json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str)[:1200])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
