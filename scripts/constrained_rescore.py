"""Constrained four-state rescore of the batch-03 candidates.

Plan-B method change: instead of asking Boltz-2 to rediscover the designed
pose from sequence, steer it with the NIM's `constraints` pocket term
(binder chain D, contacts on peptide C4-C8, max_distance 6 A). The same
constraint is applied to every peptide-bearing state (g12d, wt, decoy) so the
specificity gap is computed under identical bias; apo has no peptide, so it is
predicted unconstrained.
"""
import sys, json, time, traceback
sys.path.insert(0,"src")
from pathlib import Path
import run_pipeline, complex_analysis, design_gates, pose_filter
from nvidia_client import NvidiaClient, BOLTZ2_PATH

OUT = Path("results/constrained-rescore"); OUT.mkdir(parents=True, exist_ok=True)
CKPT = OUT/"checkpoint.jsonl"
done = set()
if CKPT.exists():
    done = {json.loads(l)["backbone"] for l in CKPT.open()}

rows=[json.loads(l)["row"] for l in open("results/batch-03/checkpoint.jsonl")]
states = run_pipeline.load_target_states(Path("data/targets"))
contacts=[{"id":run_pipeline.PEPTIDE_CHAIN,"residue_index":i} for i in range(4,9)]
pocket=[{"binder":run_pipeline.BINDER_CHAIN,"contacts":contacts,"max_distance":6.0}]
c=NvidiaClient()

def predict(polymers, constrained):
    body={"polymers":polymers,"recycling_steps":3,"sampling_steps":50,
          "diffusion_samples":1,"step_scale":1.638,"without_potentials":False}
    if constrained: body["constraints"]=pocket
    return c.post(BOLTZ2_PATH, body)

for n,r in enumerate(rows,1):
    bb=r["backbone"]
    if bb in done: continue
    seq=r["sequence"]
    rec={"backbone":bb,"strategy":r["strategy"],"sequence":seq,"calls":0}
    try:
        polymers=[dict(p) for p in states["g12d"]]+[{"id":run_pipeline.BINDER_CHAIN,
                   "molecule_type":"protein","sequence":seq}]
        p=predict(polymers, True); rec["calls"]+=1
        st=p["structures"][0]["structure"]
        m=complex_analysis.interface_metrics(st, binder_sequence=seq,
             binder_chain_hint=run_pipeline.BINDER_CHAIN)
        rec["iptm_g12d"]=run_pipeline.binder_interface_iptm(p,len(polymers)-1)
        for k in ("peptide_buried_frac","asp6_buried_frac","buried_area_peptide",
                  "buried_area_total","pos6_contacts","net_interface_charge",
                  "asp6_hbond","asp6_salt_bridge"):
            rec[k]=m.get(k)
        (OUT/f"{bb}_g12d.cif").write_text(st)
        on_peptide = (m.get("buried_area_peptide") or 0)>=25 and (m.get("pos6_contacts") or 0)>=1
        rec["pose_ok"]=on_peptide
        if on_peptide:
            for state,constrained in (("wt",True),("decoy",True),("apo",False)):
                pol=[dict(p2) for p2 in states[state]]+[{"id":run_pipeline.BINDER_CHAIN,
                      "molecule_type":"protein","sequence":seq}]
                time.sleep(8)
                pay=predict(pol, constrained); rec["calls"]+=1
                rec[f"iptm_{state}"]=run_pipeline.binder_interface_iptm(pay,len(pol)-1)
            others=[rec.get("iptm_wt"),rec.get("iptm_apo"),rec.get("iptm_decoy")]
            others=[o for o in others if o is not None]
            if others and rec.get("iptm_g12d") is not None:
                rec["specificity_gap"]=round(rec["iptm_g12d"]-max(others),4)
    except Exception as e:
        rec["error"]=f"{type(e).__name__}: {e}"[:400]
    with CKPT.open("a") as fh: fh.write(json.dumps(rec,default=str)+"\n")
    print(f"[{n}/{len(rows)}] {bb} pose_ok={rec.get('pose_ok')} iptm={rec.get('iptm_g12d')} gap={rec.get('specificity_gap')} pep={rec.get('peptide_buried_frac')} pos6={rec.get('pos6_contacts')} hb={rec.get('asp6_hbond')} err={rec.get('error','')}", flush=True)
    time.sleep(8)
print("DONE")
