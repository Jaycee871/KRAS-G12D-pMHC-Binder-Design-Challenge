import json,sys,time,statistics
sys.path.insert(0,'src')
import run_pipeline
from nvidia_client import NvidiaClient, BOLTZ2_PATH
first={json.loads(l)["backbone"]:json.loads(l) for l in open('results/constrained-rescore/checkpoint.jsonl')}
states=run_pipeline.load_target_states(__import__('pathlib').Path('data/targets'))
contacts=[{"id":"C","residue_index":i} for i in range(4,9)]
pocket=[{"binder":"D","contacts":contacts,"max_distance":6.0}]
c=NvidiaClient(); out={}
for name in ("S1-pos6-68_bb18","S1-pos6-68_bb13"):
    seq=first[name]["sequence"]; vals=[]
    for rep in range(3):
        pol=[dict(p) for p in states["g12d"]]+[{"id":"D","molecule_type":"protein","sequence":seq}]
        body={"polymers":pol,"recycling_steps":3,"sampling_steps":50,"diffusion_samples":1,
              "step_scale":1.638,"without_potentials":False,"constraints":pocket}
        p=c.post(BOLTZ2_PATH,body)
        vals.append(round(run_pipeline.binder_interface_iptm(p,len(pol)-1),4))
        time.sleep(8)
    out[name]={"g12d_constrained_repeats":vals,"sd":round(statistics.pstdev(vals),4)}
    print(name,out[name],flush=True)
json.dump(out,open('results/constrained-rescore/g12d_repeat_spread.json','w'),indent=2)
