import os, json
import pandas as pd

root = "/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/CxC_TxT/cxc_results/single_mutations/"
records = []

for d in os.listdir(root):
    folder = os.path.join(root, d)
    if not os.path.isdir(folder): 
        continue
    mut = d.replace("fold_cxc_", "").upper()  # e.g. '178D'
    confs = []
    for f in os.listdir(folder):
        if f.endswith("summary_confidences_0.json"):
            path = os.path.join(folder, f)
            with open(path) as jf:
                data = json.load(jf)
                iptm = data.get("iptm", None)
                ptm = data.get("ptm", None)
                rank = data.get("ranking_confidence", None)
                confs.append(rank or iptm or ptm)
    if confs:
        records.append({
            "Combo": mut,
            "fitness": sum(confs)/len(confs)
        })

df = pd.DataFrame(records)
df.to_csv("/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/CxC_TxT/cxc_results/single_mutations/CXC_alde_input.csv", index=False)
print(df)