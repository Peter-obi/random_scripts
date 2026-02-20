import json, pandas as pd

df = pd.read_csv("prb100_v1.csv")
manifest = {
  "n": int(len(df)),
  "bins": df["prb_bin"].value_counts().sort_index().to_dict(),
  "charge_counts": {
    "neg": int((df["formal_charge"]<0).sum()),
    "neu": int((df["formal_charge"]==0).sum()),
    "pos": int((df["formal_charge"]>0).sum()),
    "charged_total": int((df["formal_charge"]!=0).sum()),
  },
  "suite_applicability": {
    "A_elec": int(df["A_elec"].sum()),
    "A_hb": int(df["A_hb"].sum()),
    "A_desolv": int(df["A_desolv"].sum()),
  },
  "negatives_per_bin": df.groupby("prb_bin")["formal_charge"].apply(lambda s: int((s<0).sum())).sort_index().to_dict(),
  "charged_per_bin": df.groupby("prb_bin")["is_charged"].sum().astype(int).sort_index().to_dict(),
}
with open("prb100_v1_manifest.json","w") as f:
    json.dump(manifest, f, indent=2)
print("wrote prb100_v1_manifest.json")

