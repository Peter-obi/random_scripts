#!/usr/bin/env python3
"""
PRB-100 sampler (balanced across SuCOS-pocket-qcov bins) with:

1) Suite applicability balance (per bin + global):
   - A_hb (HB geometry meaningful)
   - A_desolv (desolvation traps meaningful)
   - A_elec (electrostatic sensitivity meaningful; "practical" allows polar neutrals)

2) HARD charge quotas (to avoid A_elec==everyone degeneracy):
   - per-bin: at least BIN_CHARGED_MIN charged ligands if available
   - global: at least GLOBAL_CHARGED_MIN charged ligands
   - global: at least GLOBAL_NEG_MIN anions (best effort)

This solves your loop:
- Use elec_mode=practical so Suite-1 coverage is everywhere
- Still force a real charged subset so Suite-1 isn't just "polarity"

Usage:
  python prb100_sampler.py --in annotations.csv --out prb100_v1.csv \
    --elec_mode practical \
    --bin_elec_min 8 --bin_hb_min 7 --bin_des_min 6 \
    --bin_charged_min 2 \
    --global_elec_min 80 --global_hb_min 60 --global_des_min 60 \
    --global_charged_min 30 --global_neg_min 8 \
    --ccd_cap 2000

If a bin can't meet charged/anionic minima (rare anions), it falls back gracefully
but still fills bins.
"""

from __future__ import annotations

import argparse
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

RNG_SEED_DEFAULT = 20260214

BIN_TARGETS = {
    "0-20": 12,
    "20-30": 13,
    "30-40": 12,
    "40-50": 13,
    "50-60": 12,
    "60-70": 13,
    "70-80": 12,
    "80-100": 13,
}

METAL_TOKENS = ["[Zn", "[Mg", "[Mn", "[Fe", "[Co", "[Ni", "[Cu", "[Ca", "[Na", "[K", "[Al", "[Ag", "[Au"]


def assign_bin(x: float) -> str:
    x = float(x)
    if x < 20: return "0-20"
    if x < 30: return "20-30"
    if x < 40: return "30-40"
    if x < 50: return "40-50"
    if x < 60: return "50-60"
    if x < 70: return "60-70"
    if x < 80: return "70-80"
    return "80-100"


def hac_bucket(hac: int) -> Optional[str]:
    if 10 <= hac <= 20: return "frag"
    if 21 <= hac <= 45: return "drug"
    if 46 <= hac <= 80: return "large"
    return None


def parse_mol(smiles: str) -> Optional[Chem.Mol]:
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    if any(tok in smiles for tok in METAL_TOKENS):
        return None
    return Chem.MolFromSmiles(smiles)


def mol_features(smiles: str) -> Optional[dict]:
    mol = parse_mol(smiles)
    if mol is None:
        return None
    q = Chem.GetFormalCharge(mol)
    hbd = rdMolDescriptors.CalcNumLipinskiHBD(mol)
    hba = rdMolDescriptors.CalcNumLipinskiHBA(mol)
    hb_sum = int(hbd + hba)
    return {"formal_charge": int(q), "hbd": int(hbd), "hba": int(hba), "hb_sum": hb_sum}


def charge_bucket(q: int) -> str:
    if q < 0: return "neg"
    if q > 0: return "pos"
    return "neu"


def applicability_flags(formal_charge: int, hb_sum: int, tpsa: float, elec_mode: str) -> Tuple[bool, bool, bool]:
    # Electrostatics applicability
    if elec_mode == "strict":
        A_elec = abs(formal_charge) >= 1
    elif elec_mode == "practical":
        A_elec = (abs(formal_charge) >= 1) or (tpsa >= 120.0) or (hb_sum >= 8)
    else:
        raise ValueError(f"Unknown elec_mode={elec_mode}")

    # HB geometry applicability
    A_hb = hb_sum >= 4

    # Desolvation applicability
    A_desolv = (tpsa >= 80.0) or (hb_sum >= 6)

    return A_elec, A_hb, A_desolv


@dataclass
class GlobalCounts:
    elec: int = 0
    hb: int = 0
    des: int = 0
    charged: int = 0
    neg: int = 0

    def add_row(self, row: pd.Series) -> None:
        self.elec += int(row["A_elec"])
        self.hb += int(row["A_hb"])
        self.des += int(row["A_desolv"])
        q = int(row["formal_charge"])
        if q != 0:
            self.charged += 1
        if q < 0:
            self.neg += 1

    def score_row(self, row: pd.Series, mins: Dict[str, int]) -> int:
        """
        Ranking score to greedily satisfy global minima (never blocks filling bins).
        """
        s = 0
        if self.elec < mins["elec"] and row["A_elec"]:
            s += 3
        if self.hb < mins["hb"] and row["A_hb"]:
            s += 2
        if self.des < mins["des"] and row["A_desolv"]:
            s += 2

        q = int(row["formal_charge"])
        if self.charged < mins["charged"] and q != 0:
            s += 8
        if self.neg < mins["neg"] and q < 0:
            s += 10

        # prefer rare CCD if available
        if "num_training_systems_with_similar_ccds" in row.index and pd.notna(row["num_training_systems_with_similar_ccds"]):
            x = int(row["num_training_systems_with_similar_ccds"])
            if x == 0:
                s += 1
            elif x > 1000:
                s -= 1
        return s


def sample_k(pool: pd.DataFrame, k: int, rng: random.Random, score_fn=None) -> pd.DataFrame:
    if k <= 0:
        return pool.iloc[0:0]
    if len(pool) <= k:
        return pool.copy()
    if score_fn is None:
        return pool.sample(n=k, random_state=rng.randint(0, 2**31 - 1))
    tmp = pool.copy()
    tmp["_score"] = tmp.apply(score_fn, axis=1)
    tmp = tmp.sort_values(by="_score", ascending=False)
    return tmp.head(k).drop(columns=["_score"])


def fill_bin(pool: pd.DataFrame,
             n_target: int,
             bin_mins: Dict[str, int],
             global_counts: GlobalCounts,
             global_mins: Dict[str, int],
             rng: random.Random) -> pd.DataFrame:
    """
    Phase 1: satisfy per-bin minima (best effort), including charged
    Phase 2: fill remainder using global scoring
    Phase 3: top-up random
    """
    used = set()
    selected_rows = []

    def take_mask(mask_col: str, need: int, scorer=None) -> None:
        nonlocal used, selected_rows
        if need <= 0:
            return
        cand = pool[(pool[mask_col] == True) & (~pool["system_id"].isin(used))].copy()
        if cand.empty:
            return
        picks = sample_k(cand, need, rng, score_fn=scorer)
        used |= set(picks["system_id"].tolist())
        for _, r in picks.iterrows():
            selected_rows.append(r)
            global_counts.add_row(r)

    # Prefer multi-applicable + helps global mins
    def multi_score(row: pd.Series) -> int:
        s = int(row["A_elec"]) + int(row["A_hb"]) + int(row["A_desolv"])
        s += global_counts.score_row(row, global_mins)
        return s

    # Per-bin minima (best effort)
    take_mask("is_charged", bin_mins["charged"], scorer=multi_score)
    take_mask("A_elec",      bin_mins["elec"],    scorer=multi_score)
    take_mask("A_hb",        bin_mins["hb"],      scorer=multi_score)
    take_mask("A_desolv",    bin_mins["des"],     scorer=multi_score)

    # Fill remainder by global needs
    remaining = n_target - len(selected_rows)
    if remaining > 0:
        cand = pool[~pool["system_id"].isin(used)].copy()

        def global_score(row: pd.Series) -> int:
            return global_counts.score_row(row, global_mins)

        picks = sample_k(cand, remaining, rng, score_fn=global_score)
        used |= set(picks["system_id"].tolist())
        for _, r in picks.iterrows():
            selected_rows.append(r)
            global_counts.add_row(r)

    # Top-up random if still short
    if len(selected_rows) < n_target:
        remaining = n_target - len(selected_rows)
        cand = pool[~pool["system_id"].isin(used)].copy()
        picks = sample_k(cand, remaining, rng)
        for _, r in picks.iterrows():
            selected_rows.append(r)
            global_counts.add_row(r)

    out = pd.DataFrame(selected_rows).drop_duplicates(subset=["system_id"]).reset_index(drop=True)
    return out.head(n_target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_csv", default="annotations.csv")
    ap.add_argument("--out", dest="out_csv", default="prb100_v1.csv")
    ap.add_argument("--seed", type=int, default=RNG_SEED_DEFAULT)
    ap.add_argument("--ccd_cap", type=int, default=-1)

    ap.add_argument("--elec_mode", type=str, default="practical", choices=["strict", "practical"])

    # per-bin suite minima
    ap.add_argument("--bin_elec_min", type=int, default=8)
    ap.add_argument("--bin_hb_min", type=int, default=7)
    ap.add_argument("--bin_des_min", type=int, default=6)
    ap.add_argument("--bin_charged_min", type=int, default=2)

    # global minima
    ap.add_argument("--global_elec_min", type=int, default=80)
    ap.add_argument("--global_hb_min", type=int, default=60)
    ap.add_argument("--global_des_min", type=int, default=60)
    ap.add_argument("--global_charged_min", type=int, default=30)
    ap.add_argument("--global_neg_min", type=int, default=8)

    args = ap.parse_args()
    rng = random.Random(args.seed)

    df = pd.read_csv(args.in_csv)

    required = [
        "system_id",
        "ligand_smiles",
        "ligand_is_proper",
        "num_proper_ligand_chains",
        "ligand_num_heavy_atoms",
        "ligand_tpsa",
        "sucos_shape_pocket_qcov",
    ]
    for c in required:
        if c not in df.columns:
            raise RuntimeError(f"Missing required column: {c}")

    # Hard filters (only what your metrics need)
    df = df[df["ligand_is_proper"] == True].copy()
    df = df[df["num_proper_ligand_chains"] == 1].copy()
    df = df[df["ligand_smiles"].notna()].copy()
    df = df[df["ligand_tpsa"].notna()].copy()
    df = df[df["sucos_shape_pocket_qcov"].notna()].copy()

    # HAC window
    df["hac"] = df["ligand_num_heavy_atoms"].astype(int)
    df["hac_bucket"] = df["hac"].apply(hac_bucket)
    df = df[df["hac_bucket"].notna()].copy()

    # Optional CCD cap
    if args.ccd_cap is not None and args.ccd_cap >= 0 and "num_training_systems_with_similar_ccds" in df.columns:
        df["num_training_systems_with_similar_ccds"] = df["num_training_systems_with_similar_ccds"].fillna(0).astype(int)
        df = df[df["num_training_systems_with_similar_ccds"] <= args.ccd_cap].copy()

    # RDKit features
    feats = []
    keep_idx = []
    for idx, smi in zip(df.index, df["ligand_smiles"].tolist()):
        f = mol_features(smi)
        if f is None:
            continue
        feats.append(f)
        keep_idx.append(idx)

    df = df.loc[keep_idx].copy()
    feat_df = pd.DataFrame(feats, index=df.index)
    df = pd.concat([df, feat_df], axis=1)

    df["charge_bucket"] = df["formal_charge"].apply(charge_bucket)
    df["is_charged"] = df["formal_charge"].apply(lambda q: q != 0)

    # Applicability flags
    A_elec, A_hb, A_des = [], [], []
    for _, row in df.iterrows():
        e, h, d = applicability_flags(
            formal_charge=int(row["formal_charge"]),
            hb_sum=int(row["hb_sum"]),
            tpsa=float(row["ligand_tpsa"]),
            elec_mode=args.elec_mode,
        )
        A_elec.append(e)
        A_hb.append(h)
        A_des.append(d)
    df["A_elec"] = A_elec
    df["A_hb"] = A_hb
    df["A_desolv"] = A_des

    # Bins
    df["prb_bin"] = df["sucos_shape_pocket_qcov"].apply(assign_bin)

    # Print feasibility
    print("\nAvailability per bin (after filters):")
    print(df["prb_bin"].value_counts().sort_index())

    print("\nApplicability availability per bin:")
    by = df.groupby("prb_bin")[["A_elec", "A_hb", "A_desolv", "is_charged"]].sum().astype(int)
    by["n"] = df.groupby("prb_bin").size()
    print(by.sort_index())

    # Selection
    bin_mins = {"elec": args.bin_elec_min, "hb": args.bin_hb_min, "des": args.bin_des_min, "charged": args.bin_charged_min}
    global_mins = {"elec": args.global_elec_min, "hb": args.global_hb_min, "des": args.global_des_min,
                   "charged": args.global_charged_min, "neg": args.global_neg_min}

    global_counts = GlobalCounts()
    used_ids = set()
    parts: List[pd.DataFrame] = []

    for b, n in BIN_TARGETS.items():
        pool = df[df["prb_bin"] == b].copy()
        if pool.empty:
            raise RuntimeError(f"Bin {b} is empty after filters.")

        pool = pool[~pool["system_id"].isin(used_ids)].copy()
        pool = pool.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)

        sel = fill_bin(pool, n, bin_mins, global_counts, global_mins, rng)
        if len(sel) < n:
            raise RuntimeError(f"Could not fill bin {b}: needed {n}, got {len(sel)} (pool {len(pool)})")

        used_ids |= set(sel["system_id"].tolist())
        parts.append(sel)

    out = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["system_id"]).reset_index(drop=True)
    if len(out) != 100:
        raise RuntimeError(f"Final size {len(out)} != 100")

    # Save compact columns
    keep_cols = [
        "system_id",
        "entry_pdb_id",
        "ligand_instance_chain",
        "ligand_smiles",
        "sucos_shape_pocket_qcov",
        "prb_bin",
        "formal_charge",
        "charge_bucket",
        "is_charged",
        "hbd",
        "hba",
        "hb_sum",
        "ligand_tpsa",
        "hac",
        "hac_bucket",
        "A_elec",
        "A_hb",
        "A_desolv",
        "ligand_molecular_weight",
        "ligand_num_rot_bonds",
        "ligand_num_rings",
        "ligand_num_pocket_residues",
        "num_training_systems_with_similar_ccds",
        "cluster",
    ]
    keep_cols = [c for c in keep_cols if c in out.columns]
    out[keep_cols].to_csv(args.out_csv, index=False)

    # Summary
    print("\nWrote:", args.out_csv)
    print("\nprb_bin")
    print(out["prb_bin"].value_counts().sort_index())
    print("\nCharge counts:",
          int((out["formal_charge"] < 0).sum()),
          int((out["formal_charge"] == 0).sum()),
          int((out["formal_charge"] > 0).sum()))
    print("\nCharged total:", int(out["is_charged"].sum()))
    print("\nApplicability totals:",
          "A_elec=", int(out["A_elec"].sum()),
          "A_hb=", int(out["A_hb"].sum()),
          "A_desolv=", int(out["A_desolv"].sum()))
    print("\nCharged per bin:")
    print(out.groupby("prb_bin")["is_charged"].sum().astype(int).sort_index())
    print("\nNegatives per bin:")
    print(out.groupby("prb_bin")["formal_charge"].apply(lambda s: int((s < 0).sum())).sort_index())


if __name__ == "__main__":
    main()

