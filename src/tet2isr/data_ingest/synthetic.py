"""Synthetic ground-truth world for the smoke test.

EVERYTHING written by this module is simulated.  Gene IDs are prefixed SYN so that no
synthetic number can be mistaken for a statement about a real gene.  The generator plants:
  * a shared, sequence-driven ISR dTE response (uORF-reinitiation up, TOP-like down);
  * modulation of that response by one basal cell-state program ("tc_availability");
  * a stressor-specific pleiotropic component (ER-like stress hits secretory mRNAs);
  * TET2-KO basal RNA changes, a cell-composition shift (single-cell study) and, unless
    ``null_world`` is set, a KO shift in the tc_availability program that creates a true I_TE.
A hidden WT/KO x ISR paired RNA/Ribo experiment is written to ``tet2_isr_holdout`` and is
only ever read by the validation step.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PROGRAMS = ["tc_availability", "ribosome_biogenesis", "proliferation", "myeloid_bias",
            "secretory", "isr_uorf_targets", "top_mrnas", "apoptosis_like"]


def _nb(rng, mu, disp=0.03):
    mu = np.maximum(mu, 1e-9)
    r = 1.0 / disp
    return rng.negative_binomial(r, r / (r + mu))


def _counts(rng, log2_abund, depth, disp=0.03):
    p = 2.0 ** log2_abund
    return _nb(rng, depth * p / p.sum(), disp)


def generate(out: Path, seed: int = 7, n_genes: int = 1500, null_world: bool = False,
             ko_tc_shift: float = -0.6) -> None:
    rng = np.random.default_rng(seed)
    out = Path(out)
    G = n_genes
    gid = np.array([f"SYN{i:05d}" for i in range(G)])

    # ---- genes: features, families, programs -------------------------------------------
    fam = np.arange(G)
    n_par = G // 10  # 10% of genes are paralogs that copy another gene's features
    par = rng.choice(G, n_par, replace=False)
    src = rng.choice(np.setdiff1d(np.arange(G), par), n_par)
    reinit = (rng.random(G) < 0.06).astype(float)
    top = ((rng.random(G) < 0.05) & (reinit == 0)).astype(float)
    n_uorf = rng.poisson(0.6, G).astype(float) + reinit
    utr5 = rng.normal(0, 1, G)
    gc = rng.normal(0, 1, G)
    cds = rng.normal(0, 1, G)
    for arr in (reinit, top, n_uorf, utr5, gc, cds):
        arr[par] = arr[src]
    fam[par] = fam[src]
    base = rng.normal(6, 2, G)  # log2 basal abundance
    genes = pd.DataFrame(dict(gene_id=gid, family_id=[f"FAM{f:05d}" for f in fam],
                              seq_uorf_reinit=reinit, seq_top=top, seq_n_uorf=n_uorf,
                              seq_utr5_len_z=utr5, seq_gc_z=gc, seq_cds_len_z=cds,
                              epi_promoter_acc=base + rng.normal(0, 1.5, G)))
    free = np.where((reinit == 0) & (top == 0))[0]
    rng.shuffle(free)
    sizes = dict(tc_availability=15, ribosome_biogenesis=40, proliferation=40, myeloid_bias=60,
                 secretory=120, apoptosis_like=30)
    prog, i = {}, 0
    for k, n in sizes.items():
        prog[k] = free[i:i + n]
        i += n
    prog["isr_uorf_targets"] = np.where(reinit == 1)[0]
    prog["top_mrnas"] = np.where(top == 1)[0]
    programs = pd.DataFrame([(k, gid[j]) for k, v in prog.items() for j in v],
                            columns=["program", "gene_id"])

    # ---- true translation rules ---------------------------------------------------------
    u = 1.5 * reinit - 0.6 * top + 0.10 * (n_uorf - reinit) - 0.10 * utr5 + rng.normal(0, 0.25, G)
    a = -0.7 * reinit + 0.25 * top + rng.normal(0, 0.05, G)
    u[par] = u[src] + rng.normal(0, 0.05, n_par)
    a[par] = a[src]
    pleio = np.zeros(G)
    pleio[prog["secretory"]] = -0.5

    def true_dte(intensity, z_tc, is_pleio, noise_sd=0.10):
        return intensity * (u + a * z_tc) + is_pleio * intensity * pleio + rng.normal(0, noise_sd, G)

    # ---- ISR compendium -----------------------------------------------------------------
    ctx_spec = [  # study, cell_type, stressor, class, kinase, intensity, z_tc, z_ribi, z_prolif
        ("SYN_ISR_A", "epithelial", "thapsigargin", "pleiotropic", "PERK", 1.0, 0.6, 0.3, 0.8),
        ("SYN_ISR_A", "epithelial", "arsenite", "canonical", "HRI", 1.4, 0.6, 0.3, 0.8),
        ("SYN_ISR_B", "fibroblast", "thapsigargin", "pleiotropic", "PERK", 0.8, -0.2, -0.5, -0.6),
        ("SYN_ISR_B", "fibroblast", "aa_starvation", "canonical", "GCN2", 0.7, -0.2, -0.5, -0.6),
        ("SYN_ISR_C", "esc", "arsenite", "canonical", "HRI", 1.2, 1.0, 0.9, 1.0),
        ("SYN_ISR_C", "esc", "aa_starvation", "canonical", "GCN2", 0.6, 1.0, 0.9, 1.0),
        ("SYN_ISR_D", "myeloid_line", "thapsigargin", "pleiotropic", "PERK", 0.5, -0.8, 0.1, 0.4),
        ("SYN_ISR_D", "myeloid_line", "thapsigargin_hi", "pleiotropic", "PERK", 1.5, -0.8, 0.1, 0.4),
        ("SYN_ISR_E", "erythroid_line", "hri_activator", "canonical", "HRI", 1.0, -0.5, 0.6, 0.2),
        ("SYN_ISR_F", "neuronal", "aa_starvation", "canonical", "GCN2", 0.9, 0.2, -0.8, -1.0),
        ("SYN_ISR_F", "neuronal", "pkr_activator", "canonical", "PKR", 1.1, 0.2, -0.8, -1.0),
    ]
    ct_base, ct_te0, contexts, samples, cols, truth = {}, {}, [], [], {}, {}
    for (study, ct, stressor, cls, kin, inten, ztc, zr, zp) in ctx_spec:
        if ct not in ct_base:
            b = base + rng.normal(0, 0.7, G)
            b[prog["tc_availability"]] += ztc
            b[prog["ribosome_biogenesis"]] += zr
            b[prog["proliferation"]] += zp
            ct_base[ct] = b
        cid = f"{study}|{ct}|{stressor}"
        contexts.append(dict(context_id=cid, study_id=study, species="synthetic", cell_type=ct,
                             stressor=stressor, stressor_class=cls, kinase=kin, intensity=inten,
                             duration_h=2.0))
        d_te = true_dte(inten, ztc, cls == "pleiotropic")
        d_rna = inten * 0.8 * rng.normal(0, 0.15, G)  # early time point: little RNA change
        truth[cid] = d_te
        te0 = ct_te0.setdefault(ct, rng.normal(0, 0.8, G))  # basal TE is cell-type specific
        for cond in ("control", "stress"):
            for rep in range(1, 4):
                bio = rng.normal(0, 0.12, G)
                rna = ct_base[ct] + bio + (d_rna if cond == "stress" else 0)
                ribo = rna + te0 + (d_te if cond == "stress" else 0) + rng.normal(0, 0.08, G)
                for assay, lg, depth in (("rna", rna, 4e6), ("ribo", ribo, 2e6)):
                    sid = f"{cid}|{cond}|{assay}|r{rep}"
                    cols[sid] = _counts(rng, lg, depth)
                    samples.append(dict(sample_id=sid, context_id=cid, study_id=study,
                                        condition=cond, assay=assay, bio_rep=f"{cid}|{cond}|r{rep}"))

    (out / "isr").mkdir(parents=True, exist_ok=True)
    (out / "truth").mkdir(exist_ok=True)
    genes.to_parquet(out / "genes.parquet", index=False)
    programs.to_csv(out / "programs.tsv", sep="\t", index=False)
    pd.DataFrame(contexts).to_csv(out / "isr/contexts.tsv", sep="\t", index=False)
    pd.DataFrame(samples).to_csv(out / "isr/samples.tsv", sep="\t", index=False)
    pd.DataFrame(cols, index=pd.Index(gid, name="gene_id")).to_parquet(out / "isr/counts.parquet")
    pd.DataFrame(truth, index=pd.Index(gid, name="gene_id")).to_parquet(out / "truth/dte_by_context.parquet")

    # ---- TET2 basal transcriptomes ------------------------------------------------------
    ko_fc = rng.normal(0, 0.08, G)
    ko_fc[prog["myeloid_bias"]] += 0.9
    ko_fc[prog["proliferation"]] += 0.4
    ko_fc[prog["apoptosis_like"]] -= 0.5
    sub = rng.choice(prog["isr_uorf_targets"], len(prog["isr_uorf_targets"]) // 3, replace=False)
    ko_fc[sub] += 0.8  # some ISR targets are more abundant in KO: pure substrate effect
    shift = 0.0 if null_world else ko_tc_shift
    ko_fc[prog["tc_availability"]] += shift
    hspc = base + rng.normal(0, 0.7, G)
    hspc[prog["tc_availability"]] += -0.3
    mpp_delta = np.zeros(G)
    mpp_delta[prog["proliferation"]] += 1.0
    mpp_delta[prog["myeloid_bias"]] += 0.6

    # study 1: bulk, "mouse-like" IDs needing an ortholog join
    d1 = out / "tet2/SYN_TET2_bulk"
    d1.mkdir(parents=True, exist_ok=True)
    mid = np.array([f"MSYN{i:05d}" for i in range(G)])
    fc1 = ko_fc + rng.normal(0, 0.15, G)
    smp, cols = [], {}
    for geno in ("WT", "TET2"):
        for r in range(1, 5):
            sid = f"bulk_{geno}_{r}"
            cols[sid] = _counts(rng, hspc + (fc1 if geno == "TET2" else 0) + rng.normal(0, 0.15, G), 5e6)
            smp.append(dict(sample_id=sid, study_id="SYN_TET2_bulk", species="synthetic_mouse",
                            genotype=geno, allele_type="engineered_KO", compartment="HSC",
                            donor=f"animal_{geno}_{r}"))
    pd.DataFrame(smp).to_csv(d1 / "samples.tsv", sep="\t", index=False)
    pd.DataFrame(cols, index=pd.Index(mid, name="gene_id")).to_parquet(d1 / "counts.parquet")
    omap = pd.DataFrame(dict(source_gene_id=mid, gene_id=gid))
    dup = omap.sample(40, random_state=seed).assign(gene_id=lambda d: rng.choice(gid, 40))
    pd.concat([omap, dup]).to_csv(d1 / "orthologs.tsv", sep="\t", index=False)  # 40 one-to-many

    # study 2: single cell, two compartments, KO shifts composition toward MPP
    d2 = out / "tet2/SYN_TET2_sc"
    d2.mkdir(parents=True, exist_ok=True)
    fc2 = ko_fc + rng.normal(0, 0.15, G)
    smp, cells, cols = [], [], {}
    for geno in ("WT", "TET2"):
        for r in range(1, 4):
            sid = f"sc_{geno}_{r}"
            donor = rng.normal(0, 0.15, G)
            smp.append(dict(sample_id=sid, study_id="SYN_TET2_sc", species="synthetic_human",
                            genotype=geno, allele_type="CRISPR_edit", compartment="mixed",
                            donor=f"donor_{geno}_{r}"))
            p_mpp = 0.65 if geno == "TET2" else 0.40
            for c in range(120):
                comp = "MPP" if rng.random() < p_mpp else "HSC"
                lg = hspc + donor + (fc2 if geno == "TET2" else 0) + (mpp_delta if comp == "MPP" else 0)
                cid = f"{sid}_c{c:03d}"
                cols[cid] = _counts(rng, lg, 2e4, disp=0.2)
                cells.append(dict(cell_id=cid, sample_id=sid, compartment=comp))
    pd.DataFrame(smp).to_csv(d2 / "samples.tsv", sep="\t", index=False)
    pd.DataFrame(cells).to_csv(d2 / "cells.tsv", sep="\t", index=False)
    pd.DataFrame(cols, index=pd.Index(gid, name="gene_id")).astype("int32").to_parquet(d2 / "counts.parquet")

    # ---- hidden genotype x ISR translatome (validation only) ---------------------------
    d3 = out / "tet2_isr_holdout"
    d3.mkdir(exist_ok=True)
    inten, z_wt = 1.0, -0.3
    te0 = rng.normal(0, 0.8, G)
    dte_k = {"WT": true_dte(inten, z_wt, False), "TET2": true_dte(inten, z_wt + shift, False)}
    smp, cols = [], {}
    for geno in ("WT", "TET2"):
        for cond in ("control", "stress"):
            for r in range(1, 4):
                rna = hspc + (ko_fc if geno == "TET2" else 0) + rng.normal(0, 0.12, G)
                ribo = rna + te0 + (dte_k[geno] if cond == "stress" else 0) + rng.normal(0, 0.08, G)
                for assay, lg, depth in (("rna", rna, 4e6), ("ribo", ribo, 2e6)):
                    sid = f"hold_{geno}_{cond}_{assay}_{r}"
                    cols[sid] = _counts(rng, lg, depth)
                    smp.append(dict(sample_id=sid, genotype=geno, condition=cond, assay=assay,
                                    bio_rep=f"{geno}_{cond}_{r}"))
    pd.DataFrame(smp).to_csv(d3 / "samples.tsv", sep="\t", index=False)
    pd.DataFrame(cols, index=pd.Index(gid, name="gene_id")).to_parquet(d3 / "counts.parquet")
    pd.DataFrame(dict(gene_id=gid, true_dTE_WT=dte_k["WT"], true_dTE_TET2=dte_k["TET2"],
                      true_I_TE=inten * a * shift, true_u=u, true_a=a,
                      true_ko_log2fc=ko_fc)).to_parquet(out / "truth/genotype_truth.parquet", index=False)
    (out / "SYNTHETIC_DATA_README.txt").write_text(
        f"All files in this directory are SIMULATED (seed={seed}, null_world={null_world}).\n"
        "They demonstrate that the pipeline runs and recovers planted effects. They contain\n"
        "no biological measurements and support no conclusion about TET2 or the ISR.\n")

