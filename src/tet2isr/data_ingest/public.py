"""Adapters: downloaded GEO processed files -> the data contract (data/processed/public).

Common gene space: human Ensembl gene IDs (unversioned), protein-coding genes with an
Ensembl-canonical transcript. Mouse data are mapped through one-to-one MGI homology classes;
ambiguous classes are dropped and counted.  Every sample-to-condition assignment below was
read from the GEO sample records saved by `make audit-data` (results/audit/).
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .orthologs import one_to_one

MIN_DEPTH = 3e5
RAW = Path("data/raw")
REF = RAW / "reference"

# Curated sets. Symbols recalled from the cited papers, not copied from their tables: the data
# check in results/public/reference_set_check.csv shows whether each behaves as expected.
CURATED = {
    # translationally induced / eIF2a-P-resistant mRNAs (Sidrauski 2015 PMID:25719440; Andreev 2015 PMID:25621764)
    "isr_uorf_targets": ["ATF4", "ATF5", "DDIT3", "PPP1R15A", "PPP1R15B", "IFRD1", "SLC35A4", "UCP2", "PTP4A1"],
    "eif2_ternary_complex": ["EIF2S1", "EIF2S2", "EIF2S3", "EIF2B1", "EIF2B2", "EIF2B3", "EIF2B4", "EIF2B5",
                             "EIF5", "EIF1", "EIF1AX"],
    "atf4_targets": ["ASNS", "PSAT1", "PHGDH", "PSPH", "SHMT2", "MTHFD2", "SLC7A11", "SLC7A5", "SLC3A2",
                     "SLC7A1", "CHAC1", "TRIB3", "DDIT4", "SESN2", "ATF3", "VEGFA", "XPOT", "STC2", "CEBPB",
                     "NUPR1", "ALDH1L2", "PCK2", "CTH", "SLC1A5", "GPT2", "PYCR1"],
    "bcl2_family": ["BCL2", "MCL1", "BCL2L1", "BCL2L2", "BCL2A1", "BCL2L11", "BBC3", "PMAIP1", "BAX", "BAK1",
                    "BID", "BAD", "BMF", "HRK", "BIK"],
}
HALLMARKS = ["E2F_TARGETS", "G2M_CHECKPOINT", "MYC_TARGETS_V1", "MTORC1_SIGNALING", "UNFOLDED_PROTEIN_RESPONSE",
             "HEME_METABOLISM", "APOPTOSIS", "P53_PATHWAY", "OXIDATIVE_PHOSPHORYLATION",
             "INTERFERON_ALPHA_RESPONSE", "TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE"]


# ------------------------------------------------------------------ reference
def build_genes() -> pd.DataFrame:
    """Sequence features of the Ensembl-canonical transcript of each protein-coding gene."""
    tx = {}
    with gzip.open(REF / "Homo_sapiens.GRCh38.113.gtf.gz", "rt") as fh:
        for line in fh:
            if "Ensembl_canonical" not in line or 'transcript_biotype "protein_coding"' not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("exon", "CDS", "stop_codon"):
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            t = tx.setdefault(a["transcript_id"], dict(gene_id=a["gene_id"], symbol=a.get("gene_name", ""),
                                                       strand=f[6], exon=[], cds=[]))
            t["exon" if f[2] == "exon" else "cds"].append((int(f[3]), int(f[4])))
    seqs, name = {}, None
    with gzip.open(REF / "Homo_sapiens.GRCh38.cdna.all.fa.gz", "rt") as fh:
        for line in fh:
            if line[0] == ">":
                name = line[1:].split(".")[0]
                name = name if name in tx else None
                if name:
                    seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    rows = []
    for tid, t in tx.items():
        if not t["cds"] or tid not in seqs or not t["symbol"]:
            continue
        seq = "".join(seqs[tid]).upper()
        ex = sorted(t["exon"], reverse=t["strand"] == "-")
        lo, hi = min(c[0] for c in t["cds"]), max(c[1] for c in t["cds"])
        first = hi if t["strand"] == "-" else lo            # genomic coordinate of first coding nt
        last = lo if t["strand"] == "-" else hi
        pos, start, end = 0, None, None
        for a, b in ex:
            if a <= first <= b:
                start = pos + (b - first if t["strand"] == "-" else first - a)
            if a <= last <= b:
                end = pos + (b - last if t["strand"] == "-" else last - a) + 1
            pos += b - a + 1
        if start is None or end is None or end <= start or pos != len(seq):
            continue
        u5, cds, u3 = seq[:start], seq[start:end], seq[end:]
        gc = lambda s: (s.count("G") + s.count("C")) / len(s) if s else np.nan
        n_uaug = n_uorf = n_ouorf = 0
        for m in re.finditer("(?=ATG)", u5):
            n_uaug += 1
            i = m.start()
            stop = next((j for j in range(i, len(seq) - 2, 3) if seq[j:j + 3] in ("TAA", "TAG", "TGA")), None)
            if stop is not None and stop + 3 <= start:
                n_uorf += 1                                   # uORF terminating inside the 5'UTR
            elif (start - i) % 3:
                n_ouorf += 1                                  # out-of-frame, overlapping the main start
        rows.append(dict(gene_id=t["gene_id"], symbol=t["symbol"], transcript_id=tid,
                         seq_log_utr5_len=np.log10(len(u5) + 1), seq_utr5_gc=gc(u5),
                         seq_log_n_uaug=np.log2(n_uaug + 1), seq_log_n_uorf=np.log2(n_uorf + 1),
                         seq_overlapping_uorf=float(n_ouorf > 0),
                         seq_top_like_start=float(bool(re.match("C[CT]{4}", u5))),
                         seq_log_cds_len=np.log10(len(cds)), seq_cds_gc=gc(cds),
                         seq_log_utr3_len=np.log10(len(u3) + 1)))
    g = pd.DataFrame(rows).drop_duplicates("gene_id").drop_duplicates("symbol")
    g["seq_utr5_gc"] = g.seq_utr5_gc.fillna(g.seq_utr5_gc.mean())
    # heuristic family: HGNC root symbol (trailing digits / letter suffix removed). Not a paralogy call.
    g["family_id"] = g.symbol.str.replace(r"^([A-Z]{2,}?)\d+[A-Z]{0,2}\d*$", r"\1", regex=True)
    return g.reset_index(drop=True)


def build_programs(genes: pd.DataFrame) -> pd.DataFrame:
    sym = genes.set_index("symbol").gene_id
    sets = {k: v for k, v in CURATED.items()}
    sets["ribosomal_proteins"] = [s for s in sym.index if re.fullmatch(r"RP[LS]\d+[AXY]?\d?|RPLP[012]|RPSA", s)]
    sets["aa_trna_synthetases"] = [s for s in sym.index if re.fullmatch(r"[A-Z]{1,2}ARS1", s)]
    for line in open(REF / "h.all.v2024.1.Hs.symbols.gmt"):
        f = line.rstrip("\n").split("\t")
        if f[0].replace("HALLMARK_", "") in HALLMARKS:
            sets[f[0]] = f[2:]
    rows = [(p, sym[s]) for p, ss in sets.items() for s in ss if s in sym.index]
    return pd.DataFrame(rows, columns=["program", "gene_id"])


def mouse_to_human(genes: pd.DataFrame, key: str = "symbol") -> tuple[pd.Series, dict]:
    """Mouse symbol (or Entrez ID) -> human ENSG via one-to-one MGI homology classes."""
    h = pd.read_csv(REF / "HOM_MouseHumanSequence.rpt", sep="\t", dtype=str)
    col = "Symbol" if key == "symbol" else "EntrezGene ID"
    mm = h[h["NCBI Taxon ID"] == "10090"][["DB Class Key", col]].rename(columns={col: "source_gene_id"})
    hs = h[h["NCBI Taxon ID"] == "9606"][["DB Class Key", "Symbol"]]
    pair = mm.merge(hs, on="DB Class Key").merge(genes[["symbol", "gene_id"]], left_on="Symbol", right_on="symbol")
    kept, info = one_to_one(pair)
    return kept.set_index("source_gene_id").gene_id, info


def _refseq_to_gene(genes: pd.DataFrame) -> pd.Series:
    """Human RefSeq transcript accession (unversioned) -> ENSG, from the MGI homology report."""
    hom = pd.read_csv(REF / "HOM_MouseHumanSequence.rpt", sep="\t", dtype=str)
    hom = hom[hom["NCBI Taxon ID"] == "9606"][["Symbol", "Nucleotide RefSeq IDs"]].dropna()
    rs = hom.assign(rs=hom["Nucleotide RefSeq IDs"].str.split(",")).explode("rs")
    rs = rs.merge(genes[["symbol", "gene_id"]], left_on="Symbol", right_on="symbol").drop_duplicates("rs")
    return rs.set_index("rs").gene_id


def _to_gene_space(df: pd.DataFrame, mapping: pd.Series, universe: pd.Index) -> pd.DataFrame:
    """Rows of df (indexed by a source ID) -> summed per target gene, reindexed to the universe."""
    df = df[df.index.isin(mapping.index)]
    out = df.groupby(mapping.reindex(df.index).values).sum()
    return out.reindex(universe).fillna(0).round().astype("int64")


# ------------------------------------------------------------------ ISR studies
def _isr_studies(genes: pd.DataFrame):
    """Yield (context dicts, sample rows, counts[gene_id x sample]) per study."""
    uni = pd.Index(genes.gene_id)
    ens = pd.Series(uni, index=uni)
    sym = genes.set_index("symbol").gene_id
    strip = lambda idx: idx.str.replace(r"\..*$", "", regex=True)
    m2h_sym, _ = mouse_to_human(genes, "symbol")
    m2h_ent, _ = mouse_to_human(genes, "entrez")

    def pack(study, cell, species, arms, counts, extra=None):
        """arms: {context_key: dict(meta..., control={assay:[cols]}, stress={assay:[cols]})}"""
        ctxs, smps, keep = [], [], {}
        for key, a in arms.items():
            cid = f"{study}|{cell}|{key}"
            ctxs.append(dict(context_id=cid, study_id=study, species=species, cell_type=cell,
                             stressor=a["stressor"], stressor_class=a["cls"], kinase=a["kinase"],
                             intensity=1.0, duration_h=a["h"], dose=a["dose"], role=a.get("role", "train"),
                             note=a.get("note", "")))
            for cond in ("control", "stress"):
                for assay, cols in a[cond].items():
                    for i, c in enumerate(cols, 1):
                        sid = f"{cid}|{cond}|{assay}|{c}"
                        keep[sid] = counts[c]
                        smps.append(dict(sample_id=sid, context_id=cid, study_id=study, condition=cond,
                                         assay=assay, bio_rep=f"{cid}|{cond}|{assay}|{i}", source_column=c))
        return ctxs, smps, pd.DataFrame(keep)

    # --- GSE65778: HEK293T, tunicamycin 1 ug/ml 1 h +/- ISRIB; legacy UCSC transcript IDs
    d = RAW / "GSE65778/extracted"
    # IDs are matched on the version-less accession (mRNA and footprint files use different UCSC
    # releases). mRNA files repeat one count across isoforms, so the gene value is the isoform maximum.
    kg = pd.concat([pd.read_csv(REF / f"kgXrefOld{i}.txt.gz", sep="\t", header=None, usecols=[0, 4, 5],
                                names=["uc", "symbol", "refseq"]) for i in (5, 6)])
    kg["uc"] = kg.uc.str.replace(r"\.\d+$", "", regex=True)
    kg = kg.drop_duplicates("uc")
    # current symbol first; renamed genes (e.g. AARS -> AARS1) are recovered through their RefSeq accession
    kg["gene_id"] = kg.symbol.map(sym).fillna(kg.refseq.map(_refseq_to_gene(genes)))
    uc = kg.dropna(subset=["gene_id"]).set_index("uc").gene_id
    cols = {}
    for f in sorted(d.glob("*.txt.gz")):
        t = pd.read_csv(f, sep="\t", header=None, names=["uc", "len", "n"], na_values="N/A")
        t["g"] = t.uc.str.replace(r"\.\d+$", "", regex=True).map(uc)
        cols[re.sub(r"^GSM\d+_|_hek|\.txt\.gz", "", f.name)] = t.dropna(subset=["g"]).groupby("g").n.max()
    c = pd.DataFrame(cols).reindex(uni).fillna(0).round().astype("int64")
    ab = lambda p: [f"{p}_a", f"{p}_b"]
    yield pack("GSE65778", "HEK293T", "human", {
        "tunicamycin_1h": dict(stressor="tunicamycin", cls="pleiotropic", kinase="PERK", h=1, dose="1 ug/ml",
                               control=dict(rna=ab("mrna_untr"), ribo=ab("ribo_untr")),
                               stress=dict(rna=ab("mrna_tm"), ribo=ab("ribo_tm"))),
        "tunicamycin_1h_with_ISRIB": dict(stressor="tunicamycin+ISRIB", cls="isr_blocked", kinase="PERK", h=1,
                                          dose="1 ug/ml + 200 nM ISRIB", role="specificity_control",
                                          note="stress effect with eIF2a-P signalling pharmacologically bypassed",
                                          control=dict(rna=ab("mrna_isrib"), ribo=ab("ribo_isrib")),
                                          stress=dict(rna=ab("mrna_tmisrib"), ribo=ab("ribo_tmisrib"))),
    }, c)

    # --- GSE227766: IMR90 cycling +/- thapsigargin 25 nM 3 h (STAR ReadsPerGene)
    d = RAW / "GSE227766/extracted"
    cols = {}
    for f in sorted(d.glob("*ReadsPerGene.out.tab.gz")):
        t = pd.read_csv(f, sep="\t", header=None, index_col=0).iloc[4:]
        cols[re.sub(r"^GSM\d+_|_?ReadsPerGene.*$", "", f.name)] = t[[2, 3]].sum().idxmax() and t[t[[2, 3]].sum().idxmax()]
    c = _to_gene_space(pd.DataFrame(cols), ens, uni)
    yield pack("GSE227766", "IMR90_fibroblast", "human", {
        "thapsigargin_3h": dict(stressor="thapsigargin", cls="pleiotropic", kinase="PERK", h=3, dose="25 nM",
                                note="Tg libraries deposited 18 months after controls: batch confounded with treatment",
                                control=dict(rna=["C1_RNA", "C2_RNA"],
                                             ribo=[f"Clean_T_Cyc_{i}_Ribo_S{i}_R1_001" for i in (1, 2, 3)]),
                                stress=dict(rna=["RNA_CTg1", "RNA_CTg2"], ribo=["Ribo_CTg1", "Ribo_CTg2"])),
    }, c)

    # --- GSE256237: HCT116 (G3BP-AID line, undepleted), arsenite 500 uM 1 h
    t = pd.read_csv(RAW / "GSE256237/GSE256237_cts_riboseq.csv.gz", index_col=1).drop(columns="Unnamed: 0")
    rs = pd.concat([_refseq_to_gene(genes), ens])
    t.index = strip(t.index)
    c = _to_gene_space(t.fillna(0), rs, uni)                 # RefSeq transcripts summed per gene
    r = lambda a, cond, n: [f"{a}_HCT_WT_{cond}_rep{i}" for i in n]
    yield pack("GSE256237", "HCT116", "human", {
        "arsenite_1h": dict(stressor="arsenite", cls="pleiotropic", kinase="HRI", h=1, dose="500 uM",
                            control=dict(rna=r("RNA", "Unstressed", (1, 2)), ribo=r("RPF", "Unstressed", (1, 2))),
                            stress=dict(rna=r("RNA", "oxidative", (1, 2, 3, 4)), ribo=r("RPF", "oxidative", (1, 2, 3, 4)))),
    }, c)

    # --- GSE193945: mouse SCC cells, bortezomib 100 nM 4.5 h; eIF2a-S51A cells as specificity control
    cl = lambda s: re.sub(r"_(aln|Reseq).*$", "", s)
    rp = pd.read_csv(RAW / "GSE193945/GSE193945_raw_counts_ribosome_profiling.csv.gz", index_col=0).rename(columns=cl)
    rn = pd.read_csv(RAW / "GSE193945/GSE193945_raw_counts_totalRNA_reseqMerge.csv.gz", index_col=0).rename(columns=cl)
    t = rp.join(rn, how="outer").fillna(0)
    t.index = t.index.astype(str)
    c = _to_gene_space(t, m2h_ent, uni)
    S = lambda n: [f"Sample_{i}" for i in n if f"Sample_{i}" in c.columns]
    T = lambda n: [f"tot_{i}" for i in n]
    yield pack("GSE193945", "mouse_SCC", "mouse", {
        "bortezomib_4.5h": dict(stressor="bortezomib", cls="pleiotropic", kinase="HRI/GCN2 (unresolved)", h=4.5,
                                dose="100 nM", control=dict(rna=T((1, 2, 3, 4)), ribo=S((1, 2, 3, 4))),
                                stress=dict(rna=T((9, 10, 11, 12)), ribo=S((9, 10, 11, 12)))),
        "bortezomib_4.5h_eIF2a_S51A": dict(stressor="bortezomib(S51A)", cls="isr_blocked", kinase="none", h=4.5,
                                           dose="100 nM", role="specificity_control",
                                           note="eIF2a cannot be phosphorylated",
                                           control=dict(rna=T((5, 6, 7, 8)), ribo=S((5, 6, 7, 8))),
                                           stress=dict(rna=T((13, 14, 15, 16)), ribo=S((13, 14, 15, 16)))),
    }, c)

    # --- GSE282509: HeLa (control line), tunicamycin 1 ug/ml 16 h (late)
    t = pd.read_csv(RAW / "GSE282509/GSE282509_raw_counts_matrix_eIF2A_KO_Roiuk.txt.gz", sep="\t")
    c = _to_gene_space(t.drop(columns=["Transcript stable ID", "GENE"]).groupby("Gene stable ID").sum(), ens, uni)
    r = lambda a, cond, n: [f"{a}_HeLa_control_cells_{cond}_replica_{i}" for i in n]
    yield pack("GSE282509", "HeLa", "human", {
        "tunicamycin_16h": dict(stressor="tunicamycin", cls="pleiotropic", kinase="PERK", h=16, dose="1 ug/ml",
                                note="late time point: includes ATF4-dependent transcription and feedback",
                                control=dict(rna=r("Total_RNA", "dmso", (1, 2, 3)), ribo=r("Riboseq", "dmso", (1, 2, 3))),
                                stress=dict(rna=r("Total_RNA", "TM", (1, 2)), ribo=r("Riboseq", "TM", (1, 2)))),
    }, c)

    # --- GSE291029: RPE-1 TP53-KO, WEE1 inhibitor AZD1775 24 h (GCN2-dependent ISR per the study)
    t = pd.read_csv(RAW / "GSE291029/GSE291029_counts.csv.gz", index_col=0)
    t.index = t.index.str.replace(r"\..*$", "", regex=True)
    md = pd.read_csv(RAW / "GSE291029/GSE291029_metadata.csv.gz")
    c = _to_gene_space(t, ens, uni)
    pick = lambda ty, tr: list(md[(md.type == ty) & (md.treatment == tr)].sort_values("replicate").CB)
    yield pack("GSE291029", "RPE1_TP53KO", "human", {
        "WEE1i_24h": dict(stressor="WEE1_inhibitor", cls="pleiotropic", kinase="GCN2", h=24, dose="AZD1775",
                          note="late; total RNA by VASA-seq; ISR is secondary to replication stress",
                          control=dict(rna=pick("VASA", "DMSO"), ribo=pick("RPF", "DMSO")),
                          stress=dict(rna=pick("VASA", "WEE1i"), ribo=pick("RPF", "WEE1i"))),
    }, c)

    # --- GSE119365: mouse fetal-liver basophilic erythroblasts, iron deficiency; Hri-/- as control
    d = RAW / "GSE119365/extracted"
    cols = {re.sub(r"^GSM\d+_|\.txt\.gz", "", f.name): pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
            for f in sorted(d.glob("*.txt.gz"))}
    c = _to_gene_space(pd.DataFrame(cols).fillna(0), m2h_sym, uni)
    r = lambda a, g, fe, n: [f"{a}_seq_{g}_{fe}_r{i}" for i in n]
    yield pack("GSE119365", "mouse_erythroblast", "mouse", {
        "iron_deficiency": dict(stressor="iron_deficiency", cls="canonical", kinase="HRI", h=np.nan,
                                dose="maternal iron-deficient diet (chronic)",
                                note="chronic in vivo HRI activation; replicates are pooled littermate fetal livers",
                                control=dict(rna=r("rna", "wt", "plus", (1, 2, 3)), ribo=r("ribo", "wt", "plus", (1, 2))),
                                stress=dict(rna=r("rna", "wt", "minus", (1, 2, 3)), ribo=r("ribo", "wt", "minus", (1, 2)))),
        "iron_deficiency_Hri_KO": dict(stressor="iron_deficiency(Hri-KO)", cls="isr_blocked", kinase="none",
                                       h=np.nan, dose="maternal iron-deficient diet (chronic)",
                                       role="specificity_control", note="Hri-/- erythroblasts",
                                       control=dict(rna=r("rna", "ko", "plus", (1, 2, 3)), ribo=r("ribo", "ko", "plus", (1, 2))),
                                       stress=dict(rna=r("rna", "ko", "minus", (1, 2, 3, 4)), ribo=r("ribo", "ko", "minus", (1, 2)))),
    }, c)

    # --- GSE129757: LN308 glioma cells, tunicamycin 2 h and 6 h (thapsigargin Ribo-seq is unreplicated: excluded)
    d = RAW / "GSE129757/extracted"
    cols = {re.sub(r"^GSM\d+_LN308-|\.txt\.gz", "", f.name): pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
            for f in sorted(d.glob("*.txt.gz"))}
    t = pd.DataFrame(cols)
    t.index = strip(t.index)
    c = _to_gene_space(t[~t.index.str.startswith("__")].fillna(0), ens, uni)
    R = (1, 2, 3)
    arm = lambda h, bc: dict(stressor="tunicamycin", cls="pleiotropic", kinase="PERK", h=h, dose="not stated in GEO sample record",
                             control=dict(rna=[f"untreated-R{i}" for i in R], ribo=[f"untreated-GCCAAT-R{i}" for i in R]),
                             stress=dict(rna=[f"{h}h-Tunicamycin-R{i}" for i in R],
                                         ribo=[f"{h}h-Tunicamycin-{bc}-R{i}" for i in R]))
    yield pack("GSE129757", "LN308_glioma", "human", {"tunicamycin_2h": arm(2, "AGGAAT"),
                                                      "tunicamycin_6h": arm(6, "GAATGA")}, c)

    # --- GSE254349 (RNA) + GSE254631 (Ribo): U-2OS G3BP1/2-KO + G3BP1-WT rescue, arsenite 200 uM 2 h
    cols = {}
    for acc in ("GSE254349", "GSE254631"):
        for f in sorted((RAW / acc / "extracted").glob("*.csv.gz")):
            cols[re.sub(r"^GSM\d+_|\.csv\.gz", "", f.name)] = pd.read_csv(f, index_col=1)["Read_Count"]
    t = pd.DataFrame(cols)
    t.index = strip(t.index)
    c = _to_gene_space(t.groupby(level=0).sum(), ens, uni)
    r = lambda a, tr: [f"{a}_G3BP1_{tr}_Rep{i}" for i in (1, 2)]
    yield pack("GSE254631", "U2OS", "human", {
        "arsenite_2h": dict(stressor="arsenite", cls="pleiotropic", kinase="HRI", h=2, dose="200 uM",
                            note="G3BP1/2-KO cells re-expressing G3BP1-WT; RNA series GSE254349",
                            control=dict(rna=r("RNAseq", "water"), ribo=r("Riboseq", "water")),
                            stress=dict(rna=r("RNAseq", "NaAs"), ribo=r("Riboseq", "NaAs"))),
    }, c)


# ------------------------------------------------------------------ TET2 studies
def _tet2_studies(genes: pd.DataFrame):
    """Yield (study_id, samples, counts in source IDs, ortholog map or None)."""
    m2h, info = mouse_to_human(genes, "symbol")
    omap = pd.DataFrame(dict(source_gene_id=m2h.index, gene_id=m2h.values))
    row = lambda sid, study, sp, g, allele, comp, donor: dict(sample_id=sid, study_id=study, species=sp, genotype=g,
                                                              allele_type=allele, compartment=comp, donor=donor)
    # GSE232384: mouse SMART-seq of 200 sorted cells per population; Tet2-KO vs WT animals
    t = pd.read_csv(RAW / "GSE232384/GSE232384_TET2_project_counts.csv.gz", index_col=0)
    t = t.groupby("gene").sum()
    smp = []
    for c in t.columns:
        m = re.fullmatch(r"(T2K|WT)(\d)(CLP|CMP|GMP|LT|MEP|ST)", c)
        smp.append(row(c, "GSE232384", "mouse", "TET2" if m[1] == "T2K" else "WT", "germline/conditional KO (see paper)",
                       {"LT": "LT-HSC", "ST": "ST-HSC"}.get(m[3], m[3]), f"{m[1]}{m[2]}"))
    yield "GSE232384", pd.DataFrame(smp), t.round().astype("int64"), omap

    # GSE212534: mouse 10x scRNA-seq of sorted donor HSC and MPP3, 13 months post-transplant.
    # Each library is one sorted population: summed to one pseudobulk per library.
    import h5py
    from scipy import sparse
    cols, smp = {}, []
    for f in sorted((RAW / "GSE212534/extracted").glob("*.h5")):
        code = re.search(r"_([CK])([HM])(\d)_", f.name)
        with h5py.File(f) as h:
            m = h["matrix"]
            x = sparse.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:])
            s = pd.Series(np.asarray(x.sum(axis=1)).ravel(), index=[n.decode() for n in m["features/name"][:]])
        sid = f"{code[1]}{code[2]}{code[3]}"
        cols[sid] = s.groupby(level=0).sum()
        smp.append(row(sid, "GSE212534", "mouse", "TET2" if code[1] == "K" else "WT", "KO, transplanted donor cells",
                       "HSC" if code[2] == "H" else "MPP3", f"{code[1]}{code[3]}"))
        smp[-1]["n_cells"] = int(x.shape[1])
    yield "GSE212534", pd.DataFrame(smp), pd.DataFrame(cols).fillna(0).astype("int64"), omap

    # GSE183675: mouse LSK SMART-seq, young and old; only library-size-normalised values are deposited.
    t = pd.read_csv(RAW / "GSE183675/GSE183675_WTy-WTo-T2KOy-T2KOo.normalize.txt.gz", sep=" ", index_col=0)
    t = t.groupby(level=0).sum()
    smp = [row(c, "GSE183675", "mouse", "TET2" if c.startswith("T2KO") else "WT", "KO",
               "LSK_young" if c[-2] == "y" else "LSK_old", c) for c in t.columns]
    yield "GSE183675", pd.DataFrame(smp), t.round().astype("int64"), omap

    # GSE287524: human CD34+ HSPCs (and CD14+ monocytes) from long-term culture, CRISPR TET2-mutant vs WT,
    # +/- LPS 6 h. Only untreated libraries enter the basal contrast.
    soft = open("results/audit/GSE287524.soft.txt", errors="replace").read()
    cols, smp = {}, []
    for blk in soft.split("^SAMPLE = ")[1:]:
        title = re.search(r"^!Sample_title = (.*)$", blk, flags=re.M)[1].strip()
        fn = re.search(r"suppl/(GSM\S+genes\.results\.gz)", blk)[1]
        m = re.fullmatch(r"(CD34|CD14)_([WT])(\d)_?(LPS)?", title)
        t = pd.read_csv(RAW / "GSE287524/extracted" / fn, sep="\t", index_col=0)["expected_count"]
        sid = f"{m[1]}_{m[2]}{m[3]}_{'LPS' if m[4] else 'untreated'}"
        cols[sid] = t
        comp = ("CD34_HSPC" if m[1] == "CD34" else "CD14_monocyte") + ("_LPS" if m[4] else "")
        smp.append(row(sid, "GSE287524", "human", "TET2" if m[2] == "T" else "WT", "CRISPR-edited (mutant pool)",
                       comp, f"expt{m[3]}"))
    yield "GSE287524", pd.DataFrame(smp), pd.DataFrame(cols).round().astype("int64"), None


# ------------------------------------------------------------------ driver
def build(out="data/processed/public") -> None:
    out = Path(out)
    (out / "isr").mkdir(parents=True, exist_ok=True)
    genes = build_genes()
    genes.to_parquet(out / "genes.parquet", index=False)
    programs = build_programs(genes)
    programs.to_csv(out / "programs.tsv", sep="\t", index=False)
    ctxs, smps, cnts = [], [], []
    for c, s, n in _isr_studies(genes):
        ctxs += c
        smps += s
        cnts.append(n)
    ctx, smp = pd.DataFrame(ctxs), pd.DataFrame(smps)
    counts = pd.concat(cnts, axis=1)
    counts.index.name = "gene_id"
    counts.to_parquet(out / "isr/counts_all.parquet")
    # depth rule: every arm (condition x assay) needs >=2 libraries with >= MIN_DEPTH assigned reads
    depth = counts.sum()
    smp["assigned_reads"] = depth.reindex(smp.sample_id).values
    ok = (smp.assign(ok=smp.assigned_reads >= MIN_DEPTH).groupby(["context_id", "condition", "assay"]).ok.sum()
          .groupby("context_id").min() >= 2)
    low = ctx.context_id.map(ok).eq(False) & (ctx.role == "train")
    ctx.loc[low, "role"] = "excluded_low_depth"
    ctx.loc[low, "note"] = ctx.loc[low, "note"] + f" | excluded: an arm has <2 libraries with >={MIN_DEPTH:.0e} assigned reads"
    smp = smp[(smp.assigned_reads >= MIN_DEPTH) | ~smp.context_id.isin(ctx.context_id[ctx.role == "train"])]
    ctx.to_csv(out / "isr/contexts_all.tsv", sep="\t", index=False)
    smp.to_csv(out / "isr/samples_all.tsv", sep="\t", index=False)
    train = ctx[ctx.role == "train"]
    train.to_csv(out / "isr/contexts.tsv", sep="\t", index=False)
    s2 = smp[smp.context_id.isin(train.context_id)]
    s2.to_csv(out / "isr/samples.tsv", sep="\t", index=False)
    counts[s2.sample_id].to_parquet(out / "isr/counts.parquet")
    info = {}
    for sid, s, n, omap in _tet2_studies(genes):
        d = out / "tet2" / sid
        d.mkdir(parents=True, exist_ok=True)
        if omap is None:                                      # human Ensembl IDs
            n = n[n.index.isin(genes.gene_id)]
        else:
            omap.to_csv(d / "orthologs.tsv", sep="\t", index=False)
        n.index.name = "gene_id"
        s.to_csv(d / "samples.tsv", sep="\t", index=False)
        n[list(s.sample_id)].to_parquet(d / "counts.parquet")
        info[sid] = dict(n_samples=len(s), n_source_genes=len(n))
    _, oinfo = mouse_to_human(genes, "symbol")
    info["mouse_to_human"] = oinfo
    info["n_genes"] = len(genes)
    info["isr_contexts_by_role"] = ctx.role.value_counts().to_dict()
    (out / "build_info.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))
