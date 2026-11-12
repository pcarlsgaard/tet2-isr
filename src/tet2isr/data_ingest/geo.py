"""Manifest audit and processed-file download (contingent on network access to NCBI)."""
from __future__ import annotations

import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REQUIRED = ["accession", "url", "species", "assay", "license", "limitations", "supports_baseline_tet2",
            "supports_isr_transcription", "supports_isr_translation", "supports_direct_genotype_x_isr",
            "relevance_rank", "confounding_risk", "verification_status"]


def _get(url: str, timeout=30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "tet2isr-audit/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def audit(manifest="data/dataset_manifest.tsv", out="results/audit") -> pd.DataFrame:
    """Validate manifest structure and try to fetch each GEO series' brief text record.
    Writes the raw records (when reachable) so sample-level metadata can be reviewed and the
    'unverified' fields filled in by hand; never edits the manifest automatically."""
    m = pd.read_csv(manifest, sep="\t")
    missing = [c for c in REQUIRED if c not in m.columns]
    if missing:
        raise ValueError(f"manifest lacks columns {missing}")
    if (m.supports_direct_genotype_x_isr == "yes").any():
        print("NOTE: a source is marked as a direct genotype x ISR translatome; it must be reserved for validation.")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for acc in m.accession.dropna():
        if not str(acc).startswith("GSE"):
            continue
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=brief"
        try:
            txt = _get(url)
            (out / f"{acc}.soft.txt").write_bytes(txt)
            t = txt.decode(errors="replace")
            rows.append(dict(accession=acc, reachable=True, n_samples=t.count("^SAMPLE"),
                             sha256=hashlib.sha256(txt).hexdigest(), error=""))
        except Exception as e:  # noqa: BLE001 - any network failure is recorded, not hidden
            rows.append(dict(accession=acc, reachable=False, n_samples=None, sha256="", error=str(e)[:120]))
    res = pd.DataFrame(rows)
    res["checked_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    res.to_csv(out / "geo_reachability.tsv", sep="\t", index=False)
    n_unv = int(m.verification_status.str.contains("NOT|not reachable|unverified", case=False).sum())
    print(res.to_string(index=False))
    print(f"\n{len(m)} sources; {n_unv} with unverified fields; {int(res.reachable.sum())}/{len(res)} GEO records reachable.")
    return res


def download_processed(manifest="data/dataset_manifest.tsv", raw="data/raw") -> None:
    """Download GEO series supplementary file listings + files < max size. Never downloads FASTQ."""
    m = pd.read_csv(manifest, sep="\t")
    log = []
    for acc in m.accession.dropna():
        if not str(acc).startswith("GSE"):
            continue
        d = Path(raw) / acc
        d.mkdir(parents=True, exist_ok=True)
        base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/"
        try:
            listing = _get(base + "filelist.txt").decode()
            (d / "filelist.txt").write_text(listing)
            for line in listing.splitlines()[1:]:
                parts = line.split("\t")
                if len(parts) < 4 or parts[1].endswith(".tar") and int(parts[3]) > 2e9:
                    continue
                name = parts[1]
                blob = _get(base + name, timeout=600)
                (d / name).write_bytes(blob)
                log.append(dict(accession=acc, file=name, bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(),
                                url=base + name, downloaded_utc=datetime.now(timezone.utc).isoformat(timespec="seconds")))
        except Exception as e:  # noqa: BLE001
            log.append(dict(accession=acc, file=None, error=str(e)[:160]))
    Path(raw, "download_log.json").write_text(json.dumps(log, indent=1))
    ok = sum(1 for x in log if x.get("file"))
    print(f"downloaded {ok} files; {sum(1 for x in log if x.get('error'))} accessions failed. See {raw}/download_log.json")
    if ok == 0:
        raise SystemExit("BLOCKED: no processed files could be downloaded (NCBI hosts unreachable). "
                         "Allow www.ncbi.nlm.nih.gov and ftp.ncbi.nlm.nih.gov, then re-run.")
