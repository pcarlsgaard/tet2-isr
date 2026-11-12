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
    import re
    for acc in dict.fromkeys(a for x in m.accession.dropna() for a in re.findall(r"GSE\d+", str(x))):
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
    n_unv = int(m.verification_status.str.contains("not re-read|not downloaded|unverified", case=False).sum())
    print(res.to_string(index=False))
    print(f"\n{len(m)} sources; {n_unv} not fully verified (see verification_status); {int(res.reachable.sum())}/{len(res)} GEO records reachable.")
    return res


def download_processed(sources="configs/public_sources.tsv", raw="data/raw") -> None:
    """Download the processed supplementary files listed in `sources` (never FASTQ), logging
    URL, size, SHA-256 and time. Existing files are kept and re-hashed."""
    src = pd.read_csv(sources, sep="\t")
    log = []
    for acc, name in zip(src.accession, src.file):
        d = Path(raw) / acc
        d.mkdir(parents=True, exist_ok=True)
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{name}"
        f = d / name
        try:
            if not f.exists():
                f.write_bytes(_get(url, timeout=900))
            log.append(dict(accession=acc, file=name, bytes=f.stat().st_size,
                            sha256=hashlib.sha256(f.read_bytes()).hexdigest(), url=url,
                            logged_utc=datetime.now(timezone.utc).isoformat(timespec="seconds")))
        except Exception as e:  # noqa: BLE001
            log.append(dict(accession=acc, file=name, url=url, error=str(e)[:160]))
    Path(raw, "download_log.json").write_text(json.dumps(log, indent=1))
    bad = [x for x in log if x.get("error")]
    print(f"{len(log) - len(bad)} files present; {len(bad)} failed. Log: {raw}/download_log.json")
    for x in bad:
        print("  FAILED", x["accession"], x["file"], x["error"])
    if len(bad) == len(log):
        raise SystemExit("BLOCKED: no processed files could be downloaded (NCBI unreachable).")
