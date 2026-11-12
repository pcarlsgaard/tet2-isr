from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import pipeline


def main(argv=None):
    ap = argparse.ArgumentParser(prog="tet2isr")
    ap.add_argument("stage", choices=["smoke-data", "audit-data", "download-processed", "build-public", "analyze-tet2",
                                      "analyze-isr", "fit-baselines", "evaluate", "intersect", "report", "all"])
    ap.add_argument("--config", default="configs/smoke.yml")
    a = ap.parse_args(argv)
    if a.stage == "audit-data":
        from .data_ingest.geo import audit
        audit()
        return
    if a.stage == "download-processed":
        from .data_ingest.geo import download_processed
        download_processed()
        return
    if a.stage == "build-public":
        from .data_ingest.public import build
        build()
        return
    cfg = pipeline.load_cfg(a.config)
    if a.stage != "smoke-data" and not Path(cfg["data_root"], "isr/counts.parquet").exists():
        if cfg.get("synthetic"):
            pipeline.stage_smoke_data(cfg)
        else:
            sys.exit(f"BLOCKED: {cfg['data_root']} is empty. Real-data stages are contingent on "
                     "`make download-processed build-public`.")
    stages = {"smoke-data": [pipeline.stage_smoke_data], "analyze-tet2": [pipeline.stage_tet2],
              "analyze-isr": [pipeline.stage_isr], "fit-baselines": [pipeline.stage_baselines],
              "evaluate": [pipeline.stage_evaluate], "intersect": [pipeline.stage_intersect]}
    stages["all"] = [pipeline.stage_smoke_data] * bool(cfg.get("synthetic")) + [
        pipeline.stage_isr, pipeline.stage_baselines, pipeline.stage_evaluate, pipeline.stage_intersect]
    for f in stages.get(a.stage, []):
        print(f"[tet2isr] {f.__name__}", flush=True)
        f(cfg)
    if a.stage in ("report", "all"):
        from .reporting.report import build
        p = build(cfg)
        pipeline.write_provenance(cfg)
        print(f"[tet2isr] report -> {p}")


if __name__ == "__main__":
    main()
