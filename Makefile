# CONFIG=configs/smoke.yml (synthetic, default) or configs/public.yml (needs downloaded public data)
CONFIG ?= configs/smoke.yml
PY := .venv/bin/python
RUN := PYTHONPATH=src $(PY) -m tet2isr.cli

.PHONY: intersect public build-public setup smoke-test test audit-data download-processed analyze-tet2 analyze-isr fit-baselines evaluate report clean
setup:
	python3 -m venv .venv && .venv/bin/pip install -q -r requirements.lock && .venv/bin/pip install -q -e . --no-deps
test:
	PYTHONPATH=src $(PY) -m pytest -q
smoke-test: test          ## full synthetic demo: data -> M0..M3 -> validation -> report (no downloads)
	$(RUN) all --config configs/smoke.yml
audit-data:               ## validates manifest; saves GEO sample records to results/audit/
	$(RUN) audit-data
download-processed:       ## files in configs/public_sources.tsv; never fetches FASTQ
	$(RUN) download-processed
build-public:             ## GEO processed files -> data contract (data/processed/public)
	$(RUN) build-public
public: download-processed build-public   ## full real-data run
	$(RUN) all --config configs/public.yml
analyze-isr:
	$(RUN) analyze-isr --config $(CONFIG)
analyze-tet2:
	$(RUN) analyze-tet2 --config $(CONFIG)
fit-baselines:            ## M0 + M1 (+ 2x2 example); needs analyze-isr
	$(RUN) fit-baselines --config $(CONFIG)
evaluate:                 ## M2/M3, held-out splits, ablations, leakage, genotype predictions
	$(RUN) evaluate --config $(CONFIG)
intersect:                ## TET2 transcriptome programs x stress translatomes
	$(RUN) intersect --config $(CONFIG)
report:
	$(RUN) report --config $(CONFIG)
clean:
	rm -rf results/smoke data/processed/smoke reports/research_report_synthetic_demo.html
