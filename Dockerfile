FROM python:3.13-slim
RUN apt-get update && apt-get install -y --no-install-recommends make git && rm -rf /var/lib/apt/lists/*
WORKDIR /work
COPY . .
RUN python -m venv .venv && .venv/bin/pip install --no-cache-dir -r requirements.lock && .venv/bin/pip install --no-deps -e .
CMD ["make", "smoke-test"]
# Apptainer: apptainer build tet2isr.sif docker-daemon://tet2isr:latest
