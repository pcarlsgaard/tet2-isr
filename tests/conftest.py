import pandas as pd
import pytest

from tet2isr import contracts
from tet2isr.data_ingest.synthetic import generate
from tet2isr.isr_translation import dte
from tet2isr.modeling import models

STATE = ["tc_availability", "ribosome_biogenesis", "proliferation"]


@pytest.fixture(scope="session")
def world(tmp_path_factory):
    root = tmp_path_factory.mktemp("syn")
    generate(root, seed=3, n_genes=600)
    ctx, smp, cnt = contracts.load_isr(root)
    sig = dte.build_signatures(ctx, smp, cnt)
    genes = pd.read_parquet(root / "genes.parquet")
    programs = pd.read_csv(root / "programs.tsv", sep="\t")
    c = ctx.set_index("context_id").join(models.context_states(cnt, smp, programs, STATE))
    return dict(root=root, ctx=c, smp=smp, cnt=cnt, sig=sig, genes=genes, programs=programs)
