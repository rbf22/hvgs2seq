# hvgs2seq
Simple package that will take as input a list of hgvs variants and output the final sequence. It will try to be cleaver and deal with insertions, deletions, missense, and other variants. It will hopefully deal with NMD and splicing variants as well.

Great—scoping to **one gene / one transcript** lets us make `hgvs2seq` much leaner and more deterministic. Below is a trimmed architecture, data model, and an agent-ready task plan tailored to a **single, fixed transcript** (e.g., a chosen RefSeq `NM_...` or Ensembl `ENST...`). I’ll also show concrete API/CLI shapes and edge-case handling so engineering can ship a focused v0.1 quickly, then layer NMD and splicing in v0.2.

---

# Scope & assumptions (single transcript)

- A **single canonical transcript** is provided up-front (e.g., `NM_XXXX.Y`), along with the **reference assembly** (e.g., GRCh38).  
- All variants refer to this transcript or can be projected to it (we’ll enforce that).  
- We will return edited **cDNA**, **mRNA**, and **protein** sequences for this transcript only.  
- **Phasing** is optional. If supplied, we return haplotype-specific sequences (max 2); if not, we apply variants in a defined order and warn on collisions.  
- **Splicing** (v0.2+): predictions are annotations against this one transcript; optional plugin (e.g., SpliceAI) to suggest exon changes; by default we **don’t** alter exon structure automatically.  
- **NMD** (v0.2+): rule-based call on the edited transcript for this one isoform.

We’ll still reuse proven bricks (biocommons/hgvs + SeqRepo/UTA) for parsing/projection/data, but our logic and I/O are streamlined because there’s no multi-transcript arbitration.

---

# Minimal architecture (single transcript)

```
hgvs2seq/
  __init__.py
  config.py                 # transcript + assembly config
  parse.py                  # parse/validate HGVS strings
  project.py                # project variants -> chosen transcript coordinates
  refseq.py                 # fetch reference sequences (SeqRepo/UTA)
  apply/
    single.py               # apply one normalized variant -> string
    batch.py                # merge/order/apply many variants (haplotypes optional)
  consequence/
    cds.py                  # rebuild CDS, translate to protein, detect consequence
    nmd.py                  # (v0.2) NMD rules for this transcript
  splicing/
    spliceai.py             # (v0.2) optional splice prediction wrapper
  io/
    fasta.py, jsonio.py     # outputs with provenance
  cli.py
```

### Data model (typed pydantic-style)

```python
class TranscriptConfig(BaseModel):
    transcript_id: str        # NM_/ENST
    gene_symbol: str | None
    assembly: str             # e.g. GRCh38
    strand: int               # +1/-1
    exons: list[tuple[int,int]]  # genomic coordinates (1-based, inclusive)
    cds_start_c: int | None   # cDNA index of CDS start (ATG), e.g., c.1
    cds_end_c: int | None     # cDNA index of last CDS base

class VariantIn(BaseModel):
    hgvs: str                 # original string
    phase_group: int | None   # same int => phased together (optional)

class VariantNorm(BaseModel):
    hgvs_c: str               # normalized HGVS on transcript (c.)
    kind: Literal["sub","del","ins","delins","dup","inv","fs","mnv"]
    c_start: int
    c_end: int
    alt: str                  # inserted/alt bases normalized
    meta: dict                # original string, notes

class EditPlan(BaseModel):
    haplotypes: list[list[VariantNorm]] # one list per haplotype
    policy: Literal["order_by_pos","reject_overlaps","left_shift"]
    warnings: list[str]

class SequenceBundle(BaseModel):
    cdna_ref: str
    cdna_edited: list[str]        # one per haplotype
    mrna_ref: str                 # same as cdna_ref if no UTR intron nuance
    mrna_edited: list[str]
    protein_ref: str | None
    protein_edited: list[str] | None
    annotations: dict             # consequence, frameshift/stop, splice/NMD, etc.
    provenance: dict              # versions, transcript, assembly, variants applied
```

---

# Single-transcript pipeline

1) **Config load**  
   - Load `TranscriptConfig` (JSON/YAML), or derive from UTA once and cache. Store exon structure, CDS bounds, strand.

2) **Parse & validate** (simple!)  
   - Accept list of HGVS strings. Parse/normalize with `biocommons/hgvs`. Hard-reject if not projectable to the chosen transcript.  
   - Project all to **c. coordinates** for that transcript (we operate in c. space for edits).

3) **Plan & apply edits**
   - Group variants by `phase_group` → haplotypes, else all in one.  
   - Enforce **non-overlap** within a haplotype (or resolve by position policy).  
   - Apply variants to **cDNA reference** string: implement precise index math for `sub/ins/del/delins/dup/mnv`.  
   - Produce `cdna_edited` per haplotype.

4) **Consequence & protein**
   - Rebuild **CDS** from edited cDNA (use `cds_start_c`..`cds_end_c`), translate to protein, and annotate effect: missense, synonymous, frameshift, stop-gain/loss, start-loss.  
   - Emit codon/AA diffs with positions (e.g., p.Gly12Val, p.Trp26Ter, etc.).

5) **(v0.2) Splicing annotations**  
   - For variants within ±50 bp of exon boundaries (configurable), run SpliceAI (or load precomputed scores).  
   - Annotate Δ acceptor/donor scores; **do not** change exon structure automatically. Provide an API hook allowing a caller to tell us “skip exon 7” or “extend exon 5 by N bases” and then we recompute cDNA/mRNA/protein.

6) **(v0.2) NMD rules** (single-transcript)  
   - If a PTC appears **>50–55 nt upstream** of the final exon–exon junction → **NMD-likely**, except: last exon PTC, intronless transcript, very short last exon, near-start exceptions. Produce a structured rationale (junction index, distances).

7) **Outputs**  
   - FASTA: `>transcript|haplotype=1`, edited cDNA (and protein).  
   - JSON: `SequenceBundle` with annotations + provenance (reference versions, tool versions, policies).

---

# API & CLI (single transcript)

### Python

```python
from hgvs2seq import apply_variants, load_config

cfg = load_config("NM_000000.0.json")  # TranscriptConfig
variants = [
    {"hgvs": "NM_000000.0:c.76A>T", "phase_group": 1},
    {"hgvs": "NM_000000.0:c.123_124del", "phase_group": 1},
    {"hgvs": "NM_000000.0:c.301+1G>A"}  # unphased, splicing-site
]

bundle = apply_variants(cfg, variants, output={"cdna": True, "protein": True})
print(bundle.provenance)
print(bundle.annotations["haplotype_1"]["protein_consequence"])
```

### CLI

```bash
hgvs2seq \
  --transcript NM_000000.0 \
  --assembly GRCh38 \
  --variants variants.txt \
  --out edited.json --fasta edited.fasta \
  --policy order_by_pos \
  --emit protein,cdna
```

Where `variants.txt` contains one HGVS per line (optional `#phase=1` comment).

---

# Edge cases we now simplify (thanks to single transcript)

- **Transcript selection**: gone—we already know the isoform.  
- **Coordinate ambiguity**: every input must project to this transcript; otherwise we reject with a clear error.  
- **UTR edits**: we can still apply them; protein changes only if they affect start/stop or frameshifts into CDS.  
- **Overlaps**: within a haplotype, either reject or deterministically order (config).  
- **Intronic/splice**: annotate in v0.1; simulate exon changes only if user instructs (v0.2 with plugin).  
- **Phasing**: if provided, return one edited sequence per haplotype; if not, just one.

---

# Testing strategy (tight, single transcript)

- **Golden cases on a fixed toy transcript** (short exons across 3 exons on + and – strands):  
  - SNV → missense/synonymous.  
  - In-frame del/ins/delins; frame-shift del/ins causing early PTC.  
  - `dup` and repeat-region edge behavior.  
  - UTR mutations: Kozak/start-loss; stop-loss into 3′ UTR.  
  - Splice-site SNV at +1/–1 (annotated only).  
  - NMD rule examples: PTC in penultimate exon vs last exon.  
- **Property tests** (Hypothesis): round-trip apply vs reverse-edit; index math near edges.  
- **Performance**: batch of 10–100 variants on one transcript; ensure near-linear time.

---

# Agent-ready task plan (single transcript)

## Phase 0 — Bootstrap
1. **Create repo scaffold (single-isoform)**  
   - Output: package layout above, `pyproject.toml`, CI, `ruff + mypy` config, basic `README`.

2. **Transcript config loader**  
   - Implement `load_config()` from JSON/YAML.  
   - If `--auto`, fetch from UTA once (exons, CDS bounds, strand), cache to JSON.  
   - Unit tests: serialize/deserialize fidelity.

## Phase 1 — Parse & project (single transcript)
3. **HGVS parse/normalize for one transcript**  
   - Parse strings → `VariantNorm` in **c.** coordinates for the provided transcript.  
   - Reject anything that cannot project to that transcript.  
   - Tests with valid/invalid HGVS forms.

## Phase 2 — Apply engine (cDNA only)
4. **Single-variant applier (cDNA)**  
   - Implement exact index math for `sub/ins/del/delins/dup/mnv`.  
   - Include robust bounds checking & helpful errors.

5. **Batch applier + phasing**  
   - Group by `phase_group`; enforce non-overlap or apply `order_by_pos`.  
   - Return `cdna_edited` per haplotype + warnings.

## Phase 3 — Consequence computation
6. **CDS extraction & translation**  
   - From edited cDNA, extract CDS using `cds_start_c..cds_end_c`, translate with standard genetic code.  
   - Emit consequence classification (syn/mis/fs/stopgain/stoploss/startloss) and amino-acid diffs.

7. **Provenance & outputs**  
   - FASTA emitters for cDNA/protein; JSON `SequenceBundle` with full metadata and decisions.

## Phase 4 — (Optional v0.2) Splicing & NMD
8. **Splice annotation plugin (one transcript)**  
   - SpliceAI wrapper (optional dependency).  
   - For variants within configurable windows of exon boundaries, compute/report Δscores and nearest splice sites.  
   - Provide **manual exon override API** (e.g., `apply_splicing(decision={"skip_exons":[7]})`).

9. **NMD rule engine (one transcript)**  
   - Implement standard NMD checks on edited transcript: PTC distance to last exon–exon junction, last exon exception, intronless exception, etc.  
   - Emit structured rationale (junction index, nt distances).

## Phase 5 — Docs & examples
10. **Examples**  
    - Notebook: apply a list of HGVS to `NM_…`, view cDNA/protein changes, show an NMD-likely example, and a splice-site variant with Δscores.  
    - CLI walkthrough.

---

# Suggested interfaces (thin and opinionated)

```python
def apply_variants(
    cfg: TranscriptConfig,
    variants: list[dict | VariantIn],
    *,
    policy: Literal["order_by_pos","reject_overlaps"] = "order_by_pos",
    outputs: set[str] = {"cdna","protein"},
    annotate_splicing: bool = False,
    annotate_nmd: bool = False,
    spliceai_params: dict | None = None
) -> SequenceBundle:
    ...
```

Return values will always include: edited sequences (per haplotype), consequence summaries per haplotype, and a `provenance` block with: transcript ID, assembly, exon count, cds bounds, tool/library versions, policy, and a replayable list of normalized variants (with their original HGVS).

---

# What we’re **not** doing in v0.1 (explicitly)

- No multi-transcript arbitration, no exon graph exploration, no automatic splice-altered transcript reconstruction.  
- No VCF/SPDI/VRS I/O (can come later if needed).  
- No large-SV handling beyond HGVS that maps cleanly in c. space.

---

If you’d like, I can generate:
- a **starter repo** (package skeleton + `TranscriptConfig` schema + CLI that reads a transcript JSON and a `variants.txt` and emits edited cDNA/protein); and
- a **toy transcript fixture** with 3 exons and unit tests you can run immediately.

Say the word and I’ll drop the scaffold here.
