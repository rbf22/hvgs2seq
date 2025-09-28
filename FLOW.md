Love this—let’s lock the **end-to-end flow** so `hgvs2seq` can deterministically go from “list of HGVS on one transcript” → “final mRNA(s)” → “NMD verdict” → “protein(s)”. Below is the flow I recommend, with the decision points for splicing, NMD, frameshifts, new/ lost stops, and a few extras that matter in practice.

---

# Single-transcript outcome flow

## 0) Inputs & assumptions
- **Given:** one canonical transcript (exons, strand, CDS start/stop), reference assembly, list of HGVS variants (optionally phased).
- **Goal order:**  
  1) Determine the **transcript (mRNA) sequence(s)** produced.  
  2) Determine whether the **transcript is NMD-susceptible**.  
  3) Determine the **protein sequence(s)** produced.  
  4) Provide alternatives + confidence when splicing is uncertain.

---

## 1) Parse → project → normalize
1. Parse each HGVS; project to the chosen transcript’s **c.** coordinates.  
2. Normalize (trim/left-shift) and validate bounds.  
3. Group by haplotype (if phasing provided), else single group.  
4. Resolve overlaps inside each group (policy: `reject_overlaps` or `order_by_pos`).

**Outputs:** normalized per-haplotype variant lists on **c.** coordinates + warnings.

---

## 2) Build the pre-edit reference strands
5. Fetch the transcript’s **cDNA** reference (spliced exons).  
6. Hold the **genomic exon map** (for splice distance calculations and NMD exon-junction logic).

---

## 3) Splicing impact assessment (creates transcript scenarios)
We produce **one or more splicing scenarios** per haplotype:

### 3a) Identify splicing-relevant variants
7. Flag variants in:  
   - **Essential splice sites** (±1/±2), branch point windows, polypyrimidine tract.  
   - **Near-splice windows** (configurable, e.g., ±50 nt from exon boundaries).  
   - **Deep intronic** variants (optional window for cryptic exons).  

### 3b) Score predicted splice effects (optional plugin; configurable)
8. For flagged variants, run/lookup splicing predictions (e.g., SpliceAI) and capture Δ acceptor/donor scores.  
9. Generate **candidate splicing changes** above thresholds, e.g.:
   - Exon **skipping** / **retention**  
   - **Donor/acceptor shift** → exon **extension/shortening**  
   - **Cryptic exon** creation

### 3c) Assemble transcript **scenarios**
10. Always include **Baseline scenario**: exon structure unchanged.  
11. For each high-confidence splice event, create an **Altered scenario** with the corresponding exon structure.  
12. If multiple events apply, compose them if compatible; otherwise keep separate branches.  
13. Score each scenario’s **confidence** (e.g., max Δscore, rule class, multiple-evidence bonus).

**Outputs:** for each haplotype, **one or more transcript scenarios** (exon lists with c. coordinate transforms) with confidence scores.  
If splicing plugin is off: only Baseline scenario.

---

## 4) Apply sequence edits (per haplotype × scenario)
14. For each scenario, generate the **pre-mRNA exon model** and its **cDNA** sequence.  
15. Apply the normalized c. variants to that cDNA (indels, delins, dup, SNVs, MNVs, inversions supported).  
16. Produce the **edited cDNA (mRNA)** for the scenario.

**Note:** Indels that land in introns in Baseline but become exonic in an Altered scenario (or vice versa) are naturally handled because we apply to the scenario’s cDNA.

**Outputs:** per haplotype × scenario → **mRNA_edited**.

---

## 5) Determine CDS and translate (frameshifts, new stops, stop loss)
For each mRNA_edited:

### 5a) (Re)locate start and stop
17. Start codon:
   - If reference start codon is intact → use annotated start.  
   - If **start-loss** (mutated ATG or strong Kozak disruption if you choose to model it), search for **downstream in-frame ATG** (configurable Kozak scan). If found, set **ALT start**; else **no CDS** (protein=None).

18. Stop codon:
   - If reference stop still in frame and intact → keep.  
   - If variant causes **stop-loss**, scan in-frame to the next stop; if none before poly(A) tail (we don’t know exact tail), treat as **nonstop** (see decay below).

### 5b) Translate & classify
19. Extract CDS (start..stop or to end if no stop).  
20. Translate:
   - Detect **frameshift** (indel not multiple of 3 within CDS) → protein changes downstream until first stop → **stop-gain** likely (PTC).  
   - Detect **missense/synonymous** for in-frame SNVs/MNVs.  
   - Detect **in-frame indels** (codon gain/loss) with precise AA diff.  
   - Detect **stop-gain** (new PTC introduced by SNV/frameshift).  
   - Detect **stop-loss** (lost stop; extended C-terminus length).  
   - Optionally annotate **uORFs** if 5′UTR creates upstream ATG with good Kozak (see extras).

**Outputs:** per mRNA_edited → **protein_edited** (or None), and consequence summary (frameshift, stop-gain, stop-loss, missense, synonymous, start-loss with reinit, etc.), with HGVS.p where unambiguous.

---

## 6) NMD decision (per mRNA_edited)
We evaluate **NMD** on the translated mRNA (whether or not we output protein—NMD is a property of the transcript with a PTC):

21. If **no PTC** (i.e., normal stop in last exon) → **NMD-unlikely**.  
22. If **PTC present**, compute:
   - Positions of **exon–exon junctions** (EJCs) in the edited transcript scenario.  
   - Distance from PTC to the **last exon–exon junction**.

23. Apply canonical rules (configurable):
   - **NMD-likely** if PTC is **> 50–55 nt upstream** of the last EJC.  
   - **NMD-escape** if PTC in **last exon** or within **≤50 nt** upstream of the last EJC.  
   - **Intronless transcript** → NMD-escape.  
   - Optional refinements (enable/disable): very short last exon exception, long 3′UTR nuances, near-start PTC quirk.

24. Record a **structured rationale**: exon index, nt distances, rule fired.

**Outputs:** per mRNA_edited → **nmd_status** ∈ {likely, escape, not_applicable}, with rationale.

> Policy: We still compute the **theoretical protein** for PTC cases, but mark it as **“not expected in vivo due to NMD-likely”** unless user requests otherwise.

---

## 7) Special decay modes (optional annotations)
25. **Nonstop decay (NSD):** if **stop-loss** and no in-frame stop codon before transcript end → flag **NSD-risk** (no canonical stop).  
26. **No-go decay (NGD):** not directly inferable from sequence; we won’t model.  
27. **uORFs in 5′UTR:** if new upstream ATG with good Kozak and downstream stop before main start → annotate **uORF created** (may reduce translation of main ORF; informational).

---

## 8) Choose “primary” outcome vs alternatives
28. If splicing plugin **off** → only Baseline scenario → that’s the **primary outcome**.  
29. If plugin **on** and multiple scenarios exist:
   - Choose the **highest-confidence** splicing scenario as **Primary**.  
   - Emit others as **Alternates** with confidence scores and clear differences (exon changes, CDS shifts, NMD status, protein diffs).  
   - Allow user override to pin a scenario (e.g., `--assume-exon7-skip`).

---

## 9) Outputs & provenance
30. For each **haplotype × scenario**:
   - **mRNA (FASTA)**, **protein (FASTA)** (or None), and a **JSON** block with:  
     - `transcript_outcome`: exon set, cDNA sequence digest  
     - `nmd`: status + rationale  
     - `protein_outcome`: AA sequence digest, HGVS.p, effect classification  
     - `extras`: stop-loss extension length, frameshift location, uORF notes, NSD-risk  
     - `splicing_evidence`: Δscores, decisions, thresholds  
     - `provenance`: versions, reference build, transcript ID, policies, normalized HGVS list, phasing, warnings  
     - `confidence`: per scenario overall confidence

31. If multiple haplotypes: label outputs `haplotype=1/2`.

---

# Decision tree (compact pseudocode)

```python
for hap in haplotypes:
    scenarios = [baseline_splice_model()]
    if splicing_enabled:
        scenarios += infer_splice_altered_models(hap.variants, thresholds)

    for sc in scenarios:
        mrna = apply_variants_to_cdna(sc, hap.variants)

        cds = locate_cds(mrna, start_policy, stop_policy)
        if not cds:
            protein = None
            ptc = None
        else:
            protein, consequence, ptc = translate_and_classify(mrna, cds)

        nmd_status, rationale = nmd_check(sc.exon_junctions, ptc)

        extras = analyze_special_cases(mrna, cds, ptc)  # stop-loss length, NSD-risk, uORFs

        emit(hap, sc, mrna, protein, consequence, nmd_status, rationale, extras)
```

---

# What this flow covers (your checklist)

- ✅ **What transcript will be created**  
  - Baseline (no splice change) and any **splice-altered** transcripts (skip, extend, cryptic, retention), each with cDNA sequence.

- ✅ **If it will be destroyed by NMD**  
  - Standard 50–55 nt rule, last-exon exception, intronless exception, with structured rationale.

- ✅ **What protein will be created**  
  - Handles **frameshifts** (indels), **new stop codons** (PTC), **stop codon loss** (C-term extension, NSD-risk), **missense/synonymous**, **start-loss with reinit**.

- ➕ **Extras you might want**  
  - uORF creation in 5′UTR, Kozak context changes (affect translation initiation likelihood), CDS re-definition when exons shift, per-scenario **confidence** scoring, and user overrides for splicing decisions.

---

# Practical defaults (tunable)
- **Splice window:** ±50 nt from exon boundaries.  
- **SpliceAI thresholds (if used):** e.g., `Δ ≥ 0.2` = candidate, `Δ ≥ 0.5` = high-confidence.  
- **Start re-init scan:** downstream in-frame ATG within first ~150 nt, Kozak “strong” only (A/G at –3, G at +4).  
- **NMD:** classic rule; report distance in nt and which junction counted.  
- **Protein reporting under NMD-likely:** include theoretical AA seq but mark `expected_absent=True`.

---

# Minimal outputs we should persist for downstream tools

```json
{
  "primary": {
    "haplotype": 1,
    "scenario_id": "baseline",
    "mrna_fasta": "…",
    "protein_fasta": "… or null …",
    "nmd": {"status": "likely|escape|not_applicable", "ptc_c": 432, "distance_to_last_ejc_nt": 118, "rule": "50nt"},
    "consequence": {"type": "frameshift_stop_gain", "first_aa_change": "p.Gly12Valfs*7"},
    "extras": {"stop_loss_extension_nt": 0, "nsd_risk": false, "uorf_created": false},
    "confidence": 0.92
  },
  "alternates": [ … ],
  "provenance": {
    "transcript": "NM_XXXX.Y",
    "assembly": "GRCh38",
    "exon_count": 12,
    "cds_bounds_c": [1, 1458],
    "tools": {"hgvs2seq": "0.1.0", "hgvs": "x.y", "spliceai": "z (if used)"},
    "policy": {"overlap": "order_by_pos", "start_reinit": "strict_kozak"},
    "variants_normalized": [ "NM_XXXX.Y:c.76A>T", "NM_XXXX.Y:c.123_124del", … ],
    "warnings": [ … ]
  }
}
```

---

If this flow looks good, I can turn it into:
- a short **spec doc** for engineers (state machine + data schemas), and
- a **unit-test checklist** that guarantees every branch (frameshift, stop-gain, stop-loss, NMD yes/no, splice skip/extend) is exercised on a toy transcript.
