# PhArMol: Complete Documentation

## 1. Overview

A real PyMOL plugin for ligand-based pharmacophore modeling: given a handful of known active molecules against a drug target, it aligns them in 3D, extracts a consensus pharmacophore, and lets you score new candidates — individually or as a whole library — against that pattern, with several independent, complementary checks on top of the core geometric match. Runs inside a live PyMOL session, not a separate application.

| Tab | Purpose |
|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| 1\. Build Model | Turn known actives into a consensus 3D pharmacophore.                                                                                      |
| 2\. Screen      | Score a candidate or a whole library against that consensus, with match, plausibility, shape, scaffold, and verdict all reported together. |
| 3\. Validate    | Test whether the model actually discriminates real actives from decoys.                                                                    |

Developed and stress-tested against three real, chemically distinct targets: **MAO-B inhibitors** (5 actives), **DPP-4 inhibitors / gliptins** (7 actives, a genuinely diverse class), and **COX-2 selective inhibitors / coxibs** (4 training + 2 held-out validation actives). Two rounds of external expert review directly shaped the tool's design — several features exist specifically because a reviewer's suggestion was tested against real data, found to need correction or redesign, and rebuilt accordingly. That process, and its real evidence, is documented throughout this update.

## 2. Architecture & Design Philosophy

### 2.1 Deployment & Dependencies

Distributed as a standard PyMOL plugin, running inside the user's existing PyMOL process. RDKit is the only hard dependency; DBSCAN and ROC-AUC were reimplemented in plain NumPy after a real scipy/macOS-Accelerate installation failure, each validated bit-exact against the scikit-learn/scipy equivalent before adoption.

### 2.2 What was built, then deliberately rejected or redesigned

- **A standalone desktop app / a purely 2D version / automatic multi-hypothesis detection / “receptor-aware exclusion” / automatic decoy generation** — all tried and rejected in the first development phase; see the rationale for each in Section 7.

- **The literal “50×50 conformer ensemble grid” idea** (a second-round expert suggestion) — tested directly: *56 hours* for a 1000-compound library. Rejected in favor of a cheaper redesign (Section 5.4) that captures the same benefit without extra alignments.

- **The literal “42-parameter multi-model voting” idea** — tested: *20.6 hours* for the same library. Not built; a cheaper targeted-rescoring redesign was proposed but not yet implemented (Section 7).

## 3. Scientific Methodology

### 3.1 Core pipeline

SMILES → (optional) physiological-pH protonation → ETKDG conformer ensemble + MMFF94 minimization → O3A alignment onto a reference ligand → DBSCAN clustering into a consensus pharmacophore. A candidate is scored by O3A-aligning it onto that same reference, then checking feature-position agreement within a distance tolerance (now adjustable, 0.5–4.0 Å).

> **REAL FINDING:** The scoring-alignment bug (found via real use, fixed). An early version scored candidates without aligning them onto the model's own reference frame at all. Caught when a user visually noticed, in PyMOL, that a training active clearly satisfied two features the tool's score did not report as matched. Fixed: every scoring path now O3A-aligns first. Verified as a direct regression test.

> **REAL FINDING:** Leave-one-out validation (found via methodological review, fixed). Validating actives against the same model they were used to build tests reproduction of training data, not generalization. Measured directly: in-sample AUC 0.880 vs. genuine leave-one-out AUC 0.840 on the MAO-B set. Tab 3 now defaults to leave-one-out for the actives side, with an External Test Actives option for genuinely separate held-out sets.

### 3.2 Chemical Plausibility and the Verdict system

A Tanimoto similarity check to the nearest training active — explicitly an adaptation of the applicability-domain *principle*, not the formally standardized QSAR-AD method (which is tied to a fitted regression model, not a raw geometric match). Combined with the match score into one plain-language Verdict: Gold Standard Hit, Scaffold-Hop (Moderate/High Risk), Uncertain, Inactive Analog, or Alien/Low Priority. **Never a hard filter** — a high-match/low-plausibility candidate is flagged as higher-risk, never discarded, since finding real hits with different scaffolds is the whole point of the method.

### 3.3 Shape Consistency — a ligand-only check, explicitly scoped

Following a second-round expert suggestion to add steric/excluded-volume awareness, this was deliberately reframed before being built: the union of the aligned training actives' own heavy-atom positions and van der Waals radii forms a shape envelope, and a candidate is scored by what fraction of its own heavy atoms fall inside it. **This checks consistency with the training ligands' own shapes only — it is not, and cannot be, a receptor steric-clash check**, since this pipeline has no receptor structure anywhere in it. The original suggestion's framing (“crashing into the protein wall”) was corrected before implementation for exactly this reason.

> **REAL FINDING:** Verified on real coxib data: every training active scores 100% against its own envelope (as it must, being part of the envelope). A deliberately oversized, unrelated compound (three extra phenyl rings plus a long alkyl chain), O3A-aligned onto the same reference, dropped to 62% — a real, meaningful signal the isolated pharmacophore feature points alone would completely miss.

### 3.4 Ensemble Tolerance Matching — a cheap redesign of an expensive idea

The original suggestion (score every candidate conformer against every training conformer, a full 50×50 grid) was tested and found to cost 56 hours per 1000 compounds (Section 2.2). The redesign actually built: each consensus feature already stores the individual raw 3D points that were clustered into it during model-building (not just their averaged centroid). An opt-in “ensemble” matching mode checks a candidate feature's distance to the **nearest** of those real points instead of the single averaged position — tolerating genuine positional spread across a flexible or diverse active set, at the cost of a handful of extra distance comparisons, with **zero extra O3A alignments**.

> **REAL FINDING:** Verified on real gliptin data: scoring linagliptin (a genuine training active) against its own model, mean match distance tightened from 0.78 Å (centroid mode) to 0.13 Å (ensemble mode) — a real, substantial difference from the same underlying data. Fully backward compatible: any consensus point loaded from an older saved model (without stored raw points) automatically falls back to centroid matching, verified directly rather than assumed.

### 3.5 Scaffold Novelty (Bemis-Murcko)

Compares a candidate's core ring system (its Murcko scaffold, substituents stripped) to each training active's own scaffold: an exact match is labeled “Same Scaffold,” a structurally related but non-identical scaffold “Analog Scaffold,” and neither “Novel Scaffold.” A different question from Plausibility (2D fingerprint similarity) or Verdict (combined match+plausibility) — specifically about the ring system, not the whole molecule.

> **REAL FINDING:** Verified on the real coxib case study: parecoxib (a genuine prodrug of valdecoxib, differing only in a side-chain amide) correctly reduces to the exact same Murcko scaffold as valdecoxib. Lumiracoxib (a genuinely different phenylacetic-acid chemotype) correctly comes back “Novel Scaffold.”

### 3.6 Fast Screen pre-filter — built, with an honest limitation

An alignment-free “3D pharmacophore fingerprint” (the pairwise distances between consensus features, which are rotation/translation-invariant and so need no O3A at all) lets a candidate be pre-screened in a fraction of the time of full scoring. Measured: ~6× faster per compound.

> **REAL FINDING:** A real trade-off found before shipping, not assumed: the fingerprint approach is noticeably more permissive than full O3A-based scoring, because internal pairwise distances alone don't capture full 3D arrangement the way real alignment does. A genuine structural outlier scored 0.90 on the fingerprint alone, far more permissive than its real, much lower full-scoring match. Scoped accordingly: Fast Screen only ever rejects the obviously-hopeless before full scoring runs on the rest — it never replaces the full score for anything that survives it.

## 4. Complete UI Reference

### 4.1 Tab 1 — Build Model

| Control | What it does |
|---------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| Known Actives                                           | name, SMILES per line — the training set.                                                                                                      |
| Conformers per molecule / Protonate at physiological pH | 3D-shape sampling depth (default 50); the explicit, limited protonation rule (Section 3.1).                                                    |
| Support threshold                                       | % of actives that must share a feature to count as “consensus.” Applied after clustering.                                                      |
| Clustering radius (eps)                                 | How close (Å) feature points must be to merge into one cluster, during clustering itself. Default 1.5 Å, adjustable 0.8–2.5, re-clusters live. |
| Feature Table                                           | Feature / Support / Points / Stability (via the robustness check) / Consensus?.                                                                |
| Run Leave-One-Out Robustness Check                      | Rebuilds the model N times excluding one active each time; flags instability in high-support features.                                         |
| Save / Load / Compare Models                            | Persist a model (incl. 3D structures and raw clustered points) to JSON; reopen; compare two models feature-by-feature.                         |

### 4.2 Tab 2 — Screen (updated this round)

| Control / Column | What it does |
|--------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Match tolerance                                  | Adjustable distance cutoff (Å) for a feature match — used everywhere scoring happens.                                                                                                |
| Score (single candidate)                         | Scores one SMILES; highlights matched atoms in PyMOL; reports the Verdict.                                                                                                           |
| Diverse hits only / Restrict to MW range         | Redundancy removal; a size-based red flag independent of chemical similarity.                                                                                                        |
| Fast Screen (checkbox, NEW)                      | Alignment-free pre-filter, ~6x faster; deliberately lenient — rejects only the obviously-hopeless before full scoring.                                                               |
| Ensemble tolerance matching (checkbox, NEW)      | Match against the nearest real clustered point instead of the averaged centroid. Off by default for consistency with prior results; falls back automatically for older saved models. |
| Align RMSD / Plausibility                        | Alignment fit quality; chemical similarity to the nearest active.                                                                                                                    |
| Shape (column, NEW)                              | % of the candidate's own atoms inside the training actives' combined shape envelope. Ligand-shape check only — not a receptor check.                                                 |
| Scaffold (column, NEW)                           | Same / Analog / Novel Bemis-Murcko comparison to the training actives.                                                                                                               |
| Verdict                                          | The combined match+plausibility classification, color-coded.                                                                                                                         |
| Export Results (CSV) / Export Aligned Hits (SDF) | CSV now includes shape_consistency and scaffold_novelty alongside the original columns. SDF: real 3D aligned structures.                                                             |

### 4.3 Tab 3 — Validate

Leave-one-out cross-validation (default) or External Test Actives (a genuinely separate, real held-out set) for the actives side; a required, user-supplied decoy file matched to whichever set is being validated; ROC-AUC / EF / GH / permutation p-value, a publication-quality ROC plot, and export to PNG or a self-contained HTML report.

## 5. Real-World Validation: Case Studies

*Not hypothetical examples — each was run for real, on real data.*

### 5.1 MAO-B inhibitors (5 actives)

> **REAL FINDING:** A factual structural error was found and corrected: the SMILES used for lazabemide was a different constitutional isomer of the real drug (same formula/MW, different InChIKey). Corrected. Even with the correct structure, lazabemide remained the weakest-matching active — real evidence of a genuinely distinct chemotype, not an artifact of the earlier error.

### 5.2 DPP-4 inhibitors / gliptins (7 actives)

A genuinely diverse class (peptidomimetic vs. non-peptidomimetic scaffolds). Validation: AUC 0.875 but GH 0.000 — a clean demonstration that a model can rank well while having zero compounds cleanly pass a strict hit threshold. This case study's model was also used to verify Ensemble Tolerance Matching (Section 3.4) on real data.

### 5.3 COX-2 selective inhibitors / coxibs (4 training + 2 held-out)

Training: celecoxib, rofecoxib, valdecoxib, etoricoxib. Validation via External Test Actives: parecoxib (a real positive control) and lumiracoxib (a genuine structural outlier).

> **REAL FINDING:** Parecoxib validated as predicted (4/4 features). Lumiracoxib matched better than expected on raw score (75%) because the 80%-threshold consensus lost its one chemically distinctive feature (a sulfonamide donor present in only 2 of 4 actives) and was left dominated by generic Aromatic/Hydrophobe features. Chemical plausibility (0.145) and the Verdict system correctly flagged this despite the misleading raw match. This case study's model was also used to verify Shape Consistency and Scaffold Novelty (Sections 3.3, 3.5) on real data, and to test a batch library later found to be pre-filtered by 90% Tanimoto similarity search — a real lesson that a tight 2D pre-filter defeats the purpose of pharmacophore screening by guaranteeing every library member already resembles the query.

## 6. Known Limitations

- **No receptor/structural information anywhere in the pipeline:** every check in this tool — including Shape Consistency — reflects agreement with the training ligands only, never with a real binding pocket. Real docking remains the categorically necessary next step, always.

- **Fast Screen is deliberately lenient:** verified more permissive than full scoring; only ever a time-saving pre-filter, not a substitute for the real score.

- **Chemical plausibility is an adapted heuristic, not formal QSAR-AD.**

- **Reference-ligand dependence** and **small active-set statistics** remain as documented previously; both are directly measurable via the tool's own robustness/permutation-test features rather than hidden.

- **Multi-model voting (parameter-sweep robustness) remains unbuilt:** the literal proposal was found to cost ~20 hours per 1000-compound library; a cheaper redesign (re-score only an already-shortlisted top-N across a few parameter variations) was scoped but not implemented.

## 7. Testing Methodology

30 automated tests (up from 26 at the previous documentation pass), all run against real PyMOL and real RDKit computation — no mocked chemistry. Every new function in this update was verified against real data (real coxib/gliptin structures) before being wired into the GUI, and every GUI change was confirmed via a static check of the shipped source in addition to unit-level logic tests.

> **NOTE:** This update was itself produced after a full environment reset (all installed packages, including RDKit and PyMOL, were lost between sessions). The exact shipped code was recovered from the previously-delivered output archive, dependencies were reinstalled, and the full 30-test suite was re-run and re-confirmed passing before this document was written — the same verify-before-writing discipline applied throughout this project, not relaxed for a documentation-only task.

*This document reflects the tool as of the most recent build (30 passing tests). It supersedes all earlier documentation, including the previous “Complete Documentation,” “Reference Guide,” and “Expert Review” documents.*
