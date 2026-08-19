# Changelog

All notable changes to PhArMol are documented here, in reverse chronological order.
This project developed iteratively, driven throughout by real data from real drug
targets (MAO-B, DPP-4 inhibitors, COX-2 inhibitors, HIV protease inhibitors, and
fluoroquinolone antibiotics) rather than hypothetical test cases — several entries
below describe a real bug found through actual use, not just a feature added.

## [Unreleased]

### Fixed
- **Toggle Ligands button did nothing.** `cmd.toggle()`'s first argument is a
  PyMOL representation name (e.g. `"sticks"`), not an object/selection pattern —
  `cmd.toggle("ligand_*")` silently failed with "unknown representation" every
  time. Fixed by tracking visibility state explicitly and using `cmd.enable()` /
  `cmd.disable()`, which correctly support wildcard patterns. Verified with a
  real two-click GUI test.

### Changed
- Renamed the project to **PhArMol** (from "Ligand-Based Pharmacophore Modeler").


**Four features from a second round of expert review**, each independently verified
against real data before being adopted — including one case (the fingerprint
pre-filter) where testing revealed a real trade-off worth being upfront about:

- **Scaffold Novelty** (Bemis-Murcko), Tab 2: a new "Scaffold" column classifying each
  batch-screening hit as Same/Analog/Novel relative to the training actives' own core
  ring systems — a different question than the Verdict system (chemical similarity
  overall) or Plausibility (fingerprint similarity). Verified on the real coxib case
  study: parecoxib (a genuine prodrug of valdecoxib) correctly reduces to the *exact*
  same Murcko scaffold; lumiracoxib (a genuine structural outlier) correctly comes
  back "Novel."
- **Fast Screen pre-filter**, Tab 2: an optional, alignment-free geometric pre-check
  (~6× faster than full scoring) that rejects obviously-hopeless candidates before
  spending time on full O3A-based scoring. Deliberately scoped as a lenient pre-filter
  only, never a scoring replacement — real testing found it noticeably more permissive
  than full scoring (a genuine outlier scored 0.90 on the fingerprint alone vs. its
  real, much lower full score), since internal pairwise distances alone don't capture
  full 3D arrangement the way real alignment does.
- **Shape Consistency**, Tab 2: a new "Shape" column showing what fraction of a
  candidate's own heavy atoms fall within the training actives' combined shape
  envelope. Explicitly and deliberately scoped as a **ligand**-shape check only — not
  a receptor steric-clash check, since this pipeline has no receptor structure
  anywhere in it. Verified on real data: every training active scores 100% against
  its own envelope; a deliberately oversized, unrelated compound dropped to 62%.
- **Ensemble tolerance matching**, Tab 2 (opt-in checkbox, off by default): consensus
  features now optionally store the real individual clustered point positions
  alongside their averaged centroid, and matching can check distance to the nearest
  actual point instead of collapsing everything to one fixed position — tolerating
  real positional spread across a flexible or diverse active set. A cheap redesign of
  a much more expensive "conformer ensemble" idea that was tested and found to cost
  56 hours on a 1000-compound library in its original form; this version adds no
  extra O3A alignments at all. Verified on real gliptin data (mean match distance
  0.78 → 0.13 Å for a genuine training active) and fully backward compatible — falls
  back to centroid matching automatically for any point loaded from an older saved
  model that doesn't have raw points stored.

**Combined Verdict classification (single candidate + batch screening).** Following
external expert feedback on how to properly combine the pharmacophore match score with
the chemical-plausibility score (previously two separate numbers a user had to mentally
cross-reference), scored candidates now get one plain-language verdict: **Gold Standard
Hit**, **Scaffold-Hop (Moderate Confidence)**, **Scaffold-Hop (High Risk)**,
**Uncertain**, **Inactive Analog**, or **Alien / Low Priority** — color-coded, with a
tooltip explanation. Match tiers reuse this plugin's existing thresholds (≥75% high,
≥40% medium) for consistency; plausibility tiers (≥0.5 high, ≥0.3 medium, below low)
include a deliberate, explicit "Medium" tier rather than forcing a continuous value
into an arbitrary side of a single cutoff. **Critically or deliberately never a hard
filter**: a high-match/low-plausibility candidate is labeled a higher-risk
"Scaffold-Hop" candidate, not discarded — pharmacophore modeling exists specifically
to find real hits with different scaffolds from the training actives, and a naive
plausibility cutoff would destroy that capability. Verified against the real MAO-B
false-positive case found earlier this project (a fused-ring compound matching 6/6
features via a coincidental buried amine, at only 0.14 plausibility) — correctly
classified as "Scaffold-Hop (High Risk)," not silently treated as a top hit, through
both direct unit tests and a live GUI test with real batch-screening data.

**Fixed: validation was testing actives on the same model they were used to build (a real
methodological flaw).** Previously, the Validate tab scored the exact same active compounds
that were used to build the pharmacophore against that same model — a classic "testing on
training data" problem. Since the consensus feature centroids are computed directly from
those actives' own aligned positions, this partially validates whether the model can
reproduce its own training data, not whether it generalizes. Concretely measured on the
course's 5-compound MAO-B example: in-sample AUC was 0.880 vs. a genuine leave-one-out AUC
of 0.840 — a real, if here modest, inflation that is not guaranteed to stay small for other
active sets. **Fixed**: Tab 3 now defaults to **leave-one-out cross-validation** for the
actives side — each active is scored only against a pharmacophore model rebuilt *without*
it (reusing the same rebuild machinery as the leave-one-out robustness check), while decoys
are still scored against the full model (no leakage concern there, since decoys were never
part of training either way). The old in-sample behavior remains available via an
unchecked checkbox, correctly scoped to its one legitimate use: validating a model against
a genuinely separate, curated active set that was never used to build it. Both the
in-app results and the exported HTML report clearly state which mode was used.

**Six improvements from an expert review of this plugin's documentation.** After an
external cheminformatics review of the tool's methodology and open design questions:

- **Terminology corrected**: "Applicability Domain similarity" renamed to "Chemical
  Plausibility" throughout the GUI, matching the reviewer's correction that formal
  QSAR-AD (leverage/Williams-plot machinery) is properly tied to a fitted regression
  model, not a raw geometric pharmacophore match — this tool's similarity check is an
  adaptation of the AD *principle*, not the standardized method itself.
- **Alignment RMSD** is now computed and shown alongside every match score (single
  candidate, batch screening, CSV/SDF export): a low fraction-matched with a low RMSD
  is a more credible near-miss than the same fraction-matched with a high RMSD, which
  suggests the matched features landed in tolerance somewhat by chance.
- **MW range filter** (Tab 2): an optional checkbox restricting batch-screening results
  to within ±X Da (default 50) of the training actives' own MW range. Directly
  motivated by a real finding: a 143-compound library's top pharmacophore hits were all
  found to be substantially larger than the training actives, and molecular size turned
  out to be a more diagnostic red flag than chemical-plausibility similarity alone for
  that specific failure mode.
- **Adjustable DBSCAN clustering radius (eps)**, Tab 1: previously a fixed 1.5 Å;
  now a live-updating slider (0.8–2.5 Å) next to the support threshold slider,
  re-clustering the already-aligned actives in well under a second per change (no
  re-alignment needed) — turning a hardcoded assumption into an explorable parameter.
- **Contextual guidance after the leave-one-out robustness check**: if any high-support
  ("core") feature's stability drops below 60%, the tool now explicitly suggests
  splitting the active set into more chemically consistent groups or manually choosing
  a different alignment reference — the reviewer-endorsed alternative to automatic
  multi-hypothesis detection (attempted and abandoned earlier; see below).
- **Export Aligned Hits (SDF)**, Tab 2: batch-screening results can now be exported as
  a real SDF file with each hit's actual aligned 3D conformer (not just its SMILES),
  so a promising hit can be opened directly in PyMOL to see exactly how it aligned to
  the reference.

Also fixed in this pass: the Güner-Henry score's "n/a" display (when no compound
clears the hit threshold) now explains *why* ("no compound cleared the hit threshold
at this support level — try lowering it") instead of showing a bare, unexplained value.

**Fixed: candidates were scored without proper 3D alignment (a real, consequential bug).**
Previously, `score_molecule()` independently re-embedded a fresh 3D conformer for any
candidate being scored (single scoring, batch screening, *and* validation) and compared
its raw feature positions directly against the consensus pharmacophore's centroids --
with no alignment step tying that fresh conformer's coordinate frame to the one the
model was actually built in. This happened to work reasonably well for most compounds
(RDKit's default embedding often lands in a roughly similar orientation for small,
fairly rigid molecules), but was never guaranteed, and was caught via real use of this
plugin: lazabemide's properly-aligned conformer (visible in PyMOL, built during model
construction) cleanly matched Donor and PosIonizable at ~0.2 Å each, while
`score_molecule()`'s independently re-embedded, unaligned conformer matched a
*different* pair of features (Donor + Aromatic) entirely. Every score, batch-screening
result, and validation metric computed by earlier versions of this plugin was affected
by this. **Fixed**: `score_molecule()`, `validate_pharmacophore()`, and
`screen_library()` now all accept a `reference_mol`/`reference_conf_id` (obtained via
the pharmacophore's new `get_alignment_reference()` method, available on both a
freshly-built model and one reloaded from a saved file) and properly O3A-align each
candidate conformer onto that reference before comparing feature positions -- the same
alignment the model itself was built with. Verified as a direct regression test using
the real lazabemide case (2/6 matched, wrong feature pair → 3/6 matched, correct
features recovered), both at the `core.py` level and through a real, live GUI button
click.

**Fixed: the leave-one-out robustness check was silently coupled to the display
threshold.** An earlier version only evaluated stability for features that already
cleared whatever threshold the Tab 1 slider happened to be set to -- meaning the
lower-support, more borderline features (exactly the ones a robustness check is most
useful for) were silently skipped rather than tested. `leave_one_out_stability()` now
always evaluates every feature the model finds, independent of the GUI's display
threshold; the threshold is applied separately, only when deciding which features to
show as "consensus" in the table.

**Fixed: the built-in MAO-B example used an incorrect lazabemide structure.** The
SMILES used throughout this course for lazabemide (`Clc1ccncc1C(=O)NCCN`) turned out to
be a different constitutional isomer of the real drug -- same molecular formula
(C8H10ClN3O) and molecular weight, but a different InChIKey than the genuine compound
(verified against Sigma-Aldrich and ChemSpider). The built-in "Load MAO-B Example" and
the test suite now use the correct structure (`NCCNC(=O)c1ccc(Cl)cn1`). Worth noting
directly: even with the corrected structure, lazabemide still shows up as the weakest
match among the five actives -- real evidence (not a data-entry artifact) that it
represents a genuinely distinct chemotype from the other four propargylamine-based
actives, worth treating as a candidate second binding mode rather than diluting one
forced consensus across all five.

**Statistical significance for the validation metrics (permutation test).** Validation
now reports a p-value alongside AUC: labels are shuffled 2000 times (scores held fixed)
and the real AUC is compared against that random-chance distribution. Critically, the
result also states the *smallest possible p-value your dataset size can produce* — with
5 actives and 5 decoys there are only C(10,5) = 252 distinct label arrangements, so a
p-value finer than ~1/252 isn't meaningful, and the tool says so explicitly rather than
implying false precision. Verified to correctly separate a genuinely informative test
case (p ≈ 0.0005) from a random/uninformative one (p ≈ 0.53) before being adopted.

**Leave-one-out robustness check (Tab 1).** Rebuilds the consensus pharmacophore N
times, each time leaving one active out, and reports a "Stability" column alongside
each feature's support: what fraction of those reruns still contain that same feature.
A feature at 100% support in the full model but low stability is being driven by one
or two specific compounds, not a genuinely shared pattern — this gives the
"reference-ligand dependence" limitation already documented below actual empirical
teeth instead of leaving it as a warning in prose. On the real MAO-B example, this
surfaced a real, previously-invisible result: several features at 100% support in the
full model turned out to only be 80% stable under leave-one-out testing.

**Applicability domain flag**, shown both when scoring a single candidate (Tab 2) and
in every row of a batch-screening result (new "AD Sim" column, included in the CSV
export): the candidate's fingerprint similarity to the nearest known active. A
candidate can satisfy the 3D pharmacophore geometry by coincidence while looking
nothing like any real active the model was built from — a known pitfall of
pharmacophore matching — and this makes that distinguishable at a glance rather than
silently conflating "good geometric match" with "chemically reasonable."

**Diverse hits only (Tab 2 batch screening).** A checkbox that filters an already-ranked
result list down to structurally diverse hits, keeping only the best-scoring
representative of each cluster of near-duplicate analogs (Tanimoto > 0.6) rather than
letting minor substituent variants on one scaffold dominate the top of the list.

**Compare Models... (Tab 1).** Opens a small sub-dialog to load two saved `.json`
models side by side and see which consensus features are shared vs. unique to each —
e.g. comparing a model built with 30 vs. 50 conformers, or built from a different
active subset, to see whether that choice actually changed the conclusions.

**Export Full Report (HTML), Tab 3.** A single self-contained HTML file combining
everything from a validation run: the model's build parameters, the validation set
size, all metrics (AUC, EF, GH, and the permutation p-value) with their explanatory
text, the full publication-style ROC plot (embedded inline as base64), and a complete
per-compound results table. Verified by actually rendering the generated HTML to an
image and inspecting it, not just checking the file exists.

**Diverse hits only (batch screening, Tab 2).** A checkbox that filters an already-ranked
batch-screening result list down to structurally diverse hits: walking the ranked list
best-first, a hit is kept only if it isn't too similar (Morgan/Tanimoto > 0.6) to any
hit already kept. This directly addresses the common failure mode where the top of a
ranked hit list is dominated by minor substituent variants on one scaffold (e.g. a
halogen series) rather than genuinely distinct candidates — the best-scoring member of
each near-duplicate cluster is kept, the redundant rest are dropped. Verified on
synthetic test data (a 4-member near-duplicate cluster plus 2 genuinely different
scaffolds) and again through the real GUI checkbox before shipping.

**Export Full Report (HTML), Tab 3.** A single self-contained HTML file combining
everything from a validation run: the model's build parameters (active ligand names,
conformer count, protonation setting, support threshold), the validation set size,
all three metrics (AUC, EF, GH) with their explanatory text, the full publication-style
ROC plot (embedded inline as base64 — no separate image file to lose track of), and a
complete per-compound results table (every active and decoy, sorted, with SMILES,
match count, fraction, and mean fit distance). Verified by actually rendering the
generated HTML to an image and inspecting it, not just checking the file exists.

**Physiological-pH protonation (on by default, toggleable).** Before 3D embedding,
molecules are now passed through an explicit, conservative protonation step:
aliphatic amines (not amides, not aromatic-attached — anilines are correctly left
alone, since they're too weakly basic to protonate at pH ~7.4) are protonated;
carboxylic acids are deprotonated. This is a deliberately limited rule set — not a
full pKa predictor — documented and validated in `core.py` and `test_plugin.py`
against 10 real test molecules, including all 5 built-in MAO-B actives (each
correctly comes out with formal charge +1) plus edge cases (aniline, amides,
pyridine) confirmed to be correctly left untouched.

**Default conformer count raised from 15 to 50**, and now user-adjustable directly
in the dialog (10–200). Fifteen conformers risked missing the true bioactive shape
for more flexible molecules; fifty better matches common practice for ligand-based
pharmacophore work, at a still-fast ~4 seconds for 5 molecules in testing.

**Batch library screening (Section 4).** Point the plugin at any SDF file and every
molecule in it is scored against the current consensus pharmacophore, ranked
best-first, and shown in a results table. Results can be exported to CSV, or any row
clicked to load that specific compound into PyMOL next to the pharmacophore spheres
for visual inspection — the direct link between this tool and the similarity-search
and Pharmit-generated compound libraries built earlier in the course.

**Save / Load Model.** A built pharmacophore can now be saved to a small JSON file
(Section 2, "Save Model...") and reloaded later ("Load Model...") without re-running
alignment from scratch. The saved file includes the full, unfiltered set of consensus
feature points (so the support-threshold slider still works correctly after
reloading), build metadata (source ligand names/SMILES, conformer count, whether
protonation was used), and the original aligned ligand 3D structures (as MolBlocks,
verified to preserve exact coordinates and formal charges through the round-trip) so
reloading also repopulates PyMOL with the original aligned ligands, not just bare
spheres. If you load an older or hand-edited file that only has points and no saved
ligands, the plugin falls back gracefully to showing spheres only.

**Validation (Tab 3) — EF, ROC-AUC, Güner-Henry score, with a publishable ROC plot.**
A pharmacophore model is a hypothesis to be checked, not a final answer. This tab
tests whether the current model actually discriminates the known actives from a
decoy set, mirroring the validation approach from the course's Jupyter notebook:

- **Decoys are always supplied by you, as a file** (SDF or plain-text SMILES, one per
  line) — not auto-generated. Earlier versions of this plugin included automatic
  decoy generation (a DUD-E-style property-matched background pool, and a live
  ChEMBL-inactives query specific to MAO-B). Both were removed deliberately: dedicated
  external tools built specifically for decoy generation — LIDeB Tools, DeepCoy, or
  the DUD-E decoy server — do this job better than a generator hacked into a plugin
  ever could, and this also makes the plugin genuinely target-agnostic (an earlier
  version had MAO-B's ChEMBL ID hardcoded as a default, which was a real, unexamined
  limitation for anything other than the course's own case study). Generate your
  decoys with whichever proper tool suits your target, then load the file here.
- **ROC-AUC** is computed via the Mann-Whitney U statistic, implemented directly in
  NumPy rather than requiring scikit-learn (consistent with this plugin's other
  dependency choices) — verified to match `sklearn.metrics.roc_auc_score` exactly,
  including with tied scores, before being adopted.
- **GH score** uses a deliberately meaningful hit threshold (≥75% of features
  matched, not "matched at least one feature") — an earlier, looser threshold was
  found during development of this course's teaching notebook to make GH collapse
  toward its random-chance floor by counting nearly everything as a hit.
- **The ROC curve is rendered for actual publication use**, not just on-screen
  reference — a white background regardless of the app's dark theme, gridlines, full
  tick labels at 0.2 intervals, properly labeled axes ("False Positive Rate" / "True
  Positive Rate"), a title, and a legend distinguishing the model curve from the
  random-chance diagonal. **"Export Plot (PNG)"** saves it at 1200×800 for a paper,
  poster, or slide — rendered fresh at that resolution rather than just upscaling the
  on-screen widget.

