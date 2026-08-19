# PhArMol Tutorial: A Complete Worked Example (Fluoroquinolone Antibiotics)

## 1. Why This Example

Fluoroquinolones are a large, real family of antibiotics (ciprofloxacin, levofloxacin, and relatives) that all inhibit the same bacterial enzymes, DNA gyrase and topoisomerase IV, through the same core chemical mechanism. That makes them an unusually good teaching example: there are enough real, distinct compounds to build a genuinely large training set, and — as this tutorial will show — the result reveals something important about how to interpret pharmacophore validation results in general, not just for this target.

## 2. Step 1 — The Training Set

14 real, FDA-or-internationally-approved fluoroquinolones were used as known actives, pasted directly into Tab 1's

Known Actives

box:

- ciprofloxacin, levofloxacin, ofloxacin, norfloxacin, moxifloxacin, gatifloxacin, gemifloxacin, sparfloxacin, lomefloxacin, pefloxacin, enoxacin, fleroxacin, balofloxacin, pazufloxacin

These share a conserved core (a 4-oxo-quinoline-3-carboxylic-acid ring system with a fluorine substituent) but differ meaningfully in their other substituents — different amine ring systems, different N-1 groups — giving genuine chemical diversity for the model to work with, not an artificially uniform set.

## 3. Step 2 — Building the Model

Settings used: 50 conformers per molecule, protonation on, clustering radius (eps) 1.5 Å, support threshold 80%. At this threshold the model produced

10 consensus features

— a large, tightly-shared set, reflecting how conserved the core fluoroquinolone scaffold really is across all 14 training compounds.

## 4. Step 3 — Designing a Fair Validation

To honestly test the model, 4 real fluoroquinolones that were **never used to build it** were held out and pasted into Tab 3's **External Test Actives** box:

| Compound | Role in this test |
|------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Levonadifloxacin | A close structural analog of the training compounds — expected to validate well.                                                                                             |
| Enrofloxacin     | A close structural analog — expected to validate well.                                                                                                                       |
| Orbifloxacin     | A close structural analog — expected to validate well.                                                                                                                       |
| Ibafloxacin      | Deliberately chosen as the odd one out: it lacks the separate basic-amine ring every other compound in this tutorial has, and is built on a different, fused tricyclic core. |

A decoy file of 13 real, unrelated drugs (an NSAID, ACE inhibitors, a statin, several DPP-4 inhibitors, and others) was also prepared, deliberately chosen to fall in the **same molecular-weight range** as the fluoroquinolones (roughly 265–420 Da) — this matters, because comparing against decoys of a very different size would make the validation artificially easy for reasons that have nothing to do with the actual pharmacophore.

## 5. Step 4 — Running Validation

With External Test Actives filled in, leave-one-out is automatically skipped (these compounds were never part of training, so there's nothing to “leave out”). Results:

| Metric | Result |
|---------------------|---------------------------------------------------------------------------------------------------------|
| ROC-AUC             | 1.000 (perfect — every active ranked above every decoy)                                                 |
| Enrichment Factor   | 4.25×                                                                                                   |
| Güner-Henry score   | 1.000 (every compound that passed the strict cutoff was a real active, and every real active passed it) |
| Permutation p-value | 0.0005 — the smallest value this dataset size can produce                                               |

Every one of the 4 held-out actives matched all 10 consensus features. Including Ibafloxacin — the compound deliberately chosen to be structurally different.

## 6. Step 5 — Interpreting a “Perfect” Result

A perfect validation score is worth pausing on rather than just celebrating. Why did the deliberately different compound, Ibafloxacin, still match every feature?

Because for this drug class, the shared core **is** the mechanism, not an accident of which 14 compounds happened to be chosen. The carboxylic acid, the adjacent ketone, and the aromatic ring system are what allow a fluoroquinolone to chelate the magnesium ion required to bind DNA gyrase and topoisomerase IV — every real fluoroquinolone needs this arrangement to work at all. Ibafloxacin's missing amine ring turned out to only ever have contributed one or two of the ten features; the other eight or nine came from the mechanism-defining core that every fluoroquinolone, including Ibafloxacin, necessarily shares.

> **KEY LESSON:** How well a structurally different held-out compound validates tells you something real about the target: whether the training set's shared features are mechanism-defining (essential to how the drug class actually works, so even an unusual compound will still share them) or merely incidental (a pattern specific to the particular training compounds chosen, which a genuinely different compound might not share at all). Both outcomes are informative — neither is a tool failure. A perfect score here reflects real, conserved biology, not an error.

> **NOTE:** Do not expect this result on every target. Many drug classes achieve the same biological effect through more varied chemistry, where a structural outlier genuinely may not share the training set's consensus features. A perfect validation score is a discovery about the specific target, not a guarantee the tool will always produce one.

## 7. Step 6 — Screening a Real Library

The model was then used to screen 200 real compounds from a public screening library (ZINC), with no pre-filtering. The results form a natural companion to the validation: this checks not just “does the model recognize real fluoroquinolones” but “does it correctly reject compounds that aren't.”

| Result | Count | % of library |
|-------------------------------|-----------|------------------|
| Matched ≥90% of features      | 0         | 0%               |
| Matched ≥75% of features      | 4         | 2%               |
| Verdict: Uncertain            | 104       | 52%              |
| Verdict: Alien / Low Priority | 92        | 46%              |

Not one of the 200 library compounds came close to the training/test actives' scores — the model is highly specific, not just highly sensitive. The 4 compounds that did stand out are genuinely worth a second look: real ZINC compounds sharing the quinoline-carboxylic-acid core but substituting chlorine, bromine, or a trifluoromethoxy group for the classic fluorine-and-amine-ring pattern. Each was correctly labeled a **Novel Scaffold** with low chemical plausibility (0.18–0.23) — genuinely different chemistry — but a perfect Shape score (1.00), meaning each is physically the right size and shape despite looking chemically unlike the training set. This is the Verdict system doing its job: not auto-accepted (plausibility is low), not auto-rejected (the geometry and shape are both real matches) — flagged as exactly the kind of candidate worth a closer, manual look.

## 8. Key Takeaways

- **A large, real, chemically diverse training set is achievable** when a drug class has enough distinct approved and investigational compounds — 14 real fluoroquinolones were used here without compromising on verification.

- **Property-matched decoys matter.** Matching the decoy set's molecular weight range to the actives' own range avoids an artificially easy validation.

- **A deliberately different held-out compound is the single most informative thing you can add to a validation run** — it reveals whether your training set's shared features are mechanistically essential or just incidental to the specific compounds chosen.

- **A perfect AUC is a result to interpret, not just celebrate.** Understanding why it happened (here: a conserved, mechanism-defining core) is more valuable than the number itself.

- **Validation and library screening answer different, complementary questions:** one tests recognition of real actives, the other tests specificity against everything else. Both were needed here to get the full picture.

*All compounds and results in this tutorial are real and were independently verified (via InChIKey and molecular formula cross-checking) before use. The complete training, test-active, and decoy files used here are included alongside this tutorial.*
