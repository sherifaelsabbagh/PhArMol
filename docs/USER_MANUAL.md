# PhArMol User Manual

## 1. Introduction

This plugin helps you find new candidate molecules for a drug target when you already know a handful of active compounds against it. It works in three steps: it figures out the 3D pattern of chemical features your known actives share (a **pharmacophore**), it checks whether new candidate molecules — one at a time or by the thousands — share that same pattern, and it tells you how reliable that pattern actually is.

It runs **inside PyMOL**, using PyMOL's own 3D viewer — there is no separate application to open. You will need:

- PyMOL, with the RDKit library installed in the same Python environment (Section 2 explains the recommended way to set this up)

- At least 2–3 known active compounds against your target, as SMILES strings (a standard way of writing a molecule's structure as text)

- For validation and library screening: a set of “decoy” compounds (expected inactives) and/or a library file to screen, both described later in this manual

## 2. Installation & Launching

PyMOL and RDKit need to be installed into the **same** Python environment to work together — installing RDKit into your own separate Python (or with a plain pip install) is invisible to a PyMOL that was installed some other way. The recommended path below avoids this problem entirely by installing both together from the start.

### 2.1 Recommended: create one shared environment (do this first)

This is the path most likely to work on the first try, on any operating system, with no troubleshooting. Using conda/mamba:

```
conda create -n pharm -c conda-forge python=3.10 pymol-open-source rdkit
```

```
conda activate pharm
```

```
pymol
```

Or, using a plain Python virtual environment with pip instead of conda:

```
python3 -m venv pharm-env
```

```
source pharm-env/bin/activate \# Windows: pharm-env\Scripts\activate
```

```
pip install pymol-open-source rdkit
```

```
pymol
```

Either way, PyMOL and RDKit now live in the exact same Python installation, so there is nothing further to configure. From here, install the plugin via PyMOL's **Plugin → Plugin Manager → Install New Plugin**, selecting the plugin's folder or ZIP file, then open it from the **Plugin** menu.

### 2.2 Fallback: if you already have PyMOL installed another way

If you have an existing PyMOL installation (from Homebrew, an official installer, or a Linux package) and would rather not create a fresh environment, RDKit needs to be installed into **that specific PyMOL's own bundled Python** instead. First, find out which interpreter that is by typing this into PyMOL's own command line:

```
import sys; print(sys.executable)
```

Then install RDKit into that exact interpreter (not your regular system Python):

"\<the path just printed\>" -m pip install rdkit

This works, but is more fragile than Section 2.1 — that bundled Python may use an older pip or a different Qt binding. If the plugin still cannot find RDKit after this, it will show a popup telling you exactly which interpreter PyMOL is using and the exact command to run, rather than failing with a raw error.

With RDKit available either way, you can load the plugin manually from PyMOL's command line if the Plugin Manager does not work for you:

```
import sys; sys.path.insert(0, "/path/to/pharmol_folder")
```

```
import pharmol; pharmol.run_plugin_gui()
```

> **IMPORTANT:** Always try Section 2.1 (one shared environment) first. Only use the fallback in Section 2.2 if you specifically need to keep using an existing PyMOL installation.

## 3. Overview of the Interface

The plugin window has three tabs, meant to be used in order:

| Tab | What you do here |
|-----------------|--------------------------------------------------------------------------------------|
| 1\. Build Model | Enter your known active compounds and build a 3D consensus pharmacophore from them.  |
| 2\. Screen      | Check a single candidate, or an entire compound library, against that pharmacophore. |
| 3\. Validate    | Test how well the model actually tells real actives apart from inactive compounds.   |

## 4. Quick Start

The fastest way to see the tool work end-to-end:

**1.** Open Tab 1 (**Build Model**) and click **Load MAO-B Example** — this fills the input box with 5 real, known MAO-B inhibitor drugs, purely as a working demonstration.

**2.** Click **Analyze → Load into PyMOL**. After a short wait, the aligned molecules and the consensus pharmacophore (colored spheres) will appear in PyMOL's own 3D viewer, and the feature table will fill in.

**3.** Switch to Tab 2 (**Screen**), type any SMILES into the **Screen a Candidate** box, and click **Score** to see how well it matches.

**4.** To validate the model, go to Tab 3 (**Validate**), choose a decoy file (see Section 7.3 for what this file should contain), and click **Run Validation**.

The rest of this manual documents every input, parameter, and output on all three tabs in full detail.

## 5. Tab 1: Build Model

### 5.1 Inputs

#### Known Actives

A text box where you list your known active compounds, one per line, in the format:

name, SMILES

For example:

aspirin, CC(=O)Oc1ccccc1C(=O)O

You need at least 2 compounds, though 5 or more is strongly recommended for a meaningful result. **Load MAO-B Example** fills this box with a working demonstration set, useful for learning the tool before using your own data.

### 5.2 Parameters

| Parameter | Controls | Default / Range | Guidance |
|-------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------|--------------------------------------------------------------------------------------------------------------------------------|
| Conformers per molecule       | How many different 3D shapes are generated and tested per molecule.                                                                                                                             | 50 (10–200)         | Higher = more thorough but slower. Increase for very flexible molecules; the default is fine for most cases.                   |
| Protonate at physiological pH | Whether aliphatic amines are given a positive charge and carboxylic acids a negative charge before building 3D structures, mimicking pH ~7.4 (body conditions).                                 | On                  | Leave on unless you specifically need neutral (uncharged) structures.                                                          |
| Match tolerance               | How close (in Å, a unit of length) a candidate's feature must land to a consensus point to count as a match. Used on Tabs 2 and 3 as well.                                                      | 1.8 Å (0.5–4.0)     | Lower = stricter matching. Raise it if almost nothing matches; lower it if almost everything does.                             |
| Support threshold             | The minimum percentage of your actives that must share a feature for it to count as part of the final consensus pharmacophore.                                                                  | 50%                 | Raise it for a smaller, more universally-shared set of features; lower it to include more borderline, less universal features. |
| Clustering radius (eps)       | How close (in Å) two molecules' features need to be in 3D space to be grouped together as “the same” feature by the clustering algorithm (DBSCAN — see Section 10) when building the consensus. | 1.5 Å (0.8–2.5)     | See the detailed guidance below.                                                                                               |

> **TIP:** Clustering radius (eps), in practice: if your actives are all rigid, closely related structures (e.g. fused ring systems with little flexibility), try a tighter radius like 1.2 Å — their features should genuinely overlap closely in 3D, and a tight radius avoids blurring together features that are actually in slightly different places. If your actives are flexible, acyclic, or otherwise structurally diverse, try a looser radius like 1.8 Å, since real shared features may land in a wider spread of positions across different molecules. The default, 1.5 Å, works well as a starting point for most cases.

> **TIP:** Changing Support threshold or Clustering radius updates the results immediately, without needing to re-click Analyze.

### 5.3 Action: Analyze → Load into PyMOL

Builds the model: aligns every active in 3D onto one reference active, clusters their shared features using DBSCAN (Section 10), and loads the aligned molecules plus the consensus pharmacophore spheres directly into PyMOL's viewer.

### 5.4 Outputs

#### The Feature Table

| Column | Meaning |
|------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Feature    | The chemical type: Donor/Acceptor (hydrogen bonding), Aromatic (a flat ring), Hydrophobe (a greasy region), PosIonizable/NegIonizable (likely charged group), ZnBinder (can bind a zinc ion).                                                                                                                |
| Support    | The percentage of your active compounds that share this feature.                                                                                                                                                                                                                                             |
| Points     | How many individual atom-level feature detections, across all molecules, contributed to this one entry. These are the same individual points averaged together into this feature's 3D position, and also what Ensemble tolerance matching (Section 6.6) checks against directly instead of just the average. |
| Stability  | Filled in only after running the robustness check below — how often this feature survives when the model is rebuilt with one active removed at a time.                                                                                                                                                       |
| Consensus? | A checkmark if this feature's Support meets the current Support threshold — i.e. whether it is actually part of the active pharmacophore right now.                                                                                                                                                          |

#### 3D Viewer (in PyMOL)

Your aligned active molecules (as colored sticks) and the consensus pharmacophore (as colored spheres, one color per feature type) appear directly in PyMOL's own viewport, fully rotatable and zoomable with your mouse as usual.

### 5.5 Other buttons on this tab

| Button | What it does |
|------------------------------------|---------------------------------------------------------------------------------------------------------|
| Run Leave-One-Out Robustness Check | Rebuilds the model once per active, each time leaving one out, and fills in the Stability column above. |
| Zoom to Pharmacophore              | Centers PyMOL's camera on the model.                                                                    |
| Toggle Ligands                     | Shows/hides the aligned active molecules (the spheres stay visible either way).                         |
| Save Model...                      | Saves the current model — including full 3D data — to a .json file you can reopen later.                |
| Load Model...                      | Reopens a previously saved .json model file.                                                            |
| Compare Models...                  | Opens a side-by-side comparison of two saved models.                                                    |

## 6. Tab 2: Screen

### 6.1 Screen a Candidate

Paste a single SMILES into the box and click **Score**. The candidate is loaded into PyMOL (with the atoms that matched each feature highlighted in that feature's own color), and a verdict is shown, made up of:

| Output | Meaning |
|-----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Matched X/Y           | How many of the Y active consensus features this candidate satisfied.                                                                                        |
| Mean fit distance     | The average distance (Å), across only the matched features, between the candidate and the consensus positions. Lower is a tighter fit.                       |
| Alignment RMSD        | How well the candidate's overall shape fit the reference molecule during 3D alignment (Å), via an algorithm called O3A (Section 10). Lower is more credible. |
| Chemical plausibility | How structurally similar (0–1) the candidate is to the single most similar known active.                                                                     |
| Verdict               | A plain-language summary combining the match and plausibility scores — see Section 6.4.                                                                      |

### 6.2 Batch Screen a Library

Click **Choose SDF File...** to select a library of compounds (an SDF file — a standard multi-molecule chemical file format), then set any of the options below before clicking **Screen Library**.

| Option | What it does |
|-----------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Diverse hits only                                   | Removes near-duplicate compounds from the results, keeping only the best-scoring member of each group of very similar hits.                                                                                                                                                                                                                                                                                                                                                                                             |
| Restrict to MW range of training actives ± \[X\] Da | Removes any candidate whose molecular weight falls outside your active compounds' own weight range (widened by however many Daltons you set).                                                                                                                                                                                                                                                                                                                                                                           |
| Fast Screen                                         | Speeds up screening of very large libraries. It uses a rapid 3D distance fingerprint — the pairwise distances between your consensus features — to reject obvious non-matches before running the full, more precise O3A alignment on the rest. Because this fingerprint check is deliberately quick and approximate, it only rejects candidates that clearly cannot match; it never replaces the full, precise score for anything that passes it. Recommended only for libraries of several thousand compounds or more. |
| Ensemble tolerance matching                         | A more forgiving way of measuring feature matches — see Section 6.6 for a full explanation. Off by default.                                                                                                                                                                                                                                                                                                                                                                                                             |

### 6.3 Batch Results Table

One row per compound, sorted best-first. Click any row to load that compound into PyMOL.

| Column | Meaning |
|--------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Name               | The compound's name or ID, taken from the library file.                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| Matched / Fraction | How many consensus features matched, and as a percentage.                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Mean Dist (Å)      | Average distance for the matched features.                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Align RMSD (Å)     | Overall shape-fit quality during alignment.                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Plausibility       | Structural similarity to the most similar known active (0–1).                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| Shape              | What percentage of the candidate's own atoms fit within the overall size/shape of your known active compounds — see Section 6.5 for a full explanation. As a rule of thumb, a Shape score below about 60% suggests the candidate carries significant extra bulk that none of your actives possess — a strong reason to deprioritize it, especially if it is also flagged with a ‘Scaffold-Hop (High Risk)’ verdict, since the two together mean the compound neither looks like your actives chemically nor physically. |
| Scaffold           | Whether the candidate's core ring structure is the Same as, an Analog of, or entirely Novel compared to your known actives' own ring systems — see Section 6.7.                                                                                                                                                                                                                                                                                                                                                         |
| Verdict            | See Section 6.4.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |

### 6.4 Understanding the Verdict

The Verdict combines the match score and the plausibility score into one label:

| Verdict | Meaning |
|------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Gold Standard Hit                  | Matches the pharmacophore well AND is chemically similar to a known active — the strongest kind of result.                                                                                                                                  |
| Scaffold-Hop (Moderate Confidence) | Matches well, moderately similar chemically — plausibly a real hit with a different scaffold.                                                                                                                                               |
| Scaffold-Hop (High Risk)           | Matches well but looks chemically unlike any known active. Could be a genuine new-scaffold discovery, or a coincidence. Worth further checking (e.g. size, Shape score, or docking) before trusting it — not to be discarded automatically. |
| Uncertain                          | Only a partial match either way; not strong evidence for or against.                                                                                                                                                                        |
| Inactive Analog                    | Chemically similar to a known active but missing the key 3D feature arrangement.                                                                                                                                                            |
| Alien / Low Priority               | Neither the shape nor the chemistry resembles your known actives. Lowest priority.                                                                                                                                                          |

> **IMPORTANT:** A 'Scaffold-Hop (High Risk)' verdict is a flag for extra scrutiny, not a rejection. Finding real hits with different chemical scaffolds from your known actives is one of the main reasons to use this kind of tool in the first place.

### 6.5 Understanding the Shape column in detail

The pharmacophore match only checks a handful of individual 3D points — it says nothing about the *rest* of a candidate's structure. A molecule could satisfy every consensus feature perfectly while also carrying a large substituent none of your real actives ever had, and the match score alone would never catch that. The Shape score exists specifically to catch this.

It works by first combining all of your aligned known actives into one overall 3D silhouette — treating every atom as a small sphere sized to its true atomic radius, and taking the union of all those spheres, across all actives, as the reference shape. Then, for each candidate (once it has been aligned onto the same reference), every one of its own atoms is checked: does it fall inside that combined silhouette, or does it stick out past every actual atom your training actives ever had at that position? The Shape score is simply the percentage of the candidate's atoms that stay inside.

A candidate that is exactly the same size and shape as your training actives will score close to 100%. A candidate with one bulky extra ring or a long chain none of your actives had will score visibly lower, in proportion to how much of it sticks outside the combined envelope.

> **IMPORTANT:** The Shape score reflects consistency with your training ligands' own shapes only. It is not, and cannot be, a check against a real protein binding pocket — this tool has no information about your target's actual 3D structure at all. A high Shape score means ‘shaped like my known actives,’ not ‘fits the binding site.’

### 6.6 Understanding Ensemble tolerance matching in detail

Each consensus feature is not really one single position — it is built by averaging together several individual detections, one from each of your active compounds that shares that feature (the “Points” count in the Feature Table, Section 5.4). By default, matching only checks a candidate's distance to that single *averaged* position.

Turning on Ensemble tolerance matching changes this: instead of checking distance to the average, it checks distance to **whichever of the original individual detections is closest**. In other words, a candidate is now considered a match if it resembles how *any one* of your real active compounds actually presented that feature — not only if it resembles the mathematical average of all of them combined.

This matters most when your active compounds are structurally diverse: if three of your actives place a given feature in slightly different spots (because they are different shapes overall), the averaged centroid might not sit close to any of them individually, even though each one, on its own, is a perfectly good match. Ensemble tolerance matching recovers these cases. It is turned off by default so that results stay consistent with the plugin's original, simpler behavior unless you specifically choose to use it, and it works automatically with older saved models that do not have this extra data stored (it simply falls back to the averaged position for those).

### 6.7 Understanding the Scaffold column in detail

This checks the candidate's core ring system (technically, its Bemis-Murcko scaffold — the molecule with all side chains and substituents removed, leaving only the connected ring structure) against each of your known actives' own core ring systems.

| Value | Meaning |
|-----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Same Scaffold   | The candidate's ring system is chemically identical to one of your actives' — e.g. it may be a minor variant, a different salt form, or a prodrug of a known active. |
| Analog Scaffold | The candidate's ring system is structurally related, but not identical, to one of your actives'.                                                                     |
| Novel Scaffold  | The candidate's ring system does not resemble any of your actives' — a genuinely new chemotype.                                                                      |

This is a different question from Plausibility (which compares the *whole* molecule, not just the ring system) or Verdict (which combines match and plausibility). A candidate can be a Novel Scaffold while still scoring well on Plausibility, if its substituents happen to resemble your actives' even though its core ring system does not.

### 6.8 Exporting Results

| Button | What it produces |
|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| Export Results (CSV)      | A spreadsheet with every column from the results table, for every screened compound.                                                    |
| Export Aligned Hits (SDF) | The actual 3D aligned structure of every result, as a standard SDF file you can reopen directly in PyMOL or any other molecular viewer. |

## 7. Tab 3: Validate

This tab checks whether your pharmacophore model actually tells real active compounds apart from inactive ones — not just whether it looks reasonable.

### 7.1 Choosing what to validate

#### Use leave-one-out cross-validation for actives

On by default. Each active compound is tested only against a version of the model that was rebuilt **without** it, giving a fair, honest test rather than testing the model against the very compounds it was built from. Leave this on for the normal case, where the actives you are validating are the same ones you used to build the model on Tab 1.

#### External test actives (optional)

If you have a **separate** set of known actives that were *not* used to build the model, paste them here (same name, SMILES format as Tab 1). When this box has content, the model is validated against these compounds instead, and leave-one-out is switched off automatically (it isn't needed for compounds that were never part of training).

### 7.2 Decoy file (required)

Decoys are compounds **expected to be inactive**, needed as a point of comparison. This tool does not generate decoys for you — prepare a file (SDF, or a plain text file with one SMILES per line) using a dedicated decoy-generation tool or service, ideally one that matches your decoys to your actives' general size and properties while keeping them structurally distinct. Click **Choose Decoy File...** to select it.

### 7.3 Running validation and reading the results

Click **Run Validation**. The results report four numbers:

| Metric | Meaning | How to read it |
|------------------------|-----------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ROC-AUC                | Overall ranking quality, 0 to 1.                                                                          | 0.5 = no better than random; 1.0 = perfect. Above ~0.7 is generally considered good.                                                                                                                                 |
| Enrichment Factor (EF) | How many times richer in real actives the top-ranked slice of compounds is, compared to picking randomly. | Higher is better. An EF of 5× means the top of your ranked list has 5 times more real actives than random chance would.                                                                                              |
| Güner-Henry (GH) score | A single 0–1 score combining both how many actives you found and how “clean” your hit list is.            | Above ~0.6 is usually considered good. Can be low even when AUC is good — they measure different things (see box below).                                                                                             |
| Permutation p-value    | The statistical likelihood your result could have happened by chance.                                     | Below 0.05 is typically considered statistically significant. With very few active compounds, this number cannot be very precise, and the tool will tell you the smallest possible value your dataset could produce. |

> **IMPORTANT:** ROC-AUC/EF and the GH score can genuinely disagree, and both are shown because neither tells the whole story alone. A model can rank compounds well overall (good AUC) while still having very few compounds that cleanly pass a strict pass/fail cutoff (low GH). This does not mean one number is wrong — they are answering different questions.

### 7.4 The ROC Curve

A standard plot of true positive rate against false positive rate. The dashed diagonal line represents random chance; the further your curve bows toward the upper-left corner, the better your model's ranking.

### 7.5 Exporting

| Button | What it produces |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Export Plot (PNG)         | A high-resolution image of the ROC curve alone.                                                                                                                           |
| Export Full Report (HTML) | A single, complete, shareable file: your model's settings, the validation results, the ROC plot, and a full table of every compound tested — viewable in any web browser. |

## 8. File Format Reference

### 8.1 Inputs you provide

| File / Field | Format |
|--------------------------------------------------|--------------------------------------------------------------------------------------------------------------------|
| Known Actives / External Test Actives (text box) | Plain text, one compound per line: name, SMILES                                                                    |
| Library to screen (Tab 2)                        | A standard SDF file (one or more molecules with 3D or 2D structure).                                               |
| Decoy file (Tab 3)                               | Either an SDF file, or a plain text file with one SMILES per line (optionally as “name, SMILES” or “name SMILES”). |
| Saved Model to reload                            | A .json file previously produced by this plugin's Save Model button.                                               |

### 8.2 Outputs this tool produces

| File | Contents |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Saved Model (.json)       | The full pharmacophore model: every feature (position, support, and the individual clustered points behind it), plus the aligned 3D structures of your active compounds and your build settings, so it can be exactly reloaded later. |
| Batch Results (.csv)      | One row per screened compound, with every column shown in the results table (Section 6.3).                                                                                                                                            |
| Aligned Hits (.sdf)       | The real, 3D-aligned structure of every batch result, directly reusable in PyMOL or other molecular software.                                                                                                                         |
| ROC Plot (.png)           | A print-quality image of the ROC curve.                                                                                                                                                                                               |
| Validation Report (.html) | A complete, self-contained summary of a validation run, viewable in any browser.                                                                                                                                                      |

## 9. Tips for Combining the Outputs

No single column tells the whole story on its own. A few practical combinations worth using:

- **Looking for close analogs of your known actives:** filter for Gold Standard Hit, Shape above ~70%, and Scaffold = Same or Analog.

- **Looking for genuinely novel chemotypes (a true scaffold-hop):** look at Scaffold-Hop entries specifically where Shape is *also* reasonably high (roughly 70% or above). A novel scaffold that is also a reasonable physical size/shape match is a much more credible lead than one that is both chemically unlike your actives *and* a poor shape match.

- **Deprioritizing quickly:** a low Shape score (below ~60%) combined with a ‘High Risk’ verdict is the least trustworthy combination — the candidate resembles your actives neither chemically nor physically, and likely only matched the pharmacophore points by coincidence.

- **Before trusting any single top hit:** check its Align RMSD as well — a high match fraction paired with a high RMSD suggests the matched points landed within tolerance somewhat by chance, even before looking at Shape or Plausibility at all.

## 10. Glossary

#### Pharmacophore

The 3D arrangement of chemical features (hydrogen bond donors/acceptors, charged groups, hydrophobic and aromatic regions) believed to be responsible for a molecule's biological activity.

#### SMILES

A standard, compact text format for writing a molecule's chemical structure, e.g. CCO for ethanol.

#### SDF

A standard chemical file format that can store one or many molecules, including their 3D structure and named properties.

#### Conformer

One specific 3D shape a molecule can fold into; most molecules can adopt many different conformers by rotating around their single bonds.

#### DBSCAN

The clustering algorithm this tool uses to group individual feature detections from different active molecules into one shared consensus feature, based on how close together they are in 3D space (controlled by the Clustering radius / eps parameter, Section 5.2). Unlike some other clustering methods, DBSCAN does not need to be told in advance how many groups to expect — it naturally finds groups of closely-packed points, which suits this task well since the right number of consensus features isn't known ahead of time.

#### O3A (Open3DAlign)

The algorithm used to align molecules in 3D onto a shared reference — both when building the model (aligning your actives to each other) and when scoring any candidate (aligning it onto the model's reference before checking its features).

#### RMSD

Root-mean-square deviation — a standard measure (in Å) of how different two 3D structures (or how well-aligned two molecules) are. Lower means more similar / better aligned.

#### Tanimoto similarity

A standard 0–1 measure of how structurally similar two molecules are, based on comparing their chemical fingerprints. 1.0 means identical; 0.0 means nothing in common.

#### Bemis-Murcko scaffold

A molecule's core ring system, with all side chains and substituents removed — a standard way of comparing molecules by their underlying skeleton rather than their full decoration. Used for the Scaffold column (Section 6.7).

#### Van der Waals radius

The effective “size” of an atom, used to represent it as a small sphere when computing how much physical space a molecule occupies. Used to build the shape envelope behind the Shape column (Section 6.5).

#### Decoy

A compound believed to be inactive against your target, used as a point of comparison to test whether your model can actually tell active compounds apart from inactive ones.

#### Å (Angstrom)

A unit of length equal to one ten-billionth of a meter — the natural scale for measuring distances between atoms.

*End of manual.*
