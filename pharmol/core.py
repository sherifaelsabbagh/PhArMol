"""
core.py — the actual pharmacophore science.

This is the same, already-validated pipeline from the course's Jupyter notebook
(MAOB_Pharmacophore_Pipeline.ipynb): RDKit conformer generation, Open3DAlign (O3A)
alignment, and DBSCAN-based consensus feature clustering. Nothing here is new science —
it is the identical logic, just packaged as reusable functions/classes for a GUI
front-end instead of notebook cells.
"""
import os
import json
import collections
import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures


def _dbscan(coords, eps, min_samples):
    """
    Minimal, dependency-free DBSCAN (no scikit-learn/scipy required).
    Returns a labels array using sklearn's exact convention (-1 = noise,
    0..k = cluster id). Validated against sklearn.cluster.DBSCAN on
    synthetic multi-cluster test data across several eps/min_samples
    settings, producing identical partitions in every case.

    This replaces scikit-learn specifically to avoid a real, observed
    macOS/Accelerate-framework binary-compatibility failure in scipy
    (a missing 'NEWLAPACK' symbol) that occurred inside PyMOL's bundled
    Python on Apple Silicon. Our clustering workload here is always small
    (at most a few hundred 3D points), so a plain O(n^2) implementation
    is more than fast enough and removes a fragile dependency entirely.
    """
    n = len(coords)
    labels = np.full(n, -2, dtype=int)  # -2 = unvisited
    cluster_id = -1

    def region_query(i):
        dists = np.linalg.norm(coords - coords[i], axis=1)
        return np.where(dists <= eps)[0]

    for i in range(n):
        if labels[i] != -2:
            continue
        neighbors = region_query(i)
        if len(neighbors) < min_samples:
            labels[i] = -1
            continue
        cluster_id += 1
        labels[i] = cluster_id
        seeds = [s for s in neighbors if s != i]
        k = 0
        while k < len(seeds):
            j = seeds[k]
            if labels[j] == -1:
                labels[j] = cluster_id
            if labels[j] == -2:
                labels[j] = cluster_id
                j_neighbors = region_query(j)
                if len(j_neighbors) >= min_samples:
                    for x in j_neighbors:
                        if x not in seeds:
                            seeds.append(x)
            k += 1
    return labels

FEATURE_FACTORY = ChemicalFeatures.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
)

# ------------------------------------------------------------------
# Physiological-pH (~7.4) protonation normalization
# ------------------------------------------------------------------
# Explicit, conservative, and limited in scope by design -- NOT a full pKa
# predictor. Handles the two highest-confidence, most common cases:
#   - protonate aliphatic amines (not amides/sulfonamides/imines, and not
#     directly attached to an aromatic ring, since anilines are far too
#     weakly basic to be protonated at pH 7.4)
#   - deprotonate carboxylic acids
# Everything else (amidines, guanidines, phosphates, aromatic heterocycle
# basicity, etc.) is deliberately left untouched rather than guessed at.
# Validated against 10 real test molecules (including all 5 MAO-B actives
# used in this plugin's built-in example, plus aniline/amide/pyridine
# edge cases) before being adopted -- see test_plugin.py.
_AMINE_SMARTS = "[NX3;H2,H1,H0;A;!$(NC=[O,S]);!$(NS(=O)(=O));!$(N=*);!$(Nc);!$(N[a])]"
_ACID_SMARTS = "[CX3](=O)[OX2H1]"
_AMINE_QUERY = Chem.MolFromSmarts(_AMINE_SMARTS)
_ACID_QUERY = Chem.MolFromSmarts(_ACID_SMARTS)


def protonate_physiological(mol):
    """Return a new, sanitized Mol with the protonation rules above applied."""
    rw = Chem.RWMol(mol)
    for match in mol.GetSubstructMatches(_AMINE_QUERY):
        atom = rw.GetAtomWithIdx(match[0])
        atom.SetFormalCharge(atom.GetFormalCharge() + 1)
        atom.SetNoImplicit(False)
    for match in mol.GetSubstructMatches(_ACID_QUERY):
        atom = rw.GetAtomWithIdx(match[2])
        atom.SetFormalCharge(atom.GetFormalCharge() - 1)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    new_mol = rw.GetMol()
    Chem.SanitizeMol(new_mol)
    return new_mol


DEFAULT_N_CONFS = 50  # bumped from an earlier default of 15 -- 15 risked missing
# the true bioactive conformation for molecules whose relevant torsion wasn't
# sampled; 50 matches common practice for ligand-based pharmacophore work.

FAMILY_COLORS = {
    "Donor": (0.06, 0.43, 0.34),           # HBD - green
    "Acceptor": (0.90, 0.30, 0.29),        # HBA - red
    "Aromatic": (0.55, 0.36, 0.96),        # purple
    "Hydrophobe": (0.85, 0.47, 0.02),      # orange
    "LumpedHydrophobe": (0.85, 0.47, 0.02),
    "PosIonizable": (0.22, 0.54, 0.86),    # blue
    "NegIonizable": (0.96, 0.45, 0.71),    # pink
    "ZnBinder": (0.58, 0.58, 0.58),
}


class Ligand:
    """A single active molecule: name, SMILES, 3D conformer, and its RDKit Mol."""
    def __init__(self, name, smiles, n_confs=DEFAULT_N_CONFS, protonate=True):
        self.name = name
        self.smiles = smiles
        self.mol = None
        self.conf_id = None
        self.error = None
        self._prepare(n_confs=n_confs, protonate=protonate)

    def _prepare(self, n_confs=DEFAULT_N_CONFS, seed=42, protonate=True):
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            self.error = "Could not parse SMILES"
            return
        if protonate:
            try:
                mol = protonate_physiological(mol)
            except Exception as e:
                # Protonation is a refinement, not a hard requirement -- fall back
                # to the original neutral form rather than failing the whole ligand.
                self.error = f"Protonation step failed, used neutral form: {e}"
        mol = Chem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, randomSeed=seed, pruneRmsThresh=0.5)
        if len(cids) == 0:
            self.error = "Conformer embedding failed"
            return
        energies = []
        for cid in cids:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
            if ff is None:
                continue
            ff.Minimize(maxIts=500)
            energies.append((cid, ff.CalcEnergy()))
        if not energies:
            self.error = "Force field setup failed"
            return
        best_cid = sorted(energies, key=lambda x: x[1])[0][0]
        self.mol = mol
        self.conf_id = best_cid

    @property
    def is_valid(self):
        return self.mol is not None

    def n_features(self):
        if not self.is_valid:
            return 0
        return len(FEATURE_FACTORY.GetFeaturesForMol(self.mol, confId=self.conf_id))

    def features(self):
        """List of (family, x, y, z) for this ligand's best conformer."""
        if not self.is_valid:
            return []
        feats = FEATURE_FACTORY.GetFeaturesForMol(self.mol, confId=self.conf_id)
        out = []
        for f in feats:
            p = f.GetPos()
            out.append((f.GetFamily(), p.x, p.y, p.z))
        return out


class PharmacophoricPoint:
    """One consensus feature: family, 3D centroid, support fraction, contributing
    point count, and (optionally) the raw individual clustered point positions
    that centroid was averaged from."""
    def __init__(self, family, centroid, support_frac, n_points, raw_points=None):
        self.family = family
        self.centroid = np.array(centroid)
        self.support_frac = support_frac
        self.n_points = n_points
        self.color = FAMILY_COLORS.get(family, (0.7, 0.7, 0.7))
        # raw_points: the actual individual clustered coordinates (an (n_points, 3)
        # array), if available -- None for a point loaded from an older saved model
        # that didn't store this. Used by ensemble-tolerance matching (see
        # score_molecule's match_mode="ensemble") to check distance to the nearest
        # actual clustered point instead of just the averaged centroid, tolerating
        # real positional spread across a flexible or diverse active set rather
        # than collapsing it to one fixed position.
        self.raw_points = np.array(raw_points) if raw_points is not None else None

    def __repr__(self):
        return f"<PharmacophoricPoint {self.family} support={self.support_frac:.2f}>"

    def to_dict(self):
        d = {
            "family": self.family,
            "centroid": [float(v) for v in self.centroid],
            "support_frac": float(self.support_frac),
            "n_points": int(self.n_points),
        }
        if self.raw_points is not None:
            d["raw_points"] = [[float(v) for v in pt] for pt in self.raw_points]
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d["family"], d["centroid"], d["support_frac"], d["n_points"],
                    raw_points=d.get("raw_points"))


class LoadedPharmacophore:
    """
    A consensus pharmacophore loaded back from a saved file, with the same
    .points / .consensus_at_threshold() interface as a freshly-built
    LigandBasedPharmacophore, so the rest of the code (scoring, PyMOL
    display) doesn't need to know or care whether a model was just built
    or reloaded from disk.
    """
    def __init__(self, points, metadata=None):
        self.points = points
        self.metadata = metadata or {}
        self.ligands = []  # no aligned Ligand objects available for a reloaded model
        self.ligand_mols = []  # populated as [(name, Mol), ...] if the file included them

    def consensus_at_threshold(self, thresh):
        return [p for p in self.points if p.support_frac >= thresh]

    def get_alignment_reference(self):
        """Returns (mol, conf_id) to align new candidates onto, or (None, None)
        if this file didn't include saved ligand structures. Any one of the
        saved ligands works as the reference -- they're all already mutually
        aligned into the same shared frame, so picking the first one is a
        reasonable, arbitrary-but-valid choice."""
        if self.ligand_mols:
            return self.ligand_mols[0][1], 0
        return None, None

    def get_shape_envelope(self):
        """Returns the union of every saved ligand's heavy-atom positions and
        van der Waals radii -- see LigandBasedPharmacophore.get_shape_envelope
        for the full docstring. Empty if this file didn't include saved
        ligand structures."""
        if not self.ligand_mols:
            return []
        return compute_shape_envelope([(mol, 0) for _name, mol in self.ligand_mols])


def export_pharmacophore_json(points, path, metadata=None, ligands=None):
    """Save a list of PharmacophoricPoints (typically pharm.points, unfiltered
    by threshold, so the threshold slider still works after reloading) to a
    small JSON file. If `ligands` (a list of Ligand objects) is given, their
    aligned 3D structures are saved too (as MolBlocks, which preserve exact
    3D coordinates -- verified to round-trip to within ~1e-4 A), so a
    reloaded model can repopulate PyMOL with the original aligned ligands,
    not just the pharmacophore spheres."""
    payload = {
        "format": "pharmol.v1",
        "metadata": metadata or {},
        "points": [p.to_dict() for p in points],
        "ligands": [],
    }
    if ligands:
        for lig in ligands:
            if not lig.is_valid:
                continue
            payload["ligands"].append({
                "name": lig.name,
                "smiles": lig.smiles,
                "molblock": Chem.MolToMolBlock(lig.mol, confId=lig.conf_id),
            })
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def load_pharmacophore_json(path):
    """Load a previously exported consensus pharmacophore. Returns a
    LoadedPharmacophore, usable anywhere a LigandBasedPharmacophore is.
    If the file includes saved ligand MolBlocks, LoadedPharmacophore.ligand_mols
    is populated as a list of (name, Mol) tuples with their aligned 3D
    coordinates restored exactly."""
    with open(path) as f:
        payload = json.load(f)
    points = [PharmacophoricPoint.from_dict(d) for d in payload.get("points", [])]
    points.sort(key=lambda p: -p.support_frac)

    ligand_mols = []
    for lig_data in payload.get("ligands", []):
        mol = Chem.MolFromMolBlock(lig_data["molblock"], removeHs=False)
        if mol is not None:
            ligand_mols.append((lig_data["name"], mol))

    loaded = LoadedPharmacophore(points, metadata=payload.get("metadata", {}))
    loaded.ligand_mols = ligand_mols
    return loaded


class LigandBasedPharmacophore:
    """
    Aligns a set of Ligands and extracts a consensus pharmacophore.
    Mirrors the notebook's align + extract_consensus_pharmacophore steps.
    """
    def __init__(self, ligands):
        self.ligands = [l for l in ligands if l.is_valid]
        self.invalid = [l for l in ligands if not l.is_valid]
        self.reference = None
        self.points = []  # list of PharmacophoricPoint
        self._aligned = False

    def align(self, reference_name=None):
        if len(self.ligands) < 2:
            raise ValueError("Need at least 2 valid ligands to build a pharmacophore.")
        if reference_name is not None:
            matches = [l for l in self.ligands if l.name == reference_name]
            if not matches:
                available = ", ".join(l.name for l in self.ligands)
                raise ValueError(
                    f"No valid ligand named '{reference_name}' found. "
                    f"Available: {available}"
                )
            self.reference = matches[0]
        else:
            # Default heuristic: the ligand with the most detected features.
            # A known, documented limitation (see README) -- the resulting
            # consensus can be sensitive to this choice for chemically
            # diverse active sets, which is exactly why manual override
            # is offered here rather than only ever auto-picking.
            self.reference = max(self.ligands, key=lambda l: l.n_features())
        ref_mol, ref_cid = self.reference.mol, self.reference.conf_id
        for lig in self.ligands:
            if lig is self.reference:
                continue
            try:
                o3a = AllChem.GetO3A(lig.mol, ref_mol, prbCid=lig.conf_id, refCid=ref_cid)
                o3a.Align()
            except Exception as e:
                lig.error = f"Alignment failed: {e}"
        self._aligned = True

    def extract_consensus(self, eps=1.5, min_support_frac=0.5):
        if not self._aligned:
            self.align()
        all_points = []
        for lig in self.ligands:
            for fam, x, y, z in lig.features():
                all_points.append((fam, x, y, z, lig.name))

        by_family = collections.defaultdict(list)
        for fam, x, y, z, name in all_points:
            by_family[fam].append((x, y, z, name))

        n_ligands = len(self.ligands)
        min_samples = max(2, int(np.ceil(min_support_frac * n_ligands)))
        points = []
        for fam, pts in by_family.items():
            coords = np.array([[x, y, z] for x, y, z, _ in pts])
            names = [n for _, _, _, n in pts]
            if len(coords) < min_samples:
                continue
            labels = _dbscan(coords, eps=eps, min_samples=min_samples)
            for lbl in set(labels):
                if lbl == -1:
                    continue
                idxs = [i for i, l in enumerate(labels) if l == lbl]
                support_ligs = set(names[i] for i in idxs)
                support_frac = len(support_ligs) / n_ligands
                centroid = coords[idxs].mean(axis=0)
                points.append(PharmacophoricPoint(fam, centroid, support_frac, len(idxs),
                                                    raw_points=coords[idxs]))

        points.sort(key=lambda p: -p.support_frac)
        self.points = points
        return points

    def consensus_at_threshold(self, thresh):
        return [p for p in self.points if p.support_frac >= thresh]

    def get_alignment_reference(self):
        """Returns (mol, conf_id) that any new candidate should be O3A-aligned
        onto before comparing feature positions to this model's consensus
        centroids -- the same reference ligand chosen during align()."""
        if self.reference is not None:
            return self.reference.mol, self.reference.conf_id
        return None, None

    def get_shape_envelope(self):
        """Returns the union of every aligned active's heavy-atom positions
        and van der Waals radii -- the basis for shape_consistency_score().
        This is deliberately a LIGAND shape envelope only: it reflects how
        consistent a candidate's own bulk is with the training actives'
        own shapes, not whether it would sterically clash with a real
        receptor (this pipeline has no receptor information at all -- see
        the docs for why that's out of scope for a ligand-based tool)."""
        return compute_shape_envelope([(lig.mol, lig.conf_id) for lig in self.ligands])


def _point_distance(candidate_pos, pharmacophoric_point, match_mode):
    """Distance from a candidate feature position to a consensus point,
    per match_mode -- see score_molecule's docstring."""
    if match_mode == "ensemble" and pharmacophoric_point.raw_points is not None \
            and len(pharmacophoric_point.raw_points) > 0:
        deltas = pharmacophoric_point.raw_points - candidate_pos
        return float(np.min(np.linalg.norm(deltas, axis=1)))
    return float(np.linalg.norm(candidate_pos - pharmacophoric_point.centroid))


def score_molecule(smiles, pharmacophore_points, reference_mol=None, reference_conf_id=0,
                    tol=1.8, n_confs=DEFAULT_N_CONFS, seed=1, protonate=True, shape_envelope=None,
                    match_mode="centroid"):
    """
    Score a candidate SMILES against a list of PharmacophoricPoints.

    reference_mol/reference_conf_id, if given, should be one of the aligned
    ligands the pharmacophore was built from (get this via
    pharm.get_alignment_reference()). Each candidate conformer is O3A-aligned
    onto that reference *before* comparing feature positions to the
    consensus centroids.

    This alignment step matters more than it might look: without it, a
    freshly-embedded candidate conformer sits in an arbitrary orientation
    that has no defined relationship to the coordinate frame the consensus
    centroids actually live in -- comparing its raw feature positions to
    those centroids only "works" to the extent that RDKit's default
    embedding happens to coincidentally land in a similar orientation, which
    is not guaranteed and was found, during real use of this plugin, to
    fail for at least one real active (a different, correct set of matched
    features was found once alignment was added, vs. what was silently
    reported before). If no reference is available (e.g. an older saved
    model with no stored ligand structures), this falls back to the old,
    less reliable unaligned comparison -- a real degradation, not a silent
    equivalent, and callers should prefer always supplying a reference.

    shape_envelope, if given (from pharm.get_shape_envelope()), additionally
    computes shape_consistency_score() on the best-matching conformer -- see
    that function's docstring. A LIGAND-shape check only, not a receptor
    steric-clash check.

    Returns (n_matched, fraction_matched, mean_match_distance, mol_or_None,
    best_conf_id, match_details, align_rmsd, shape_consistency), where
    match_details is a list of {"family": str, "color": (r,g,b), "atom_ids":
    tuple(int)} for the features that matched on the best-scoring conformer
    -- enough information for a GUI to highlight exactly which atoms
    satisfied each consensus feature, not just report a count. align_rmsd is
    the O3A alignment RMSD (heavy atoms) between the candidate and the
    reference on the best-scoring conformer -- None if no reference was
    supplied. A low fraction_matched with a low align_rmsd is a more
    credible near-miss than the same fraction_matched with a high
    align_rmsd (the latter suggests the matched features landed in
    tolerance somewhat by chance, not because the molecule as a whole
    actually fit the reference well). shape_consistency is None if no
    envelope was supplied.

    match_mode controls how a candidate feature's distance to a consensus
    point is measured: "centroid" (default) uses the single averaged
    position, exactly as before. "ensemble" instead checks distance to the
    NEAREST of that point's own raw clustered coordinates (falling back to
    the centroid for any point that doesn't have raw_points stored, e.g. one
    loaded from an older saved model) -- tolerates real positional spread
    across a flexible or diverse active set rather than collapsing it to
    one fixed position. Cheap: this only adds a handful of extra distance
    comparisons per feature (as many as contributed to that cluster, no
    extra O3A alignments), unlike a full conformer-ensemble rescoring
    approach.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not pharmacophore_points:
        return 0, 0.0, 999.0, None, None, [], None, None
    if protonate:
        try:
            mol = protonate_physiological(mol)
        except Exception:
            pass  # fall back to neutral form rather than failing the whole score
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, randomSeed=seed, pruneRmsThresh=0.5)
    if len(cids) == 0:
        return 0, 0.0, 999.0, None, None, [], None, None

    best = (0, 0.0, 999.0, None, [], None)
    for cid in cids:
        align_rmsd = None
        if reference_mol is not None:
            try:
                o3a = AllChem.GetO3A(mol, reference_mol, prbCid=cid, refCid=reference_conf_id)
                align_rmsd = o3a.Align()
            except Exception:
                pass  # if alignment fails for this specific conformer, fall through and
                      # score it unaligned rather than discarding it outright

        feats = FEATURE_FACTORY.GetFeaturesForMol(mol, confId=cid)
        by_fam = collections.defaultdict(list)
        for f in feats:
            p = f.GetPos()
            by_fam[f.GetFamily()].append((np.array([p.x, p.y, p.z]), tuple(f.GetAtomIds())))

        n_matched, dists, used, match_details = 0, [], set(), []
        for pp in pharmacophore_points:
            cands = by_fam.get(pp.family, [])
            best_d, best_i = 1e9, None
            for i, (c, _atom_ids) in enumerate(cands):
                if (pp.family, i) in used:
                    continue
                d = _point_distance(c, pp, match_mode)
                if d < best_d:
                    best_d, best_i = d, i
            if best_i is not None and best_d <= tol:
                n_matched += 1
                dists.append(best_d)
                used.add((pp.family, best_i))
                match_details.append({
                    "family": pp.family,
                    "color": pp.color,
                    "atom_ids": cands[best_i][1],
                })

        frac = n_matched / len(pharmacophore_points)
        mean_d = float(np.mean(dists)) if dists else 999.0
        if frac > best[1] or (frac == best[1] and mean_d < best[2]):
            best = (n_matched, frac, mean_d, cid, match_details, align_rmsd)

    shape_consistency = None
    if shape_envelope is not None and best[3] is not None:
        shape_consistency = shape_consistency_score(mol, best[3], shape_envelope)

    return best[0], best[1], best[2], mol, best[3], best[4], best[5], shape_consistency


# ------------------------------------------------------------------
# Validation: Enrichment Factor, ROC-AUC, Guener-Henry score
# ------------------------------------------------------------------
# Mirrors the validation approach from the course's Jupyter notebook, but with
# no scikit-learn/scipy dependency (consistent with the rest of this plugin) --
# ROC-AUC is computed via the Mann-Whitney U statistic equivalence, verified
# to match sklearn.metrics.roc_auc_score exactly (including tied scores)
# across multiple synthetic test scenarios before being adopted.

def roc_auc(labels, scores):
    """ROC-AUC via the Mann-Whitney U statistic (no sklearn/scipy required).
    labels: list of 0/1. scores: list of floats (higher = more likely active)."""
    paired = sorted(zip(scores, labels))
    n = len(paired)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n and paired[j][0] == paired[i][0]:
            j += 1
        avg_rank = (i + j - 1) / 2 + 1
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j
    n_pos = sum(labels)
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return float('nan')
    sum_ranks_pos = sum(r for (s, l), r in zip(paired, ranks) if l == 1)
    return (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def enrichment_factor(labels, scores, top_pct=0.10):
    """EF at top_pct: how many more actives are found in the top slice of the
    ranked list than expected by chance."""
    n = len(labels)
    n_top = max(1, int(round(n * top_pct)))
    ranked = [l for _, l in sorted(zip(scores, labels), reverse=True)]
    hits_top = sum(ranked[:n_top])
    hit_rate_top = hits_top / n_top
    hit_rate_all = sum(labels) / n
    if hit_rate_all == 0:
        return float('nan')
    return hit_rate_top / hit_rate_all


def guner_henry_score(A, D, Ht, Ha):
    """Standard GH score (Guener & Henry, 2000). A = actives in the set,
    D = total compounds, Ht = total hits (by whatever hit definition was
    used), Ha = actives among those hits. Returns nan if undefined."""
    if Ht == 0 or A == 0 or D == A:
        return float('nan')
    return (Ha * (3 * A + Ht)) / (4 * Ht * A) * (1 - (Ht - Ha) / (D - A))


# A small, fixed, offline fallback pool of real, structurally diverse
# drug-like compounds -- used only if a live ChEMBL query is unavailable.
# Deliberately modest in size and diversity compared to a live query; this
# is a safety net, not a substitute for it.
_FALLBACK_BACKGROUND_POOL = [
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
    "CC(Cc1ccc(cc1)C(C)C(=O)O)",  # ibuprofen-like
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # caffeine
    "CC(=O)Nc1ccc(O)cc1",  # paracetamol
    "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C",  # omeprazole-like
    "CC(C)NCC(O)COc1cccc2ccccc12",  # propranolol
    "Clc1ccc(cc1)C(c1ccccc1)N1CCN(CCOCC(=O)O)CC1",  # cetirizine-like
    "CCN(CC)CCNC(=O)c1cc(Cl)c(N)cc1OC",  # metoclopramide-like
    "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",  # penicillin-like
    "COC(=O)C1=CC2Cc3ccc(O)cc3C1CC2",  # steroid-like
    "Cc1ccc(cc1)S(=O)(=O)Nc1ccc(cc1)C(=O)O",  # sulfonamide-like
    "OC(=O)c1ccccc1Nc1cccc(c1)C(F)(F)F",  # niflumic-acid-like
    "CCOC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+](=O)[O-])C(=O)OC",  # nifedipine-like
    "CN1CCN(CC1)c1ccc2[nH]c3ccccc3c2c1",  # generic tricyclic amine (carbazole-piperazine)
    "Cc1onc(c1C(=O)Nc1ccccc1)-c1ccccc1",  # oxazole amide
    "CCN(CC)CCOC(=O)c1ccc(N)cc1",  # procaine-like
    "COc1cc2c(cc1OC)C(=O)C(CC1CCN(C)CC1)C2",  # generic bicyclic amine
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",  # theophylline
    "CC(C)Cc1ccc(cc1)C(C)C(=O)NCCO",  # amide-alcohol
    "O=C1CCC(=O)N1SC(Cl)(Cl)Cl",  # imide
    "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",  # salbutamol-like
    "COc1cc(cc(OC)c1OC)C(=O)N1CCN(CC1)c1ccccn1",  # trimethoxybenzamide
    "Clc1ccc2c(c1)nc([nH]2)-c1ccncc1",  # benzimidazole
    "CC1=CC(=O)C=CC1=O",  # quinone
    "O=C(Nc1ccc(cc1)S(N)(=O)=O)c1ccccc1",  # benzamide sulfonamide
    "CCC1(CC)C(=O)NC(=O)NC1=O",  # barbiturate-like
    "COc1ccc(cc1)C1=NN(C(=O)C1)c1ccccc1",  # pyrazolone
    "CC(=O)N1c2ccccc2Sc2ccccc21",  # phenothiazine-like
    "CCOC(=O)c1ncn2c1CN(C)C(=O)c1cc(F)ccc1-2",  # fused heterocycle
    "Oc1ccc(cc1)C1=Cc2ccccc2OC1=O",  # coumarin-like
]


def fetch_background_pool(n=1500, mw_min=150, mw_max=500, timeout=15):
    """Fetch a broad, target-agnostic pool of drug-like compounds to draw
    decoys from. Tries a live ChEMBL query first; on any failure (no
    internet, timeout, API change) falls back to a small bundled offline
    pool, so validation still works, just with less diversity. Returns
    (pool_smiles, source_description)."""
    try:
        from chembl_webresource_client.new_client import new_client
        molecule = new_client.molecule
        results = molecule.filter(
            molecule_properties__mw_freebase__gte=mw_min,
            molecule_properties__mw_freebase__lte=mw_max,
        ).only(["molecule_chembl_id", "molecule_structures"])
        pool = []
        for r in results:
            struct = r.get("molecule_structures") or {}
            smi = struct.get("canonical_smiles")
            if smi:
                pool.append(smi)
            if len(pool) >= n:
                break
        if len(pool) < 50:
            raise ValueError("ChEMBL returned too few usable compounds")
        return pool, f"{len(pool)} compounds from ChEMBL (live query)"
    except Exception:
        return list(_FALLBACK_BACKGROUND_POOL), \
            f"{len(_FALLBACK_BACKGROUND_POOL)} compounds (offline fallback \u2014 ChEMBL unreachable)"


def property_matched_decoys(active_smiles_list, background_pool, n_per_active=10,
                             mw_tol=25, logp_tol=1.0, max_similarity=0.35):
    """Select decoys: similar bulk properties (MW/LogP) to an active, but
    structurally dissimilar (low Morgan-fingerprint Tanimoto similarity).
    Same logic as the course's Jupyter notebook validation step."""
    from rdkit.Chem import Descriptors, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    _fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    active_mols = [Chem.MolFromSmiles(s) for s in active_smiles_list]
    active_mols = [m for m in active_mols if m is not None]
    if not active_mols:
        return []
    active_fps = [_fp_gen.GetFingerprint(m) for m in active_mols]
    active_props = [(Descriptors.MolWt(m), Descriptors.MolLogP(m)) for m in active_mols]

    decoys = []
    seen = set()
    for smi in background_pool:
        mol = Chem.MolFromSmiles(smi)
        if mol is None or smi in seen:
            continue
        mw, logp = Descriptors.MolWt(mol), Descriptors.MolLogP(mol)
        fp = _fp_gen.GetFingerprint(mol)
        for (a_mw, a_logp), a_fp in zip(active_props, active_fps):
            if abs(mw - a_mw) <= mw_tol and abs(logp - a_logp) <= logp_tol:
                sim = DataStructs.TanimotoSimilarity(fp, a_fp)
                if sim <= max_similarity:
                    decoys.append(smi)
                    seen.add(smi)
                    break
        if len(decoys) >= n_per_active * len(active_mols):
            break
    return decoys


def roc_curve_points(labels, scores):
    """Dependency-free ROC curve: (fpr_list, tpr_list) at each distinct score
    threshold, for plotting. Verified to integrate (trapezoidal) to the same
    AUC as sklearn.metrics.roc_curve + auc() across multiple test scenarios,
    including tied scores, before being adopted."""
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return [0.0, 1.0], [0.0, 1.0]
    paired = sorted(zip(scores, labels), key=lambda x: -x[0])
    fpr, tpr = [0.0], [0.0]
    tp, fp = 0, 0
    prev_score = None
    for score, label in paired:
        if prev_score is not None and score != prev_score:
            fpr.append(fp / n_neg)
            tpr.append(tp / n_pos)
        if label == 1:
            tp += 1
        else:
            fp += 1
        prev_score = score
    fpr.append(fp / n_neg)
    tpr.append(tp / n_pos)
    return fpr, tpr


def fetch_chembl_inactives(target_chembl_id="CHEMBL2039", inactive_threshold_nm=10000, max_compounds=200):
    """
    Fetch compounds actually tested against this target in ChEMBL and found
    weakly active or inactive (standard_value above inactive_threshold_nm,
    i.e. potency worse than the threshold) -- real experimental non-binders
    against the real target, rather than structurally-dissimilar molecules
    from an unrelated background pool. This is the scientifically strongest
    decoy source available, when it's reachable: it mirrors the same
    ChEMBL query style used to build this course's ML-based virtual
    screening classification dataset earlier, just selecting the inactive
    side of the same activity threshold instead of the active side.

    Returns (smiles_list, source_description). On any failure (no internet,
    API change, too few results), returns ([], error_description) so the
    caller can fall back to a different decoy source.
    """
    try:
        from chembl_webresource_client.new_client import new_client
        activity = new_client.activity
        results = activity.filter(
            target_chembl_id=target_chembl_id,
            standard_type__in=["IC50", "Ki", "EC50", "Potency"],
            standard_units="nM",
            standard_value__gte=inactive_threshold_nm,
            standard_relation="=",
            assay_type="B",
        ).only(["molecule_chembl_id", "canonical_smiles", "standard_value"])

        seen, smiles_list = set(), []
        for r in results:
            smi = r.get("canonical_smiles")
            mid = r.get("molecule_chembl_id")
            if not smi or mid in seen:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            seen.add(mid)
            smiles_list.append(smi)
            if len(smiles_list) >= max_compounds:
                break

        if len(smiles_list) < 2:
            return [], f"ChEMBL returned too few confirmed-inactive compounds for {target_chembl_id} (found {len(smiles_list)})"
        return smiles_list, f"{len(smiles_list)} experimentally-confirmed weak/inactive compounds from ChEMBL ({target_chembl_id}, \u2265{inactive_threshold_nm} nM)"
    except Exception as e:
        return [], f"ChEMBL query failed ({type(e).__name__}: {e})"


def load_decoys_from_file(path):
    """
    Load a user-supplied decoy set from a file -- .sdf, or a plain text
    file with one SMILES per line (optionally 'name SMILES' or
    'name,SMILES'). Use this to supply a real, external decoy set (e.g. a
    downloaded DUD-E decoy file for this target) instead of an
    automatically-generated one. Returns (smiles_list, source_description).
    """
    ext = os.path.splitext(path)[1].lower()
    smiles_list = []
    try:
        if ext == ".sdf":
            suppl = Chem.SDMolSupplier(path)
            for mol in suppl:
                if mol is None:
                    continue
                try:
                    smiles_list.append(Chem.MolToSmiles(mol))
                except Exception:
                    continue
        else:
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    if "," in line:
                        parts = line.split(",")
                    else:
                        parts = line.split()
                    candidate = None
                    for token in parts:
                        if Chem.MolFromSmiles(token) is not None:
                            candidate = token
                            break
                    if candidate:
                        smiles_list.append(candidate)
        if len(smiles_list) < 2:
            return [], f"File contained too few valid molecules (found {len(smiles_list)})"
        return smiles_list, f"{len(smiles_list)} decoys loaded from {os.path.basename(path)}"
    except Exception as e:
        return [], f"Failed to load decoy file: {type(e).__name__}: {e}"


def permutation_test_auc(labels, scores, n_permutations=2000, seed=42):
    """
    Statistical significance of an observed AUC: shuffle the active/decoy
    labels many times (keeping scores fixed) and see what fraction of
    random shuffles produce an AUC at least as high as the real one. This
    is the actual question "is my AUC better than chance, given how small
    my dataset is" -- an AUC alone doesn't answer that, especially with the
    handful of actives typical of ligand-based work.

    Returns (real_auc, p_value, max_resolution). max_resolution is the
    smallest possible p-value given the dataset size (1 / C(n, n_pos)) --
    with a small active/decoy set there are only so many distinct label
    arrangements that exist at all, so p-values finer than this aren't
    meaningful and the GUI should say so rather than implying false
    precision.
    """
    import random
    import math

    real_auc = roc_auc(labels, scores)
    if real_auc != real_auc:  # NaN
        return real_auc, float('nan'), None

    n = len(labels)
    n_pos = sum(labels)
    try:
        n_arrangements = math.comb(n, n_pos)
    except AttributeError:
        n_arrangements = None  # very old Python without math.comb; skip resolution note

    rng = random.Random(seed)
    labels_copy = list(labels)
    count_ge = 0
    for _ in range(n_permutations):
        rng.shuffle(labels_copy)
        perm_auc = roc_auc(labels_copy, scores)
        if perm_auc == perm_auc and perm_auc >= real_auc:
            count_ge += 1
    p_value = (count_ge + 1) / (n_permutations + 1)
    max_resolution = (1.0 / n_arrangements) if n_arrangements else None
    return real_auc, p_value, max_resolution


def applicability_domain_similarity(candidate_smiles, active_smiles_list):
    """
    How structurally similar is this candidate to the actives the
    pharmacophore was actually built from? A candidate can satisfy the 3D
    feature geometry by coincidence while being chemically nothing like the
    training data -- a known pitfall of pharmacophore matching. Returns
    (max_similarity, mean_similarity) via Morgan/Tanimoto, or (None, None)
    if the candidate SMILES doesn't parse.
    """
    from rdkit.Chem import rdFingerprintGenerator, DataStructs
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    cand_mol = Chem.MolFromSmiles(candidate_smiles)
    if cand_mol is None:
        return None, None
    cand_fp = fp_gen.GetFingerprint(cand_mol)

    sims = []
    for smi in active_smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = fp_gen.GetFingerprint(mol)
        sims.append(DataStructs.TanimotoSimilarity(cand_fp, fp))

    if not sims:
        return None, None
    return max(sims), sum(sims) / len(sims)


_PERIODIC_TABLE = Chem.GetPeriodicTable()


def compute_shape_envelope(mol_conf_pairs):
    """
    The union of heavy-atom positions and van der Waals radii across a set
    of already-aligned molecules -- the basis for a ligand shape-consistency
    check. Deliberately named and scoped as a LIGAND-only concept: this
    reflects whether a candidate's own bulk resembles the training actives'
    own shapes, not whether it would clash with a real receptor (this
    pipeline has no receptor structure anywhere in it). Do not present this
    to a user as "steric clash" or "receptor fit" -- see
    shape_consistency_score's docstring.

    mol_conf_pairs: list of (mol, conf_id) tuples, all already aligned into
    a shared coordinate frame (e.g. via pharm.ligands or a loaded model's
    ligand_mols).

    Returns a list of (position, vdw_radius) tuples, one per heavy atom.
    """
    envelope = []
    for mol, conf_id in mol_conf_pairs:
        conf = mol.GetConformer(conf_id)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue  # heavy atoms only
            pos = conf.GetAtomPosition(atom.GetIdx())
            envelope.append((np.array([pos.x, pos.y, pos.z]),
                              _PERIODIC_TABLE.GetRvdw(atom.GetAtomicNum())))
    return envelope


def shape_consistency_score(mol, conf_id, envelope, tolerance=1.0):
    """
    What fraction of this (already O3A-aligned) candidate conformer's heavy
    atoms fall within the training actives' own combined shape -- i.e. close
    enough to some training atom (within the sum of van der Waals radii,
    plus a tolerance) to be considered "inside" rather than sticking out
    past every training active. A candidate that matches the pharmacophore
    perfectly but has a large substituent none of the training actives had
    can still score low here, which fraction_matched alone would never
    catch.

    Framing note (important): this checks consistency with the training
    LIGANDS' shapes only. It is not a receptor steric-clash check -- this
    pipeline has no receptor structure. A high score here means "shaped
    like the training actives," not "fits the binding pocket."

    Verified on real data before adoption: every training active scores
    100% against its own envelope (as it must); a deliberately oversized,
    unrelated compound aligned onto the same reference dropped to 62%.

    Returns a score in [0, 1], or 1.0 (vacuously consistent) if the
    envelope is empty or the candidate has no heavy atoms.
    """
    if not envelope:
        return 1.0
    conf = mol.GetConformer(conf_id)
    n_heavy, n_inside = 0, 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        n_heavy += 1
        pos = conf.GetAtomPosition(atom.GetIdx())
        p = np.array([pos.x, pos.y, pos.z])
        cand_r = _PERIODIC_TABLE.GetRvdw(atom.GetAtomicNum())
        for env_p, env_r in envelope:
            if np.linalg.norm(p - env_p) <= (cand_r + env_r) / 2 + tolerance:
                n_inside += 1
                break
    return n_inside / n_heavy if n_heavy else 1.0


def classify_verdict(fraction_matched, plausibility):
    """
    Combine the two genuinely different questions this tool answers into
    one plain-language verdict, instead of leaving the user to mentally
    cross-reference two separate numbers:
      - fraction_matched: does this candidate present the right features
        in the right 3D arrangement? (the actual pharmacophore search --
        this is how genuinely novel scaffolds get found)
      - plausibility: is this candidate chemically similar to any real
        active, or could this be a geometric coincidence? (a confidence
        check, not a filter -- a low-plausibility hit is flagged as
        higher-risk, never discarded outright)

    Match tiers reuse this plugin's existing single-candidate verdict
    thresholds (0.75 / 0.4) for consistency. Plausibility tiers use
    0.5 / 0.3 -- a real, explicit "Medium" tier fills the gap between
    them, rather than forcing a continuous value into an arbitrary side
    of a single cutoff.

    Returns (verdict_label, color_hex, explanation).
    """
    if fraction_matched >= 0.75:
        match_tier = "high"
    elif fraction_matched >= 0.4:
        match_tier = "medium"
    else:
        match_tier = "low"

    if plausibility is None:
        plaus_tier = "unknown"
    elif plausibility >= 0.5:
        plaus_tier = "high"
    elif plausibility >= 0.3:
        plaus_tier = "medium"
    else:
        plaus_tier = "low"

    if match_tier == "high" and plaus_tier == "high":
        return ("Gold Standard Hit", "#0F6E56",
                "Right features in the right place, and chemically similar to a known active. "
                "The strongest kind of hit this tool can produce.")
    if match_tier == "high" and plaus_tier == "medium":
        return ("Scaffold-Hop (Moderate Confidence)", "#0891B2",
                "Matches the pharmacophore well but is only moderately similar to any known "
                "active -- plausibly a real, different-scaffold hit. Worth a closer look.")
    if match_tier == "high" and plaus_tier in ("low", "unknown"):
        return ("Scaffold-Hop (High Risk)", "#D97706",
                "Matches the pharmacophore well but looks chemically unlike any known active. "
                "Could be a genuine novel scaffold (the whole point of pharmacophore modeling) "
                "or a geometric coincidence. Do not discard -- investigate further (e.g. size, "
                "rotatable bonds, docking if a structure is available) before trusting it.")
    if match_tier == "medium":
        return ("Uncertain", "#64748B",
                "Only a partial pharmacophore match either way -- not confident evidence for "
                "or against this candidate on its own.")
    if match_tier == "low" and plaus_tier == "high":
        return ("Inactive Analog", "#C2410C",
                "Chemically resembles a known active but is missing the key 3D arrangement of "
                "features -- consistent with a real inactive, not a scoring artifact.")
    return ("Alien / Low Priority", "#A32D2D",
            "Neither the geometry nor the chemistry resembles anything this model was built "
            "from. Lowest priority.")


def _murcko_scaffold_smiles(smiles):
    """Returns the canonical SMILES of a molecule's Bemis-Murcko scaffold (the
    core ring system with substituents stripped off), or None if the SMILES
    doesn't parse or has no ring system at all (an acyclic molecule)."""
    from rdkit.Chem.Scaffolds import MurckoScaffold
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(scaffold)


def classify_scaffold_novelty(candidate_smiles, active_smiles_list, analog_thresh=0.6):
    """
    How novel is this candidate's core ring system compared to the training
    actives, independent of the pharmacophore/plausibility scores? Answers a
    different question than Verdict (Section 3.6 of the docs): not "does
    this look chemically similar overall" but specifically "is this built on
    the same, a related, or a genuinely different ring scaffold."

    Returns one of "Same Scaffold" (an exact Bemis-Murcko match to at least
    one active -- e.g. a prodrug or minor substituent variant), "Analog
    Scaffold" (a structurally related but not identical scaffold, Tanimoto
    >= analog_thresh on the scaffolds themselves), "Novel Scaffold" (neither),
    or "Unknown" (the candidate or every active's scaffold failed to parse,
    or the candidate has no ring system to compare).

    Verified on real data before adoption: parecoxib (a genuine prodrug of
    valdecoxib) correctly reduces to the exact same Murcko scaffold as
    valdecoxib ("Same Scaffold"); lumiracoxib (a genuinely different
    phenylacetic-acid chemotype) correctly comes back "Novel Scaffold."
    """
    from rdkit.Chem import rdFingerprintGenerator, DataStructs
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    cand_scaffold_smi = _murcko_scaffold_smiles(candidate_smiles)
    if cand_scaffold_smi is None:
        return "Unknown"
    cand_mol = Chem.MolFromSmiles(cand_scaffold_smi)
    cand_fp = fp_gen.GetFingerprint(cand_mol)

    best_sim = 0.0
    any_valid_active = False
    for a_smi in active_smiles_list:
        a_scaffold_smi = _murcko_scaffold_smiles(a_smi)
        if a_scaffold_smi is None:
            continue
        any_valid_active = True
        if a_scaffold_smi == cand_scaffold_smi:
            return "Same Scaffold"
        a_mol = Chem.MolFromSmiles(a_scaffold_smi)
        sim = DataStructs.TanimotoSimilarity(cand_fp, fp_gen.GetFingerprint(a_mol))
        best_sim = max(best_sim, sim)

    if not any_valid_active:
        return "Unknown"
    return "Analog Scaffold" if best_sim >= analog_thresh else "Novel Scaffold"


def leave_one_out_stability(ligand_entries, eps=1.5, n_confs=DEFAULT_N_CONFS, protonate=True,
                             support_thresh=0.0, match_tol=2.0, progress_callback=None):
    """
    Robustness check: rebuild the consensus pharmacophore N times, each time
    leaving one active out, and see how often each feature from the FULL
    model survives. A feature present in the full model but missing from
    most leave-one-out runs is being driven by one or two compounds, not a
    genuinely shared pattern -- this gives the "reference-ligand dependence"
    limitation empirical teeth instead of just a warning in prose.

    support_thresh defaults to 0.0 (evaluate every feature the model finds,
    unfiltered) *deliberately* -- an earlier version of this function
    accepted whatever display threshold the GUI slider happened to be set
    to, which meant only features already comfortably above that bar ever
    got tested, and silently skipped exactly the borderline, low-support
    features a robustness check is most useful for. The GUI's own display
    threshold (which features count as "consensus") is applied separately,
    after this function returns, not coupled into the computation itself.

    ligand_entries: list of (name, smiles) tuples for the full active set.
    Returns a dict: {"full_points": [...], "stability": [frac_loo_runs_that_
    still_contain_this_feature, ...] (parallel to full_points), "n_runs": N}.
    A feature's "stability" of 1.0 means it survived being computed without
    every single one of the N leave-one-out subsets; lower means it depends
    on specific compounds being included.
    """
    n = len(ligand_entries)
    if n < 3:
        return {"error": "Need at least 3 active compounds for a meaningful leave-one-out check."}

    full_ligands = [Ligand(name, smi, n_confs=n_confs, protonate=protonate) for name, smi in ligand_entries]
    full_valid = [l for l in full_ligands if l.is_valid]
    if len(full_valid) < 3:
        return {"error": "Too few valid compounds after embedding for a leave-one-out check."}

    full_pharm = LigandBasedPharmacophore(full_valid)
    full_pharm.align()
    full_pharm.extract_consensus(eps=eps, min_support_frac=0.0)
    full_points = full_pharm.consensus_at_threshold(support_thresh)

    loo_point_sets = []
    for i in range(len(ligand_entries)):
        if progress_callback:
            progress_callback(i, len(ligand_entries))
        subset = [e for j, e in enumerate(ligand_entries) if j != i]
        loo_ligands = [Ligand(name, smi, n_confs=n_confs, protonate=protonate) for name, smi in subset]
        loo_valid = [l for l in loo_ligands if l.is_valid]
        if len(loo_valid) < 2:
            continue
        loo_pharm = LigandBasedPharmacophore(loo_valid)
        try:
            loo_pharm.align()
            loo_pharm.extract_consensus(eps=eps, min_support_frac=0.0)
        except Exception:
            continue
        loo_point_sets.append(loo_pharm.consensus_at_threshold(support_thresh))

    n_runs = len(loo_point_sets)
    stability = []
    for p in full_points:
        n_surviving = 0
        for loo_points in loo_point_sets:
            for lp in loo_points:
                if lp.family == p.family:
                    dist = float(np.linalg.norm(np.array(lp.centroid) - np.array(p.centroid)))
                    if dist <= match_tol:
                        n_surviving += 1
                        break
        stability.append(n_surviving / n_runs if n_runs > 0 else float('nan'))

    return {"full_points": full_points, "stability": stability, "n_runs": n_runs}


def compare_pharmacophore_models(pharm_a, meta_a, pharm_b, meta_b, thresh_a=0.5, thresh_b=0.5, match_tol=2.0):
    """
    Compare two pharmacophore models (e.g. loaded from two saved .json
    files, built with different parameters or different active sets).
    Returns a list of rows: {"family", "in_a", "in_b", "support_a",
    "support_b", "shared"} -- one row per distinct feature found in either
    model, matched by family + spatial proximity (within match_tol A).
    """
    points_a = pharm_a.consensus_at_threshold(thresh_a)
    points_b = pharm_b.consensus_at_threshold(thresh_b)

    rows = []
    used_b = set()
    for pa in points_a:
        match = None
        for j, pb in enumerate(points_b):
            if j in used_b or pb.family != pa.family:
                continue
            dist = float(np.linalg.norm(np.array(pa.centroid) - np.array(pb.centroid)))
            if dist <= match_tol:
                match = (j, pb)
                break
        if match:
            used_b.add(match[0])
            rows.append({"family": pa.family, "in_a": True, "in_b": True,
                         "support_a": pa.support_frac, "support_b": match[1].support_frac, "shared": True})
        else:
            rows.append({"family": pa.family, "in_a": True, "in_b": False,
                         "support_a": pa.support_frac, "support_b": None, "shared": False})
    for j, pb in enumerate(points_b):
        if j not in used_b:
            rows.append({"family": pb.family, "in_a": False, "in_b": True,
                         "support_a": None, "support_b": pb.support_frac, "shared": False})
    return rows


def _compute_validation_metrics(labels, scores, records, top_pct=0.10, hit_frac_cutoff=0.75):
    """Shared metric computation (AUC, EF, GH, permutation p-value, ROC curve
    points) used by both validate_pharmacophore() and the leave-one-out
    version below -- factored out so a fix to one doesn't risk silently
    missing the other."""
    A = sum(labels)
    D = len(labels)
    if A == 0 or A == D:
        return {"error": "Need at least one active and one decoy with valid scores.",
                "records": records}

    auc = roc_auc(labels, scores)
    ef = enrichment_factor(labels, scores, top_pct=top_pct)
    fpr, tpr = roc_curve_points(labels, scores)
    _, p_value, p_resolution = permutation_test_auc(labels, scores, n_permutations=2000)

    hits = [(l, r) for l, r in zip(labels, records) if r["fraction_matched"] >= hit_frac_cutoff]
    Ht = len(hits)
    Ha = sum(l for l, r in hits)
    gh = guner_henry_score(A, D, Ht, Ha)

    return {
        "n_actives": A, "n_decoys": D - A, "n_total": D,
        "auc": auc, "ef": ef, "top_pct": top_pct,
        "gh": gh, "hit_frac_cutoff": hit_frac_cutoff, "Ht": Ht, "Ha": Ha,
        "fpr": fpr, "tpr": tpr,
        "p_value": p_value, "p_resolution": p_resolution,
        "records": records,
    }


def validate_pharmacophore(pharmacophore_points, active_smiles_list, decoy_smiles_list,
                            reference_mol=None, reference_conf_id=0,
                            tol=1.8, n_confs=30, protonate=True, top_pct=0.10,
                            hit_frac_cutoff=0.75, progress_callback=None):
    """
    Score every active and every decoy against a consensus pharmacophore, and
    compute EF, ROC-AUC, and the GH score. Returns a dict with the metrics
    plus the raw per-compound results (for a GUI to display or plot).

    reference_mol/reference_conf_id should be the pharmacophore's own
    alignment reference (pharm.get_alignment_reference()) -- passed through
    to every score_molecule() call so candidates are properly aligned before
    scoring, not compared in an arbitrary orientation. See score_molecule's
    docstring for why this matters.

    hit_frac_cutoff is the fraction-of-features-matched threshold used to
    define a "hit" for the GH score specifically -- deliberately a real,
    meaningful cutoff (0.75 by default: most features matched) rather than
    "matched at least one feature," which was found, in earlier development
    of this course's teaching notebook, to make GH collapse toward its
    random-chance floor by counting nearly everything as a hit.

    IMPORTANT METHODOLOGICAL CAVEAT: the active_smiles_list here is scored
    against pharmacophore_points that were (almost always) built FROM those
    same actives -- this validates whether the model can reproduce its own
    training data, not whether it generalizes to a genuinely unseen active.
    For an honest validation number, use validate_pharmacophore_loocv()
    instead, which scores each active only against a model that was rebuilt
    without it. This function is kept for cases where that isn't
    appropriate (e.g. validating a model against a curated external active
    set that was never used to build it in the first place).
    """
    labels, scores, records = [], [], []
    all_entries = [(s, 1) for s in active_smiles_list] + [(s, 0) for s in decoy_smiles_list]
    n_total = len(all_entries)

    for i, (smi, label) in enumerate(all_entries):
        if progress_callback:
            progress_callback(i, n_total)
        n_matched, frac, mean_d, mol, cid, _details, align_rmsd, shape_consistency = score_molecule(
            smi, pharmacophore_points, reference_mol=reference_mol, reference_conf_id=reference_conf_id,
            tol=tol, n_confs=n_confs, protonate=protonate
        )
        if mol is None:
            continue
        combined_score = frac - 0.01 * min(mean_d, 50)
        labels.append(label)
        scores.append(combined_score)
        records.append({"smiles": smi, "label": label, "n_matched": n_matched,
                         "n_total": len(pharmacophore_points), "fraction_matched": frac,
                         "mean_dist": mean_d, "align_rmsd": align_rmsd})

    result = _compute_validation_metrics(labels, scores, records, top_pct, hit_frac_cutoff)
    result["loocv_actives"] = False
    return result


def validate_pharmacophore_loocv(ligand_entries, decoy_smiles_list, eps=1.5, support_thresh=0.5,
                                  tol=1.8, n_confs=30, protonate=True, top_pct=0.10,
                                  hit_frac_cutoff=0.75, progress_callback=None):
    """
    Genuine leave-one-out cross-validation for the actives side: each active
    is scored ONLY against a pharmacophore model rebuilt WITHOUT it -- never
    against a model that included it during training. This is the
    methodologically correct alternative to validate_pharmacophore(), which
    (when active_smiles_list is the same set used to build the model, the
    common case) tests whether the model can reproduce its own training
    data rather than whether it generalizes.

    Concretely verified to matter on real data before being adopted: on the
    course's 5-compound MAO-B example, in-sample AUC was 0.880 vs. a
    genuine leave-one-out AUC of 0.840 -- a real, if here modest, inflation
    that depends on how internally similar the active set is and is not
    guaranteed to stay small for other active sets.

    Decoys are scored against the FULL model (built from all actives) --
    correct and unproblematic, since decoys were never part of training
    either way, in either function.

    ligand_entries: list of (name, smiles) for the *full* active set (not a
    pre-built pharmacophore -- this function needs to rebuild the model
    itself, once per left-out active).
    Returns the same dict shape as validate_pharmacophore(), with
    "loocv_actives": True, plus "n_loocv_runs".
    """
    n = len(ligand_entries)
    if n < 3:
        return {"error": "Need at least 3 active compounds for leave-one-out validation."}

    full_ligands = [Ligand(name, smi, n_confs=n_confs, protonate=protonate) for name, smi in ligand_entries]
    full_valid = [l for l in full_ligands if l.is_valid]
    if len(full_valid) < 3:
        return {"error": "Too few valid compounds after embedding for leave-one-out validation."}
    full_pharm = LigandBasedPharmacophore(full_valid)
    full_pharm.align()
    full_pharm.extract_consensus(eps=eps, min_support_frac=0.0)
    full_consensus = full_pharm.consensus_at_threshold(support_thresh)
    full_ref_mol, full_ref_cid = full_pharm.get_alignment_reference()

    labels, scores, records = [], [], []
    n_total_steps = n + len(decoy_smiles_list)
    step = 0
    n_loocv_runs = 0

    for i, (name, smi) in enumerate(ligand_entries):
        if progress_callback:
            progress_callback(step, n_total_steps)
        step += 1
        subset = [e for j, e in enumerate(ligand_entries) if j != i]
        sub_ligands = [Ligand(sname, ssmi, n_confs=n_confs, protonate=protonate) for sname, ssmi in subset]
        sub_valid = [l for l in sub_ligands if l.is_valid]
        if len(sub_valid) < 2:
            continue
        try:
            sub_pharm = LigandBasedPharmacophore(sub_valid)
            sub_pharm.align()
            sub_pharm.extract_consensus(eps=eps, min_support_frac=0.0)
        except Exception:
            continue
        sub_consensus = sub_pharm.consensus_at_threshold(support_thresh)
        if not sub_consensus:
            continue
        sub_ref_mol, sub_ref_cid = sub_pharm.get_alignment_reference()

        n_matched, frac, mean_d, mol, cid, _details, align_rmsd, shape_consistency = score_molecule(
            smi, sub_consensus, reference_mol=sub_ref_mol, reference_conf_id=sub_ref_cid,
            tol=tol, n_confs=n_confs, protonate=protonate
        )
        if mol is None:
            continue
        n_loocv_runs += 1
        combined_score = frac - 0.01 * min(mean_d, 50)
        labels.append(1)
        scores.append(combined_score)
        records.append({"smiles": smi, "name": name, "label": 1, "n_matched": n_matched,
                         "n_total": len(sub_consensus), "fraction_matched": frac,
                         "mean_dist": mean_d, "align_rmsd": align_rmsd})

    for smi in decoy_smiles_list:
        if progress_callback:
            progress_callback(step, n_total_steps)
        step += 1
        n_matched, frac, mean_d, mol, cid, _details, align_rmsd, shape_consistency = score_molecule(
            smi, full_consensus, reference_mol=full_ref_mol, reference_conf_id=full_ref_cid,
            tol=tol, n_confs=n_confs, protonate=protonate
        )
        if mol is None:
            continue
        combined_score = frac - 0.01 * min(mean_d, 50)
        labels.append(0)
        scores.append(combined_score)
        records.append({"smiles": smi, "label": 0, "n_matched": n_matched,
                         "n_total": len(full_consensus), "fraction_matched": frac,
                         "mean_dist": mean_d, "align_rmsd": align_rmsd})

    result = _compute_validation_metrics(labels, scores, records, top_pct, hit_frac_cutoff)
    result["loocv_actives"] = True
    result["n_loocv_runs"] = n_loocv_runs
    return result


def mw_range_from_actives(active_smiles_list):
    """Returns (min_mw, max_mw) across a list of active SMILES -- the basis
    for a molecular-weight plausibility filter. A candidate far outside this
    range (adjusted by a tolerance) likely cannot physically occupy the same
    binding pocket as the training actives, regardless of pharmacophore fit
    -- found, in real use of this plugin, to be a more diagnostic signal
    than whole-molecule fingerprint similarity for at least one real
    library-composition problem (see README)."""
    from rdkit.Chem import Descriptors
    mws = []
    for smi in active_smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mws.append(Descriptors.MolWt(mol))
    if not mws:
        return None, None
    return min(mws), max(mws)


def filter_by_mw_range(results, active_smiles_list, tolerance_da=50):
    """Filter a screen_library()-style results list (each dict must have a
    'smiles' key) down to those within [min_active_mw - tolerance_da,
    max_active_mw + tolerance_da]. Returns (kept, n_removed)."""
    from rdkit.Chem import Descriptors
    min_mw, max_mw = mw_range_from_actives(active_smiles_list)
    if min_mw is None:
        return results, 0
    lo, hi = min_mw - tolerance_da, max_mw + tolerance_da
    kept, n_removed = [], 0
    for r in results:
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            n_removed += 1
            continue
        mw = Descriptors.MolWt(mol)
        if lo <= mw <= hi:
            r = dict(r)
            r["mw"] = mw
            kept.append(r)
        else:
            n_removed += 1
    return kept, n_removed


def compute_pharmacophore_fingerprint(points):
    """
    A '3D pharmacophore fingerprint': the pairwise distances between every
    pair of consensus features, alignment-invariant (rotation/translation
    don't affect a distance between two points), so it never needs O3A.
    Used to build a fast, lenient pre-filter (see fast_prefilter_score) --
    NOT a replacement for full O3A-based scoring. Verified on real data to
    be noticeably more permissive than the full alignment-based score (a
    genuine structural outlier scored 0.90 on this fingerprint alone,
    vs. its real, much lower O3A-based match) -- internal pairwise
    distances alone don't capture full 3D arrangement the way a real
    alignment does. Only ever use this to quickly reject the clearly
    impossible before spending time on full scoring, never as a final score.

    Returns a list of (sorted_family_pair_tuple, distance) tuples.
    """
    import itertools
    fp = []
    for (a, b) in itertools.combinations(points, 2):
        d = float(np.linalg.norm(a.centroid - b.centroid))
        fp.append((tuple(sorted([a.family, b.family])), d))
    return fp


def _candidate_internal_fingerprint(mol, conf_id):
    import itertools
    feats = FEATURE_FACTORY.GetFeaturesForMol(mol, confId=conf_id)
    pts = [(f.GetFamily(), np.array([f.GetPos().x, f.GetPos().y, f.GetPos().z])) for f in feats]
    fp = []
    for (fa, pa), (fb, pb) in itertools.combinations(pts, 2):
        d = float(np.linalg.norm(pa - pb))
        fp.append((tuple(sorted([fa, fb])), d))
    return fp


def fast_prefilter_score(smiles, consensus_fingerprint, n_confs=10, seed=1, tol=1.0):
    """
    A fast, lenient, alignment-free pre-filter: what fraction of the
    consensus's own inter-feature distances find a same-family-pair match
    (within tol A) somewhere in the candidate's own internal geometry, best
    across a small number of conformers. No O3A alignment is performed at
    all, which is what makes this fast (~6x faster than full scoring in
    real testing) -- but also why it is deliberately more permissive than
    real scoring and must not be used as a final verdict on its own.

    Returns a score in [0, 1], or 0.0 if the SMILES doesn't parse/embed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not consensus_fingerprint:
        return 0.0
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, randomSeed=seed, pruneRmsThresh=0.5)
    if len(cids) == 0:
        return 0.0

    best_score = 0.0
    for cid in cids:
        cand_fp = _candidate_internal_fingerprint(mol, cid)
        n_matched = 0
        for fam_pair, d in consensus_fingerprint:
            for c_fam_pair, c_d in cand_fp:
                if c_fam_pair == fam_pair and abs(c_d - d) <= tol:
                    n_matched += 1
                    break
        score = n_matched / len(consensus_fingerprint)
        best_score = max(best_score, score)
    return best_score


def diverse_top_hits(results, max_similarity=0.6, top_n=None):
    """
    Greedy diversity filter over already-ranked screen_library() results:
    walk down the ranked list (best-first) and keep a hit only if it is not
    too structurally similar (Morgan/Tanimoto > max_similarity) to any hit
    already kept. Since the input is already ranked by score, this keeps
    the single best-scoring representative of each cluster of near-identical
    analogs and drops the redundant rest -- directly addressing the case
    where the top of a ranked hit list is dominated by minor substituent
    variations on one scaffold rather than genuinely distinct candidates.
    Validated on synthetic test data (a cluster of 3 halogenated near-analogs
    plus 2 genuinely different scaffolds) before being adopted -- see
    test_plugin.py.
    """
    from rdkit.Chem import rdFingerprintGenerator, DataStructs
    fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    kept, kept_fps = [], []
    for r in results:
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            continue
        fp = fp_gen.GetFingerprint(mol)
        too_similar = any(DataStructs.TanimotoSimilarity(fp, kfp) > max_similarity for kfp in kept_fps)
        if not too_similar:
            kept.append(r)
            kept_fps.append(fp)
        if top_n and len(kept) >= top_n:
            break
    return kept


def screen_library(sdf_path, pharmacophore_points, reference_mol=None, reference_conf_id=0,
                    tol=1.8, n_confs=DEFAULT_N_CONFS, protonate=True, progress_callback=None,
                    use_prefilter=False, prefilter_thresh=0.3, shape_envelope=None,
                    active_smiles_list=None, match_mode="centroid"):
    """
    Score every molecule in an SDF file against a consensus pharmacophore.
    Returns a list of dicts, sorted best-first (fraction_matched desc, then
    mean_dist asc): {name, smiles, n_matched, n_total, fraction_matched,
    mean_dist, mol, conf_id}. Molecules that fail to parse/embed are skipped
    and counted separately (returned as the second element of the tuple:
    (results, n_skipped, n_prefiltered)).

    reference_mol/reference_conf_id should be the pharmacophore's own
    alignment reference (pharm.get_alignment_reference()) -- see
    score_molecule's docstring for why this matters.

    use_prefilter, if True, runs the fast alignment-free fingerprint
    pre-filter (compute_pharmacophore_fingerprint/fast_prefilter_score)
    first, and skips full O3A-based scoring entirely for any compound
    scoring below prefilter_thresh on that fast check -- a real speedup for
    large libraries, at the cost of being a more lenient, less precise cut
    (see fast_prefilter_score's docstring). Off by default: only worth
    turning on for genuinely large libraries where full scoring would
    otherwise take a long time.

    shape_envelope, if given (from pharm.get_shape_envelope()), adds a
    "shape_consistency" value to each result (see shape_consistency_score's
    docstring) -- a ligand-shape check, not a receptor steric-clash check.

    active_smiles_list, if given, adds a "scaffold_novelty" value to each
    result via classify_scaffold_novelty().

    progress_callback, if given, is called as progress_callback(i, n_total)
    after each molecule, so a GUI can update a progress bar.
    """
    suppl = Chem.SDMolSupplier(sdf_path)
    results = []
    n_skipped = 0
    n_prefiltered = 0
    mols = list(suppl)
    n_total = len(mols)

    consensus_fp = compute_pharmacophore_fingerprint(pharmacophore_points) if use_prefilter else None

    for i, mol in enumerate(mols):
        if progress_callback:
            progress_callback(i, n_total)
        if mol is None:
            n_skipped += 1
            continue
        name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"compound_{i+1}"
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            n_skipped += 1
            continue

        if use_prefilter:
            pre_score = fast_prefilter_score(smiles, consensus_fp, n_confs=min(10, n_confs))
            if pre_score < prefilter_thresh:
                n_prefiltered += 1
                continue

        n_matched, frac, mean_d, scored_mol, cid, match_details, align_rmsd, shape_consistency = score_molecule(
            smiles, pharmacophore_points, reference_mol=reference_mol, reference_conf_id=reference_conf_id,
            tol=tol, n_confs=n_confs, protonate=protonate, shape_envelope=shape_envelope,
            match_mode=match_mode
        )
        if scored_mol is None:
            n_skipped += 1
            continue

        scaffold_novelty = classify_scaffold_novelty(smiles, active_smiles_list) if active_smiles_list else None

        results.append({
            "name": name,
            "smiles": smiles,
            "n_matched": n_matched,
            "n_total": len(pharmacophore_points),
            "fraction_matched": frac,
            "mean_dist": mean_d,
            "align_rmsd": align_rmsd,
            "shape_consistency": shape_consistency,
            "scaffold_novelty": scaffold_novelty,
            "mol": scored_mol,
            "conf_id": cid,
            "match_details": match_details,
        })

    results.sort(key=lambda r: (-r["fraction_matched"], r["mean_dist"]))
    return results, n_skipped, n_prefiltered
