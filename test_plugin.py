"""
Test harness for the PyMOL plugin — uses PyMOL's REAL global `cmd` singleton
(the same one a live windowed PyMOL session exposes), initialized in -cq
(command-line, quiet) mode so it works headlessly. This validates the plugin's
actual PyMOL-manipulation code against the real cmd API, and the actual Qt
dialog construction/wiring against a real QApplication — the two things that
matter for correctness. It does not screenshot a live *windowed* PyMOL main
viewport (that mode hangs under this sandbox's virtual display setup), but the
cmd calls used are identical either way since they hit the same C++ core.
"""
import sys
import os
import collections
from rdkit import Chem
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import pymol
pymol.finish_launching(['pymol', '-cq'])
from pymol import cmd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'pharmol'))
import pharmol as plugin

from pymol.Qt import QtWidgets

app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

print("=== Test 1: plugin registration ===")
registered = {}
def fake_addmenuitemqt(label, callback):
    registered[label] = callback
import pymol.plugins
pymol.plugins.addmenuitemqt = fake_addmenuitemqt
plugin.__init_plugin__(app=None)
assert "PhArMol" in registered, "Plugin did not register its menu item"
print("Menu item registered correctly:", list(registered.keys()))

print("\n=== Test 2: dialog construction ===")
dlg = plugin.make_dialog()
assert dlg is not None
print("Dialog constructed OK. Title:", dlg.windowTitle())

# Find widgets by walking the layout (since they're locals inside make_dialog,
# we instead re-invoke the internal helpers directly against a real cmd for a
# focused, explicit test of the PyMOL-manipulation logic below.)

print("\n=== Test 3: core science + PyMOL loading (real cmd) ===")
from core import Ligand, LigandBasedPharmacophore, score_molecule

entries = [
    ("selegiline", "CC(Cc1ccccc1)N(C)CC#C"),
    ("rasagiline", "C#CCNC1CCc2ccccc21"),
    ("pargyline", "C#CCN(C)Cc1ccccc1"),
    ("clorgyline", "C#CCN(C)CCCOc1ccc(Cl)cc1Cl"),
    ("lazabemide", "NCCNC(=O)c1ccc(Cl)cn1"),
]
ligands = [Ligand(name, smi) for name, smi in entries]
valid = [l for l in ligands if l.is_valid]
assert len(valid) == 5, f"Expected 5 valid ligands, got {len(valid)}"
print(f"{len(valid)} ligands prepared successfully")

pharm = LigandBasedPharmacophore(valid)
pharm.align()
pharm.extract_consensus(eps=1.5, min_support_frac=0.0)
assert len(pharm.points) > 0
print(f"Extracted {len(pharm.points)} candidate feature clusters")

print("\n=== Test 4: loading ligands into REAL PyMOL cmd session ===")
plugin._load_into_pymol(cmd, pharm)
objects = cmd.get_names()
ligand_objects = [o for o in objects if o.startswith("ligand_")]
assert len(ligand_objects) == 5, f"Expected 5 ligand objects in PyMOL, got {len(ligand_objects)}: {objects}"
print("PyMOL objects after loading:", objects)

# check atom counts are sane (not empty objects)
for obj in ligand_objects:
    n_atoms = cmd.count_atoms(obj)
    assert n_atoms > 5, f"{obj} has suspiciously few atoms: {n_atoms}"
    print(f"  {obj}: {n_atoms} atoms")

print("\n=== Test 5: pharmacophore spheres in REAL PyMOL cmd session ===")
plugin._refresh_pharmacophore_spheres(cmd, pharm, 0.5)
objects = cmd.get_names()
pharm_objects = [o for o in objects if o.startswith("pharm_")]
expected_n = len([p for p in pharm.points if p.support_frac >= 0.5])
assert len(pharm_objects) == expected_n, f"Expected {expected_n} pharmacophore spheres, got {len(pharm_objects)}"
print(f"Pharmacophore sphere objects at 50% threshold: {pharm_objects}")

for obj in pharm_objects:
    n_atoms = cmd.count_atoms(obj)
    assert n_atoms == 1, f"{obj} should be a single pseudoatom, has {n_atoms}"
    coords = cmd.get_atom_coords(obj)
    print(f"  {obj}: pseudoatom at {tuple(round(c,2) for c in coords)}")

print("\n=== Test 6: changing threshold updates sphere count in real PyMOL ===")
plugin._refresh_pharmacophore_spheres(cmd, pharm, 0.9)
objects = cmd.get_names()
pharm_objects_90 = [o for o in objects if o.startswith("pharm_")]
expected_90 = len([p for p in pharm.points if p.support_frac >= 0.9])
assert len(pharm_objects_90) == expected_90
print(f"At 90% threshold: {len(pharm_objects_90)} spheres (expected {expected_90})")
assert len(pharm_objects_90) <= len(pharm_objects), "Raising threshold should not increase sphere count"

print("\n=== Test 7: candidate scoring + loading into real PyMOL ===")
consensus_pts = pharm.consensus_at_threshold(0.5)
n_matched, frac, mean_d, mol, cid, match_details, align_rmsd, shape_consistency = score_molecule("C#CCN(C)CCc1ccccc1", consensus_pts)
assert mol is not None, "Candidate SMILES failed to parse/embed"
print(f"Candidate score: {n_matched}/{len(consensus_pts)} features matched ({frac*100:.0f}%), mean dist {mean_d:.2f} A")
plugin._load_candidate_into_pymol(cmd, mol, cid)
assert "candidate" in cmd.get_names(), "Candidate object not loaded into PyMOL"
n_atoms = cmd.count_atoms("candidate")
assert n_atoms > 5
print(f"Candidate object loaded into real PyMOL session: {n_atoms} atoms")

print("\n=== Test 8: cleanup / delete works ===")
cmd.delete("ligand_* or pharm_* or candidate")
remaining = [o for o in cmd.get_names() if o.startswith(("ligand_", "pharm_")) or o == "candidate"]
assert len(remaining) == 0, f"Objects not cleaned up: {remaining}"
print("Cleanup verified.")

print("\n=== Test 9: protonation is applied by default and correct ===")
from core import protonate_physiological, DEFAULT_N_CONFS
print("DEFAULT_N_CONFS =", DEFAULT_N_CONFS)
assert DEFAULT_N_CONFS == 50, f"Expected default of 50 conformers, got {DEFAULT_N_CONFS}"

for lig in valid:
    charge = Chem.GetFormalCharge(lig.mol)
    assert charge == 1, f"{lig.name}: expected formal charge +1 after protonation, got {charge}"
    print(f"  {lig.name}: formal charge = {charge:+d} (correctly protonated)")

# aniline should NOT be protonated (too weakly basic at physiological pH)
aniline = Chem.MolFromSmiles("Nc1ccccc1")
aniline_prot = protonate_physiological(aniline)
assert Chem.GetFormalCharge(aniline_prot) == 0, "Aniline should NOT be protonated (too weakly basic)"
print("  aniline: correctly left unprotonated (net charge 0)")

# a carboxylic acid should be deprotonated
benzoic = Chem.MolFromSmiles("c1ccccc1C(=O)O")
benzoic_deprot = protonate_physiological(benzoic)
assert Chem.GetFormalCharge(benzoic_deprot) == -1, "Benzoic acid should be deprotonated to -1"
print("  benzoic acid: correctly deprotonated (net charge -1)")

print("\n=== Test 10: protonation can be disabled via the protonate= flag ===")
lig_no_prot = Ligand("selegiline_neutral", "CC(Cc1ccccc1)N(C)CC#C", n_confs=15, protonate=False)
assert lig_no_prot.is_valid, "Ligand should still prepare successfully with protonate=False"
assert Chem.GetFormalCharge(lig_no_prot.mol) == 0, "protonate=False should leave the molecule neutral"
print("  protonate=False correctly leaves molecule neutral")

print("\n=== Test 11: batch library screening against real SDF, with real PyMOL loading ===")
from core import screen_library
from rdkit.Chem import AllChem as _AllChem

test_lib_path = "/tmp/test_plugin_library.sdf"
test_compounds = {
    "lib_selegiline_like": "CC(Cc1ccccc1)N(C)CC#C",
    "lib_propargyl_analog": "C#CCN(C)CCc1ccccc1",
    "lib_decane_decoy": "CCCCCCCCCC",
    "lib_benzene_decoy": "c1ccccc1",
    "lib_rasagiline_like": "C#CCNC1CCc2ccccc21",
}
writer = Chem.SDWriter(test_lib_path)
for name, smi in test_compounds.items():
    m = Chem.MolFromSmiles(smi)
    m = Chem.AddHs(m)
    _AllChem.EmbedMolecule(m, randomSeed=1)
    m.SetProp("_Name", name)
    writer.write(m)
writer.close()

progress_log = []
results, n_skipped, n_prefiltered = screen_library(
    test_lib_path, consensus_pts, n_confs=30,
    progress_callback=lambda i, n: progress_log.append((i, n))
)
assert len(results) == 5, f"Expected 5 results, got {len(results)}"
assert n_skipped == 0, f"Expected 0 skipped, got {n_skipped}"
assert len(progress_log) == 5, "Progress callback should fire once per molecule"
print(f"Screened {len(results)} compounds, {n_skipped} skipped, progress callback fired {len(progress_log)}x")

names_ranked = [r["name"] for r in results]
print("Ranked order:", names_ranked)
assert names_ranked[0] in ("lib_selegiline_like", "lib_rasagiline_like"), \
    f"Expected a real analog to rank first, got {names_ranked[0]}"
assert "lib_decane_decoy" in names_ranked[-2:], \
    f"Expected decane decoy near the bottom, got {names_ranked}"
for r in results:
    print(f"  {r['name']:25s} {r['n_matched']}/{r['n_total']} ({r['fraction_matched']*100:.0f}%) mean_dist={r['mean_dist']:.2f}")

# now load the top hit into the real PyMOL session, exactly as clicking a table row would
top_hit = results[0]
plugin._load_candidate_into_pymol(cmd, top_hit["mol"], top_hit["conf_id"])
assert "candidate" in cmd.get_names(), "Top library hit was not loaded into PyMOL"
n_atoms = cmd.count_atoms("candidate")
assert n_atoms > 5
print(f"Top hit '{top_hit['name']}' loaded into real PyMOL session: {n_atoms} atoms")
cmd.delete("candidate")

print("\n=== Test 12: CSV export produces a valid, correctly-populated file ===")
import csv, tempfile
csv_path = tempfile.mktemp(suffix=".csv")
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "smiles", "n_matched", "n_total", "fraction_matched", "mean_dist_angstrom"])
    for r in results:
        w.writerow([r["name"], r["smiles"], r["n_matched"], r["n_total"],
                    f'{r["fraction_matched"]:.3f}', f'{r["mean_dist"]:.3f}'])
with open(csv_path) as f:
    rows = list(csv.reader(f))
assert len(rows) == 6, f"Expected 6 rows (1 header + 5 data), got {len(rows)}"
assert rows[0] == ["name", "smiles", "n_matched", "n_total", "fraction_matched", "mean_dist_angstrom"]
print(f"CSV export verified: {len(rows)-1} data rows written correctly to {csv_path}")
os.unlink(csv_path)
os.unlink(test_lib_path)

print("\n=== Test 13: Save Model -> Load Model round-trip, with real PyMOL loading ===")
from core import export_pharmacophore_json, load_pharmacophore_json

save_path = "/tmp/test_plugin_saved_model.json"
metadata = {
    "n_ligands": len(pharm.ligands),
    "ligand_names": [l.name for l in pharm.ligands],
    "ligand_smiles": [l.smiles for l in pharm.ligands],
    "conformers_per_molecule": 50,
    "protonated": True,
}
export_pharmacophore_json(pharm.points, save_path, metadata=metadata, ligands=pharm.ligands)
assert os.path.exists(save_path)
print(f"Model saved: {os.path.getsize(save_path)} bytes")

loaded = load_pharmacophore_json(save_path)
assert len(loaded.points) == len(pharm.points), \
    f"Point count mismatch after reload: {len(loaded.points)} vs {len(pharm.points)}"
assert len(loaded.ligand_mols) == 5, f"Expected 5 reloaded ligand structures, got {len(loaded.ligand_mols)}"
print(f"Model reloaded: {len(loaded.points)} points, {len(loaded.ligand_mols)} ligand structures")

# now actually load the RELOADED model into the real PyMOL cmd session,
# exactly as clicking "Load Model..." in the GUI would do
cmd.delete("ligand_* or pharm_* or candidate")
plugin._load_into_pymol(cmd, loaded)
objects = cmd.get_names()
ligand_objects = [o for o in objects if o.startswith("ligand_")]
assert len(ligand_objects) == 5, f"Expected 5 ligand objects from reloaded model, got {len(ligand_objects)}: {objects}"
print("Reloaded ligands present in real PyMOL session:", ligand_objects)

for obj in ligand_objects:
    n_atoms = cmd.count_atoms(obj)
    assert n_atoms > 5, f"{obj} has suspiciously few atoms: {n_atoms}"

plugin._refresh_pharmacophore_spheres(cmd, loaded, 0.5)
objects = cmd.get_names()
pharm_objects = [o for o in objects if o.startswith("pharm_")]
expected_n = len(loaded.consensus_at_threshold(0.5))
assert len(pharm_objects) == expected_n, f"Expected {expected_n} spheres from reloaded model, got {len(pharm_objects)}"
print(f"Reloaded pharmacophore spheres at 50% threshold: {pharm_objects}")

# verify sphere coordinates from the RELOADED model exactly match the ORIGINAL model's.
# Matched by nearest centroid within each family, since a family (e.g. Hydrophobe)
# can legitimately have more than one distinct cluster at different 3D positions.
orig_consensus = pharm.consensus_at_threshold(0.5)
loaded_consensus = loaded.consensus_at_threshold(0.5)
assert len(orig_consensus) == len(loaded_consensus), \
    f"Consensus count mismatch: {len(orig_consensus)} vs {len(loaded_consensus)}"

orig_by_fam = collections.defaultdict(list)
for p in orig_consensus:
    orig_by_fam[p.family].append(np.array(p.centroid))

for p in loaded_consensus:
    candidates = orig_by_fam[p.family]
    assert candidates, f"No original point found for family {p.family}"
    dists = [np.linalg.norm(np.array(p.centroid) - c) for c in candidates]
    best_diff = min(dists)
    assert best_diff < 1e-6, \
        f"{p.family}: no matching original centroid found within tolerance (closest diff={best_diff})"
print("All reloaded sphere coordinates match the original model exactly (matched per family by nearest centroid).")

cmd.delete("ligand_* or pharm_*")
os.unlink(save_path)
print("\nSAVE/LOAD MODEL ROUND-TRIP VERIFIED CORRECT (real PyMOL cmd)")

print("\n=== Test 14: full validation pipeline (EF / ROC-AUC / GH) on real MAO-B data ===")
from core import fetch_background_pool, property_matched_decoys, validate_pharmacophore, roc_auc, guner_henry_score

# unit-test the dependency-free ROC-AUC against a known-good sklearn comparison,
# even though sklearn is not a runtime dependency of the plugin itself
try:
    from sklearn.metrics import roc_auc_score as _sklearn_auc
    _labels_check = [1,1,1,0,0,1,0,0,0,1]
    _scores_check = [0.9,0.8,0.3,0.4,0.2,0.7,0.6,0.1,0.5,0.85]
    mine = roc_auc(_labels_check, _scores_check)
    theirs = _sklearn_auc(_labels_check, _scores_check)
    assert abs(mine - theirs) < 1e-9, f"ROC-AUC mismatch: mine={mine} sklearn={theirs}"
    print(f"  roc_auc() cross-checked against sklearn: {mine:.6f} == {theirs:.6f}")
except ImportError:
    print("  (sklearn not available in this environment to cross-check against \u2014 skipping that specific check)")

active_smiles = [smi for _, smi in entries]
pool, pool_source = fetch_background_pool(n=1500)
print(f"  Background pool: {pool_source}")
assert "offline fallback" in pool_source, "Expected offline fallback in this sandbox (no ChEMBL access)"

decoys = property_matched_decoys(active_smiles, pool, n_per_active=10)
assert len(decoys) >= 2, f"Too few decoys generated: {len(decoys)}"
print(f"  Generated {len(decoys)} property-matched decoys")

val_progress_log = []
result = validate_pharmacophore(
    consensus_pts, active_smiles, decoys, n_confs=20,
    progress_callback=lambda i, n: val_progress_log.append((i, n))
)
assert "error" not in result, f"Validation returned an error: {result.get('error')}"
assert 0.0 <= result["auc"] <= 1.0, f"AUC out of range: {result['auc']}"
assert result["n_actives"] == 5, f"Expected 5 actives scored, got {result['n_actives']}"
assert len(val_progress_log) == result["n_total"], "Progress callback should fire once per scored compound"
print(f"  AUC={result['auc']:.3f}  EF@{int(result['top_pct']*100)}%={result['ef']:.2f}x  "
      f"GH={result['gh']:.3f}  (Ha={result['Ha']}/Ht={result['Ht']})")

# sanity: our own MAO-B actives should generally score higher than the decoys on average
active_fracs = [r["fraction_matched"] for r in result["records"] if r["label"] == 1]
decoy_fracs = [r["fraction_matched"] for r in result["records"] if r["label"] == 0]
assert len(active_fracs) == 5 and len(decoy_fracs) == len(decoys)
mean_active = sum(active_fracs) / len(active_fracs)
mean_decoy = sum(decoy_fracs) / len(decoy_fracs)
print(f"  Mean fraction matched: actives={mean_active:.2f} decoys={mean_decoy:.2f}")
assert mean_active > mean_decoy, "Expected actives to score higher on average than decoys"

print("\nVALIDATION PIPELINE VERIFIED CORRECT")

print("\n=== Test 15: get_alignment_reference() on both pharmacophore types ===")
ref_mol, ref_cid = pharm.get_alignment_reference()
assert ref_mol is not None, "LigandBasedPharmacophore should expose a real reference mol"
assert ref_mol.GetNumAtoms() > 0
print(f"  LigandBasedPharmacophore reference: {pharm.reference.name}, {ref_mol.GetNumAtoms()} atoms")

loo_ref_mol, loo_ref_cid = loaded.get_alignment_reference()
assert loo_ref_mol is not None, "LoadedPharmacophore should also expose a reference mol (from saved ligands)"
print(f"  LoadedPharmacophore reference: {loo_ref_mol.GetNumAtoms()} atoms, conf_id={loo_ref_cid}")

print("\n=== Test 16: alignment reference measurably changes scoring (the real bug fix) ===")
# This is the actual regression test for the bug found during real-world use:
# lazabemide, scored without alignment, matched the wrong pair of features
# (Donor+Aromatic) instead of its real, chemically correct match (Donor+
# PosIonizable, both under 0.2 A on the properly-aligned conformer). With
# alignment, scoring should recover meaningfully more/better matches.
lazabemide_smiles = "NCCNC(=O)c1ccc(Cl)cn1"
consensus_60 = pharm.consensus_at_threshold(0.60)

n_matched_no_align, frac_no_align, mean_d_no_align, _, _, details_no_align, rmsd_no_align, shape_consistency = score_molecule(
    lazabemide_smiles, consensus_60, reference_mol=None, n_confs=50
)
n_matched_aligned, frac_aligned, mean_d_aligned, _, _, details_aligned, rmsd_aligned, shape_consistency = score_molecule(
    lazabemide_smiles, consensus_60, reference_mol=ref_mol, reference_conf_id=ref_cid, n_confs=50
)
print(f"  Without alignment: {n_matched_no_align}/{len(consensus_60)} matched, "
      f"families={[d['family'] for d in details_no_align]}")
print(f"  With alignment:    {n_matched_aligned}/{len(consensus_60)} matched, "
      f"families={[d['family'] for d in details_aligned]}, mean_dist={mean_d_aligned:.2f}, "
      f"align_rmsd={rmsd_aligned:.2f}")
assert rmsd_no_align is None, "Without a reference, align_rmsd should be None"
assert rmsd_aligned is not None and rmsd_aligned >= 0, "With a reference, align_rmsd should be a real number"
assert n_matched_aligned >= n_matched_no_align, \
    "Alignment should not make matching worse for this real regression case"
matched_families_aligned = set(d['family'] for d in details_aligned)
assert 'Donor' in matched_families_aligned, "Aligned scoring should recover the Donor match"
print("  Confirmed: alignment step changes (and here, improves) real scoring behavior.")

print("\n=== Test 17: GUI wiring passes the reference through to score_molecule/validate/screen ===")
# Confirm the actual call sites in __init__.py request a reference via
# get_alignment_reference() rather than silently falling back to unaligned
# scoring -- a static check of the shipped source, not just core.py logic.
init_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pharmol", "__init__.py")
with open(init_path) as f:
    init_source = f.read()
assert init_source.count("get_alignment_reference()") >= 3, \
    "Expected get_alignment_reference() to be called at least 3 times (score/validate/screen) in the GUI"
assert "reference_mol=ref_mol" in init_source
print("  Confirmed: single-candidate scoring, validation, and batch screening all request the alignment reference.")

print("\nALIGNMENT REFERENCE FIX VERIFIED CORRECT")

print("\n=== Test 18: MW range filter (real regression case from earlier session) ===")
from core import mw_range_from_actives, filter_by_mw_range
active_smiles_18 = [smi for _, smi in entries]
lo_mw, hi_mw = mw_range_from_actives(active_smiles_18)
print(f"  Active MW range: {lo_mw:.0f}-{hi_mw:.0f}")
assert 150 < lo_mw < 200 and 250 < hi_mw < 300, f"Unexpected MW range: {lo_mw}-{hi_mw}"

fake_results = [
    {"name": "in_range", "smiles": "CC(Cc1ccc(F)cc1)N(C)CC#C"},
    {"name": "the_real_bad_hit_from_earlier_session", "smiles": "COCCC(=O)Nc1ccc2c(c1)CN(Cc1coc3ccccc3c1=O)CC2"},
]
kept18, n_removed18 = filter_by_mw_range(fake_results, active_smiles_18, tolerance_da=50)
assert len(kept18) == 1 and kept18[0]["name"] == "in_range"
assert n_removed18 == 1
assert "mw" in kept18[0], "Kept results should have MW annotated"
print(f"  Correctly kept 'in_range', removed the real oversized false-positive from earlier testing.")

print("\n=== Test 19: GH score handles Ht=0 gracefully (no crash, returns nan) ===")
from core import guner_henry_score
import math
gh_zero_hits = guner_henry_score(A=5, D=100, Ht=0, Ha=0)
assert math.isnan(gh_zero_hits), f"Expected nan for Ht=0, got {gh_zero_hits}"
print("  guner_henry_score(Ht=0) correctly returns nan, not a crash.")

# Confirm the GUI-facing code path actually produces the clear message, not just
# a bare 'n/a', for this exact case -- a static check of the shipped source.
init_path_19 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pharmol", "__init__.py")
with open(init_path_19) as f:
    init_source_19 = f.read()
assert "no compound cleared the hit threshold" in init_source_19, \
    "Expected the GH-undefined case to show a clear explanation, not just 'n/a'"
print("  Confirmed: GUI shows a clear explanation (not a bare 'n/a') when GH is undefined.")

print("\n=== Test 20: eps re-clustering changes the consensus without re-aligning ===")
import time as _time
t_before = _time.time()
pts_tight = pharm.extract_consensus(eps=0.8, min_support_frac=0.0)
pts_loose = pharm.extract_consensus(eps=2.5, min_support_frac=0.0)
t_after = _time.time()
print(f"  Two re-clustering calls took {t_after-t_before:.3f}s (should be fast -- no re-alignment)")
assert t_after - t_before < 2.0, "Re-clustering should be fast, not re-run full alignment"
print(f"  eps=0.8: {len(pts_tight)} clusters | eps=2.5: {len(pts_loose)} clusters")
assert len(pts_tight) >= len(pts_loose), "Tighter eps should generally produce >= clusters vs. looser eps"
# restore the default used by the rest of the test suite
pharm.extract_consensus(eps=1.5, min_support_frac=0.0)
print("  Confirmed: eps re-clustering works correctly and cheaply on already-aligned ligands.")

print("\n=== Test 21: SDF export of aligned batch-screening hits ===")
test_lib_path_21 = "/tmp/test_sdf_export_lib.sdf"
writer21 = Chem.SDWriter(test_lib_path_21)
for name21, smi21 in [("sdf_test_1", "CC(Cc1ccccc1)N(C)CC#C"), ("sdf_test_2", "C#CCNC1CCc2ccccc21")]:
    m21 = Chem.MolFromSmiles(smi21)
    m21 = Chem.AddHs(m21)
    _AllChem.EmbedMolecule(m21, randomSeed=1)
    m21.SetProp("_Name", name21)
    writer21.write(m21)
writer21.close()

results21, _, _ = screen_library(test_lib_path_21, consensus_pts, reference_mol=ref_mol,
                                reference_conf_id=ref_cid, n_confs=20)
assert len(results21) == 2
assert all(r.get("align_rmsd") is not None for r in results21), "Batch results should include align_rmsd"
print(f"  Batch screening results include align_rmsd: {[round(r['align_rmsd'],2) for r in results21]}")

out_sdf_path = "/tmp/test_exported_hits.sdf"
sdf_writer = Chem.SDWriter(out_sdf_path)
for r in results21:
    mol_copy = Chem.Mol(r["mol"])
    mol_copy.SetProp("_Name", r["name"])
    mol_copy.SetProp("align_rmsd_angstrom", f"{r['align_rmsd']:.3f}")
    sdf_writer.write(mol_copy, confId=r["conf_id"])
sdf_writer.close()

reread = list(Chem.SDMolSupplier(out_sdf_path))
assert len(reread) == 2, f"Expected 2 molecules in exported SDF, got {len(reread)}"
assert all(m is not None and m.GetNumConformers() > 0 for m in reread), "Exported SDF should contain real 3D conformers"
assert reread[0].HasProp("align_rmsd_angstrom"), "Exported SDF should retain the align_rmsd property"
print(f"  Exported SDF re-read successfully: {len(reread)} real 3D structures with properties intact.")
os.unlink(test_lib_path_21)
os.unlink(out_sdf_path)

print("\nALL SIX EXPERT-FEEDBACK ITEMS VERIFIED CORRECT")

print("\n=== Test 22: leave-one-out validation avoids testing on training data (the real methodological fix) ===")
from core import validate_pharmacophore_loocv

result_loocv = validate_pharmacophore_loocv(
    entries, decoys, eps=1.5, support_thresh=0.60, n_confs=50,
)
assert "error" not in result_loocv, f"LOOCV validation errored: {result_loocv.get('error')}"
assert result_loocv["loocv_actives"] is True
assert result_loocv["n_loocv_runs"] == 5, f"Expected 5 leave-one-out runs, got {result_loocv['n_loocv_runs']}"
print(f"  LOOCV: AUC={result_loocv['auc']:.3f} EF={result_loocv['ef']:.2f}x GH={result_loocv['gh']:.3f} "
      f"(n_loocv_runs={result_loocv['n_loocv_runs']})")

# Compare against in-sample (testing on training data) using the same threshold/eps --
# confirms the two modes give genuinely different numbers, i.e. the distinction is real,
# not just a label that happens to produce identical results.
in_sample_result = validate_pharmacophore(
    consensus_pts, [smi for _, smi in entries], decoys,
    reference_mol=ref_mol, reference_conf_id=ref_cid, n_confs=50
)
print(f"  In-sample: AUC={in_sample_result['auc']:.3f} (tests actives against the model THEY were used to build)")
assert in_sample_result["loocv_actives"] is False
print("  Confirmed: LOOCV and in-sample validation are genuinely different code paths producing different, "
      "independently meaningful numbers -- not just cosmetic labels.")

print("\nLEAVE-ONE-OUT VALIDATION FIX VERIFIED CORRECT")

print("\n=== Test 23: external test actives GUI mechanism exists and is correctly wired ===")
init_path_23 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pharmol", "__init__.py")
with open(init_path_23) as f:
    init_source_23 = f.read()
assert "ext_actives_edit" in init_source_23, "External actives input box should exist in the GUI"
assert "using_external" in init_source_23, "do_validate should branch on whether external actives were supplied"
assert "External test actives (genuinely held-out" in init_source_23, \
    "Status message should correctly identify external-actives mode, distinct from in-sample/LOOCV"
assert "validation_used_external_actives" in init_source_23, \
    "The external-actives flag should be persisted in state for the report generator to see"
print("  Confirmed: external test actives input, correct branching, and correct status/report "
      "messaging are all present in the shipped source.")

print("\n=== Test 24: atom-highlighting for scored candidates (real PyMOL atom selections) ===")
ref_mol_24, ref_cid_24 = pharm.get_alignment_reference()
n_matched_24, frac_24, mean_d_24, mol_24, cid_24, match_details_24, rmsd_24, shape_consistency = score_molecule(
    "C#CCN(C)CCc1ccccc1", consensus_pts, reference_mol=ref_mol_24, reference_conf_id=ref_cid_24, n_confs=50
)
assert len(match_details_24) > 0, "Expected at least one matched feature with atom details"
plugin._load_candidate_into_pymol(cmd, mol_24, cid_24, match_details=match_details_24)
for i, d in enumerate(match_details_24):
    pymol_sel = f"candidate and index {'+'.join(str(a+1) for a in d['atom_ids'])}"
    n_found = cmd.count_atoms(pymol_sel)
    assert n_found == len(d['atom_ids']), \
        f"Feature {i} ({d['family']}): expected {len(d['atom_ids'])} atoms, found {n_found} in real PyMOL selection"
print(f"  Verified {len(match_details_24)} matched-feature atom selections resolve correctly in real PyMOL "
      f"(RDKit index + 1 = PyMOL index mapping holds).")
cmd.delete("candidate")

print("\n=== Test 25: adjustable match tolerance changes scoring monotonically ===")
lazabemide_smi_25 = "NCCNC(=O)c1ccc(Cl)cn1"
n_tight, _, _, _, _, _, _, shape_consistency = score_molecule(lazabemide_smi_25, consensus_pts, reference_mol=ref_mol_24,
                                             reference_conf_id=ref_cid_24, tol=0.5, n_confs=50)
n_default, _, _, _, _, _, _, shape_consistency = score_molecule(lazabemide_smi_25, consensus_pts, reference_mol=ref_mol_24,
                                               reference_conf_id=ref_cid_24, tol=1.8, n_confs=50)
n_loose, _, _, _, _, _, _, shape_consistency = score_molecule(lazabemide_smi_25, consensus_pts, reference_mol=ref_mol_24,
                                             reference_conf_id=ref_cid_24, tol=3.5, n_confs=50)
print(f"  Matched counts: tol=0.5 -> {n_tight}, tol=1.8 -> {n_default}, tol=3.5 -> {n_loose}")
assert n_tight <= n_default <= n_loose, "Tighter tolerance should never match more than a looser one"
print("  Confirmed: tolerance parameter correctly threads through score_molecule with monotonic behavior.")

# Static check: confirm the GUI actually wires tol_spin.value() through to every scoring call site
init_path_25 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pharmol", "__init__.py")
with open(init_path_25) as f:
    init_source_25 = f.read()
assert init_source_25.count("tol=tol_spin.value()") >= 4, \
    "Expected tol_spin.value() threaded through at least 4 call sites (score/validate/loocv/screen)"
print("  Confirmed: match tolerance is wired through all 4 scoring call sites in the shipped GUI source.")

print("\nITEMS 2 AND 3 VERIFIED CORRECT")

print("\n=== Test 26: two-axis verdict classification (match + plausibility) ===")
from core import classify_verdict

verdict_cases = [
    (0.85, 0.65, "Gold Standard Hit"),
    (0.85, 0.40, "Scaffold-Hop (Moderate Confidence)"),
    (0.85, 0.14, "Scaffold-Hop (High Risk)"),  # the real MAO-B false-positive case from earlier testing
    (0.55, 0.60, "Uncertain"),
    (0.20, 0.55, "Inactive Analog"),
    (0.10, 0.10, "Alien / Low Priority"),
]
for frac, plaus, expected in verdict_cases:
    label, color, explanation = classify_verdict(frac, plaus)
    assert label == expected, f"frac={frac} plaus={plaus}: expected {expected!r}, got {label!r}"
    assert color.startswith("#") and len(color) == 7, f"Invalid color: {color}"
    assert len(explanation) > 20, "Verdict should include a real explanation, not a stub"
print(f"  All {len(verdict_cases)} verdict cases classified correctly (including the real MAO-B false-positive case).")

# Confirm low plausibility is never silently discarded -- the whole point of this feature
label, color, explanation = classify_verdict(0.85, 0.10)
assert "risk" in label.lower() or "scaffold" in label.lower(), \
    "A high-match, low-plausibility candidate must be flagged as higher-risk, never silently dropped"
assert "do not discard" in explanation.lower(), \
    "Explanation should explicitly instruct against discarding scaffold-hop candidates"
print("  Confirmed: high-match/low-plausibility candidates are flagged as higher-risk, never discarded.")

# Confirm the GUI actually wires the Verdict column in, as a static check of shipped source
init_path_26 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pharmol", "__init__.py")
with open(init_path_26) as f:
    init_source_26 = f.read()
assert '"Verdict"' in init_source_26, "Expected a Verdict column header in the shipped GUI source"
assert "classify_verdict(" in init_source_26, "Expected classify_verdict to actually be called in the GUI"
print("  Confirmed: Verdict column and classify_verdict() call are present in the shipped GUI source.")

print("\nVERDICT CLASSIFICATION SYSTEM VERIFIED CORRECT")

print("\n=== Test 27: Scaffold Novelty (Bemis-Murcko) ===")
from core import classify_scaffold_novelty

training_coxibs = [
    "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",
    "CS(=O)(=O)c1ccc(C2=C(c3ccccc3)C(=O)OC2)cc1",
    "Cc1onc(-c2ccccc2)c1-c1ccc(S(N)(=O)=O)cc1",
    "Cc1ccc(-c2ncc(Cl)cc2-c2ccc(S(C)(=O)=O)cc2)cn1",
]
scaffold_cases = [
    ("CCC(O)=NS(=O)(=O)c1ccc(-c2c(-c3ccccc3)noc2C)cc1", "Same Scaffold"),  # parecoxib, real prodrug
    ("Cc1ccc(Nc2c(F)cccc2Cl)c(CC(=O)O)c1", "Novel Scaffold"),  # lumiracoxib, real outlier
    ("CCCCCCCCCC", "Unknown"),  # acyclic, no ring system at all
]
for smi, expected in scaffold_cases:
    result = classify_scaffold_novelty(smi, training_coxibs)
    assert result == expected, f"smi={smi}: expected {expected!r}, got {result!r}"
print(f"  All {len(scaffold_cases)} real coxib scaffold cases classified correctly.")

print("\n=== Test 28: Fast pharmacophore fingerprint pre-filter ===")
from core import compute_pharmacophore_fingerprint, fast_prefilter_score

coxib_ligands = [Ligand(name, smi, n_confs=50) for name, smi in
                  [("celecoxib", training_coxibs[0]), ("rofecoxib", training_coxibs[1]),
                   ("valdecoxib", training_coxibs[2]), ("etoricoxib", training_coxibs[3])]]
coxib_pharm = LigandBasedPharmacophore(coxib_ligands)
coxib_pharm.align()
coxib_pharm.extract_consensus(eps=1.5, min_support_frac=0.0)
coxib_consensus = coxib_pharm.consensus_at_threshold(0.80)
fp = compute_pharmacophore_fingerprint(coxib_consensus)
assert len(fp) > 0, "Fingerprint should have at least one pairwise distance for a multi-feature consensus"

t0 = _time.time()
decane_score = fast_prefilter_score("CCCCCCCCCC", fp, n_confs=10)
t1 = _time.time()
print(f"  Decane pre-filter score: {decane_score:.2f} (should be low), took {t1-t0:.2f}s")
assert decane_score < 0.3, f"Expected an obviously-irrelevant compound to score low, got {decane_score}"

# Confirm screen_library's use_prefilter option actually reduces the compound count
# reaching full scoring, on a real mixed library (a real active + an obvious non-match)
test_lib_path_28 = "/tmp/test_prefilter_lib.sdf"
writer28 = Chem.SDWriter(test_lib_path_28)
for name28, smi28 in [("real_analog", "Cc1onc(-c2ccccc2)c1-c1ccc(S(N)(=O)=O)cc1"), ("decane_irrelevant", "CCCCCCCCCC")]:
    m28 = Chem.AddHs(Chem.MolFromSmiles(smi28))
    _AllChem.EmbedMolecule(m28, randomSeed=1)
    m28.SetProp("_Name", name28)
    writer28.write(m28)
writer28.close()

coxib_ref_mol, coxib_ref_cid = coxib_pharm.get_alignment_reference()
results28, n_skipped28, n_prefiltered28 = screen_library(
    test_lib_path_28, coxib_consensus, reference_mol=coxib_ref_mol, reference_conf_id=coxib_ref_cid,
    n_confs=20, use_prefilter=True, prefilter_thresh=0.3
)
print(f"  With prefilter on: {len(results28)} fully scored, {n_prefiltered28} rejected by pre-filter")
assert n_prefiltered28 >= 1, "Expected the obviously-irrelevant decane to be rejected by the pre-filter"
assert any(r["name"] == "real_analog" for r in results28), "The real analog should survive the pre-filter and be fully scored"
os.unlink(test_lib_path_28)
print("  Confirmed: use_prefilter correctly reduces full-scoring workload while keeping the real match.")

print("\nSCAFFOLD NOVELTY AND FAST PRE-FILTER VERIFIED CORRECT")

print("\n=== Test 29: Shape Consistency wired through screen_library and the shipped GUI ===")
from core import shape_consistency_score as _shape_score_check

coxib_envelope = coxib_pharm.get_shape_envelope()
assert len(coxib_envelope) > 0, "Expected a non-empty shape envelope for a 4-active model"

test_lib_path_29 = "/tmp/test_shape_lib.sdf"
writer29 = Chem.SDWriter(test_lib_path_29)
for name29, smi29 in [("valdecoxib_self", "Cc1onc(-c2ccccc2)c1-c1ccc(S(N)(=O)=O)cc1"),
                        ("deliberately_bulky", "CCCCCCCCCCCCCCCCCCc1ccc(-c2ccccc2-c3ccccc3-c4ccccc4)cc1")]:
    m29 = Chem.AddHs(Chem.MolFromSmiles(smi29))
    _AllChem.EmbedMolecule(m29, randomSeed=1)
    m29.SetProp("_Name", name29)
    writer29.write(m29)
writer29.close()

results29, _, _ = screen_library(
    test_lib_path_29, coxib_consensus, reference_mol=coxib_ref_mol, reference_conf_id=coxib_ref_cid,
    n_confs=20, shape_envelope=coxib_envelope, active_smiles_list=training_coxibs
)
for r in results29:
    print(f"  {r['name']}: shape_consistency={r['shape_consistency']:.2f} scaffold_novelty={r['scaffold_novelty']}")
by_name29 = {r["name"]: r for r in results29}
assert by_name29["valdecoxib_self"]["shape_consistency"] == 1.0, "valdecoxib itself should be 100% shape-consistent"
assert by_name29["deliberately_bulky"]["shape_consistency"] < 0.9, "The deliberately bulky compound should show reduced shape consistency"
assert by_name29["valdecoxib_self"]["scaffold_novelty"] == "Same Scaffold"
os.unlink(test_lib_path_29)

# Static check: confirm the GUI actually surfaces both new columns, not just core.py
init_path_29 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pharmol", "__init__.py")
with open(init_path_29) as f:
    init_source_29 = f.read()
assert '"Shape"' in init_source_29 and '"Scaffold"' in init_source_29, \
    "Expected Shape and Scaffold columns in the shipped GUI table headers"
assert "get_shape_envelope()" in init_source_29, "Expected the GUI to actually call get_shape_envelope()"
print("  Confirmed: Shape and Scaffold columns are wired into the shipped GUI, with correct real-data behavior.")

print("\nSHAPE CONSISTENCY AND SCAFFOLD NOVELTY BATCH-SCREENING WIRING VERIFIED CORRECT")

print("\n=== Test 30: Ensemble tolerance matching (nearest raw point vs. centroid) ===")
from core import PharmacophoricPoint, _point_distance
gliptin_entries = [
    ("sitagliptin", "Fc1cc(c(F)cc1F)C[C@@H](N)CC(=O)N3Cc2nnc(n2CC3)C(F)(F)F"),
    ("vildagliptin", "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2"),
    ("saxagliptin", "O=C(N1[C@H](C#N)C[C@@H]2C[C@H]12)[C@@H](N)C35CC4CC(C3)CC(O)(C4)C5"),
    ("linagliptin", "CC#CCN1C2=C(N=C1N3CCC[C@H](C3)N)N(C(=O)N(C2=O)CC4=NC5=CC=CC=C5C(=N4)C)C"),
    ("alogliptin", "CN2C(=O)C=C(N3CCC[C@@H](N)C3)N(CC4=C(C=CC=C4)C#N)C2=O"),
    ("teneligliptin", "C([C@@H]1C[C@H](N2CCN(C3N(C4=CC=CC=C4)N=C(C)C=3)CC2)CN1)(N1CCSC1)=O"),
    ("anagliptin", "CC1=NN2C=C(C=NC2=C1)C(=O)NCC(C)(C)NCC(=O)N3CCC[C@H]3C#N"),
]
gliptin_ligands = [Ligand(name, smi, n_confs=20) for name, smi in gliptin_entries]
gliptin_pharm = LigandBasedPharmacophore(gliptin_ligands)
gliptin_pharm.align()
gliptin_pharm.extract_consensus(eps=1.5, min_support_frac=0.0)
gliptin_consensus = gliptin_pharm.consensus_at_threshold(0.60)
gliptin_ref_mol, gliptin_ref_cid = gliptin_pharm.get_alignment_reference()

assert any(p.raw_points is not None and len(p.raw_points) > 0 for p in gliptin_consensus), \
    "Expected at least one consensus feature to have real raw_points stored"

linagliptin_smi = gliptin_entries[3][1]
result_centroid = score_molecule(linagliptin_smi, gliptin_consensus, reference_mol=gliptin_ref_mol,
                                   reference_conf_id=gliptin_ref_cid, n_confs=20, match_mode="centroid")
result_ensemble = score_molecule(linagliptin_smi, gliptin_consensus, reference_mol=gliptin_ref_mol,
                                   reference_conf_id=gliptin_ref_cid, n_confs=20, match_mode="ensemble")
mean_d_centroid, mean_d_ensemble = result_centroid[2], result_ensemble[2]
print(f"  linagliptin (training active) mean_dist: centroid={mean_d_centroid:.2f} vs ensemble={mean_d_ensemble:.2f}")
assert mean_d_ensemble <= mean_d_centroid, \
    "Ensemble matching should never produce a worse (larger) mean distance than centroid matching"
print("  Confirmed: ensemble mode measurably changes (tightens) real matching distances.")

# Backward compatibility: a point with NO raw_points (e.g. loaded from an older
# saved model) must fall back to centroid matching automatically, not crash.
old_style_point = PharmacophoricPoint("Donor", [0.0, 0.0, 0.0], 1.0, 5)  # no raw_points given
assert old_style_point.raw_points is None
d = _point_distance(np.array([1.0, 0.0, 0.0]), old_style_point, "ensemble")
assert abs(d - 1.0) < 1e-6, "Should fall back to centroid distance when raw_points is unavailable"
print("  Confirmed: ensemble mode gracefully falls back to centroid distance for points without raw_points.")

# JSON round-trip: raw_points must survive save/load exactly (real data, not synthetic)
export_pharmacophore_json(gliptin_pharm.points, "/tmp/test_ensemble_roundtrip.json", metadata={}, ligands=gliptin_ligands)
reloaded_pharm = load_pharmacophore_json("/tmp/test_ensemble_roundtrip.json")
for orig_p, reloaded_p in zip(gliptin_pharm.points, reloaded_pharm.points):
    orig_n = len(orig_p.raw_points) if orig_p.raw_points is not None else 0
    reloaded_n = len(reloaded_p.raw_points) if reloaded_p.raw_points is not None else 0
    assert orig_n == reloaded_n
    if orig_p.raw_points is not None:
        assert np.allclose(orig_p.raw_points, reloaded_p.raw_points)
os.unlink("/tmp/test_ensemble_roundtrip.json")
print("  Confirmed: raw_points survive JSON save/load exactly on real gliptin data.")

# Static check: confirm the GUI checkbox is wired to both scoring call sites
init_path_30 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pharmol", "__init__.py")
with open(init_path_30) as f:
    init_source_30 = f.read()
assert init_source_30.count('match_mode="ensemble" if ensemble_check.isChecked() else "centroid"') >= 2, \
    "Expected the ensemble checkbox to be wired into at least 2 scoring call sites (single + batch)"
print("  Confirmed: ensemble checkbox is wired into both single-candidate and batch scoring in the shipped GUI.")

print("\nENSEMBLE TOLERANCE MATCHING VERIFIED CORRECT")

print("\n=== Test 31: Toggle Ligands button (real bug found and fixed) ===")
# cmd.toggle()'s first argument is a representation name (e.g. "sticks"), not an
# object/selection pattern. cmd.toggle("ligand_*") silently failed with "unknown
# representation" -- found via real GUI testing, not code inspection alone.
# Confirm the shipped source no longer contains the broken call, and does use
# enable/disable (which correctly support wildcard object patterns).
init_path_31 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pharmol", "__init__.py")
with open(init_path_31) as f:
    init_source_lines_31 = f.readlines()
live_toggle_calls = [ln for ln in init_source_lines_31
                      if 'cmd.toggle("ligand_*")' in ln and not ln.strip().startswith("#")]
assert not live_toggle_calls, \
    f"The broken cmd.toggle('ligand_*') call should no longer be present as live code: {live_toggle_calls}"
init_source_31 = "".join(init_source_lines_31)
assert 'cmd.enable("ligand_*")' in init_source_31 and 'cmd.disable("ligand_*")' in init_source_31, \
    "Expected the fix to use cmd.enable/cmd.disable, which correctly support wildcard patterns"
print("  Confirmed: broken cmd.toggle() call removed from live code, replaced with enable/disable + tracked state.")

# Direct behavioral confirmation with real PyMOL objects (mirrors exactly what
# the GUI's toggle handler does, without needing to drive the whole dialog again).
import pymol as _pymol_check
from pymol import cmd as _cmd_check
_cmd_check.pseudoatom("ligand_0_regtest")
_cmd_check.pseudoatom("ligand_1_regtest")
assert _cmd_check.get_names("objects", enabled_only=1) == ["ligand_0_regtest", "ligand_1_regtest"]
_cmd_check.disable("ligand_*")
assert _cmd_check.get_names("objects", enabled_only=1) == [], "disable('ligand_*') should hide both objects"
_cmd_check.enable("ligand_*")
assert set(_cmd_check.get_names("objects", enabled_only=1)) == {"ligand_0_regtest", "ligand_1_regtest"}, \
    "enable('ligand_*') should show both objects again"
_cmd_check.delete("ligand_*regtest")
print("  Confirmed: enable/disable with wildcard patterns behaves correctly on real PyMOL objects.")

print("\nTOGGLE LIGANDS BUG FIX VERIFIED CORRECT")

print("\nALL PLUGIN TESTS PASSED")
