# PhArMol

**Ligand-based pharmacophore modeling, natively inside PyMOL.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

PhArMol takes a handful of known active compounds against a drug target, builds a
consensus 3D pharmacophore from them, and lets you score new candidates — one at a
time or by the thousands — against that pattern, with built-in statistical validation
at every step. It runs **inside your already-open PyMOL session** — there is no
separate application and no embedded viewer. When you click Analyze, the aligned
ligands and pharmacophore spheres load directly into PyMOL's own live 3D viewport,
rotatable and zoomable with your mouse exactly like any other PyMOL object.

## Why PhArMol

- **Runs natively in PyMOL** — not a standalone app wrapping a separate viewer
- **Statistically honest by default** — leave-one-out cross-validation, permutation
  significance testing, and a dedicated external-test-actives workflow, so validation
  results reflect genuine generalization rather than testing a model against its own
  training data
- **A multi-signal Verdict, not just a geometric score** — combines pharmacophore
  match, chemical plausibility, physical shape consistency, and scaffold novelty into
  one plain-language judgment, deliberately never using any single signal as a hard
  filter
- **RDKit-only** — no scikit-learn, no scipy, no proprietary dependencies

## Quick start

```bash
conda create -n pharm -c conda-forge python=3.10 pymol-open-source rdkit
conda activate pharm
pymol
```

Then, in PyMOL: **Plugin → Plugin Manager → Install New Plugin**, and select this
repository's `pharmol/` folder. Open it from **Plugin → PhArMol**.

Click **Load MAO-B Example** on the Build Model tab for a working demonstration, or
follow the [full tutorial](docs/TUTORIAL.md) for a complete, realistic worked example.

Full installation options (including a fallback for an existing PyMOL install) are in
the [User Manual](docs/USER_MANUAL.md).

## Documentation

| Document | For |
|---|---|
| [User Manual](docs/USER_MANUAL.md) | Every input, parameter, and output, explained |
| [Tutorial](docs/TUTORIAL.md) | A complete worked example, start to finish, with real data |
| [Complete Documentation](docs/COMPLETE_DOCUMENTATION.md) | Architecture, methodology, and validation case studies |
| [CHANGELOG](CHANGELOG.md) | Version history |

## Repository layout

| Path | Purpose |
|---|---|
| `pharmol/__init__.py` | Plugin entry point: dialog UI and all PyMOL `cmd` integration |
| `pharmol/core.py` | The pharmacophore science: alignment, clustering, scoring, validation |
| `pharmol/launcher.py` | Powers the optional `pharmol-gui` one-command launcher |
| `test_plugin.py` | End-to-end test suite, run against PyMOL's real `cmd` API |
| `examples/` | Real, verified example datasets (training actives, test actives, decoys) |

## Testing

31 automated tests, run against real PyMOL (`-cq` mode) and real RDKit computation —
no mocked chemistry. Every scoring function, GUI wiring path, and file-format
round-trip is checked directly. Run locally with:

```bash
python test_plugin.py
```

CI runs the same suite on every push and pull request (see badge above).

## Citing PhArMol

See [CITATION.cff](CITATION.cff), or cite directly:

> *(citation details to be finalized upon first tagged release)*

## Contributing

Contributions are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT — see [LICENSE](LICENSE).
