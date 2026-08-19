# Contributing to PhArMol

Thanks for considering a contribution. This is a young project and contributions of
any size are welcome — bug reports, documentation fixes, and new features alike.

## Getting set up

Follow the [Quick Start](README.md#quick-start) to get a working environment, then:

```bash
git clone https://github.com/OWNER/pharmol.git
cd pharmol
pip install -e .
python test_plugin.py
```

The test suite runs against a real PyMOL session (`-cq` mode) and real RDKit
computation — no mocked chemistry. It should complete cleanly before and after any
change you make.

## Reporting a bug

Please include:
- The exact steps that reproduce it
- What you expected vs. what happened
- Your PyMOL/RDKit versions and how they were installed (`conda`, `pip`, Homebrew, etc.)
- The full traceback if there is one

If you've found something that produces a real, incorrect scientific result (not just
a crash), please say so explicitly and include the compounds/SMILES involved — several
real bugs in this project were found exactly this way, not through code review.

## Making a change

1. Open an issue first for anything beyond a small fix, so the approach can be
   discussed before you invest time in it.
2. Add or update a test in `test_plugin.py` that would have caught the bug, or that
   verifies the new feature against real data — not synthetic/toy data where possible.
3. Run the full test suite locally before opening a pull request.
4. Keep the PR focused on one change; separate unrelated fixes into separate PRs.

## Code style

- Match the existing style in `core.py` and `__init__.py` — plain, explicit code with
  docstrings that explain *why*, not just *what*, especially for anything scientifically
  non-obvious.
- New scientific claims or heuristics (thresholds, defaults, approximations) should be
  documented with the reasoning behind the chosen value, and verified against real data
  before being adopted, the same way existing features in this project were.

## Questions

Open an issue — there's no separate mailing list or chat for this project yet.
