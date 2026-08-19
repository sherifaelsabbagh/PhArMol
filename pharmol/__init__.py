"""
PhArMol — Pharmacophore modeling in PyMOL. A real PyMOL plugin.

Unlike a standalone app, this plugin runs *inside* your already-open PyMOL session.
When you click Analyze, the aligned ligands and consensus pharmacophore spheres are
loaded directly into PyMOL's own live 3D viewport — the one you already have open —
so you rotate/zoom it exactly like any other PyMOL object, with your mouse, in real time.

Install: PyMOL > Plugin > Plugin Manager > Install New Plugin > choose this folder
(zipped) or this file's parent directory.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

_IMPORT_ERROR = None
try:
    from core import (Ligand, LigandBasedPharmacophore, score_molecule, screen_library,
                       diverse_top_hits, filter_by_mw_range, mw_range_from_actives,
                       export_pharmacophore_json, load_pharmacophore_json,
                       fetch_background_pool, property_matched_decoys,
                       fetch_chembl_inactives, load_decoys_from_file,
                       validate_pharmacophore, validate_pharmacophore_loocv,
                       applicability_domain_similarity, classify_verdict,
                       leave_one_out_stability, compare_pharmacophore_models,
                       FAMILY_COLORS, DEFAULT_N_CONFS)
except ImportError as e:
    _IMPORT_ERROR = e
    Ligand = LigandBasedPharmacophore = score_molecule = screen_library = FAMILY_COLORS = None
    diverse_top_hits = filter_by_mw_range = mw_range_from_actives = None
    export_pharmacophore_json = load_pharmacophore_json = None
    fetch_background_pool = property_matched_decoys = validate_pharmacophore = None
    validate_pharmacophore_loocv = None
    fetch_chembl_inactives = load_decoys_from_file = None
    applicability_domain_similarity = classify_verdict = None
    leave_one_out_stability = compare_pharmacophore_models = None
    DEFAULT_N_CONFS = 50


def _compat_enum(qtmods, *paths):
    """
    Return the first resolvable enum value from a list of dotted attribute
    paths, tried in order. Needed because PyMOL bundles either PyQt5 (flat
    enum style, e.g. QtCore.Qt.Horizontal) or PyQt6 (nested scoped-enum
    style, e.g. QtCore.Qt.Orientation.Horizontal) depending on platform and
    build. Verified against both real PyQt5 and real PyQt6 before shipping.
    """
    last_err = None
    for path in paths:
        try:
            parts = path.split(".")
            cur = qtmods[parts[0]]
            for p in parts[1:]:
                cur = getattr(cur, p)
            return cur
        except AttributeError as e:
            last_err = e
            continue
    raise last_err


def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('PhArMol', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None

EXAMPLE_TEXT = """selegiline, CC(Cc1ccccc1)N(C)CC#C
rasagiline, C#CCNC1CCc2ccccc21
pargyline, C#CCN(C)Cc1ccccc1
clorgyline, C#CCN(C)CCCOc1ccc(Cl)cc1Cl
lazabemide, NCCNC(=O)c1ccc(Cl)cn1"""


def run_plugin_gui():
    global dialog

    if _IMPORT_ERROR is not None:
        _show_missing_dependency_dialog(_IMPORT_ERROR)
        return

    if dialog is None:
        dialog = make_dialog()
    dialog.show()
    dialog.raise_()


def _show_missing_dependency_dialog(error):
    """
    Shown instead of a raw traceback when RDKit (or another required package)
    isn't installed in PyMOL's own Python environment. This is a common,
    expected failure mode — different PyMOL distributions (Homebrew, the
    official installer, conda-forge, pip) each bundle a separate, isolated
    Python, and dependencies have to be installed into that *specific* one.
    """
    from pymol.Qt import QtWidgets
    exe = sys.executable
    msg = (
        "This plugin needs RDKit installed in PyMOL's own Python environment, "
        "and it wasn't found.\n\n"
        f"Missing dependency error:\n{error}\n\n"
        f"PyMOL's Python interpreter is:\n{exe}\n\n"
        "Fix — run this in a terminal (not inside PyMOL):\n\n"
        f'"{exe}" -m pip install rdkit\n\n'
        "(If that path looks like a Homebrew/system PyMOL install rather than a "
        "conda or venv environment, this is a known rough edge — see the README's "
        "Installation section for the more portable conda-forge/pip route, "
        "recommended for anyone installing this plugin fresh.)"
    )
    box = QtWidgets.QMessageBox()
    box.setWindowTitle("PhArMol — Missing Dependency")
    box.setText(msg)
    box.exec()


def _make_roc_widget(QtWidgets, QtGui, QtCore):
    """
    A publication-quality ROC curve plot using Qt's own QPainter directly —
    no matplotlib, consistent with this plugin's other dependency choices.
    White background regardless of the surrounding dark app theme (matching
    standard figure conventions for papers/posters), with gridlines, full
    tick labels, proper axis titles, a legend, and a title. Exportable to a
    real PNG file via export_png(), rendered at higher resolution than the
    on-screen widget for print quality.

    Built as a factory function (not a module-level class) because it needs
    pymol.Qt's QtWidgets/QtGui, which are only importable once PyMOL/Qt has
    actually been initialised.
    """
    class RocCurveWidget(QtWidgets.QWidget):
        def __init__(self):
            super().__init__()
            self.fpr = None
            self.tpr = None
            self.auc = None
            self.n_actives = None
            self.n_decoys = None
            self.setMinimumHeight(320)
            self.setStyleSheet("background-color: white;")

        def set_curve(self, fpr, tpr, auc, n_actives=None, n_decoys=None):
            self.fpr, self.tpr, self.auc = fpr, tpr, auc
            self.n_actives, self.n_decoys = n_actives, n_decoys
            self.update()

        def clear_curve(self):
            self.fpr = self.tpr = self.auc = None
            self.n_actives = self.n_decoys = None
            self.update()

        def _render(self, painter, w, h, scale=1.0):
            def px(v):
                return v * scale

            try:
                hint = _compat_enum(
                    {"QtGui": QtGui}, "QtGui.QPainter.RenderHint.Antialiasing", "QtGui.QPainter.Antialiasing"
                )
                painter.setRenderHint(hint)
            except Exception:
                pass

            painter.fillRect(0, 0, w, h, QtGui.QColor("white"))

            margin_l, margin_b, margin_t, margin_r = px(58), px(52), px(38), px(24)
            plot_w = w - margin_l - margin_r
            plot_h = h - margin_t - margin_b
            if plot_w <= 10 or plot_h <= 10:
                return

            def to_xy(x, y):
                return (margin_l + x * plot_w, margin_t + (1 - y) * plot_h)

            title_font = QtGui.QFont("Arial", int(11 * scale))
            title_font.setBold(True)
            axis_font = QtGui.QFont("Arial", int(9 * scale))
            tick_font = QtGui.QFont("Arial", int(8 * scale))

            # Title
            painter.setFont(title_font)
            painter.setPen(QtGui.QColor("#0F172A"))
            painter.drawText(int(margin_l), int(px(20)), "ROC Curve — Pharmacophore Validation")

            # Gridlines + tick labels at 0.2 intervals
            painter.setFont(tick_font)
            grid_pen = QtGui.QPen(QtGui.QColor("#E2E8F0"))
            grid_pen.setWidthF(1.0)
            for i in range(6):
                v = i / 5.0
                gx, _ = to_xy(v, 0)
                _, gy = to_xy(0, v)
                painter.setPen(grid_pen)
                painter.drawLine(int(gx), int(margin_t), int(gx), int(margin_t + plot_h))
                painter.drawLine(int(margin_l), int(gy), int(margin_l + plot_w), int(gy))
                painter.setPen(QtGui.QColor("#334155"))
                label = f"{v:.1f}"
                painter.drawText(int(gx - px(8)), int(margin_t + plot_h + px(16)), label)
                painter.drawText(int(margin_l - px(28)), int(gy + px(3)), label)

            # Axis box
            axis_pen = QtGui.QPen(QtGui.QColor("#0F172A"))
            axis_pen.setWidthF(1.2)
            painter.setPen(axis_pen)
            painter.drawRect(QtCore.QRectF(margin_l, margin_t, plot_w, plot_h))

            # Axis titles
            painter.setFont(axis_font)
            painter.setPen(QtGui.QColor("#0F172A"))
            painter.drawText(
                int(margin_l + plot_w / 2 - px(45)), int(h - px(10)), "False Positive Rate"
            )
            painter.save()
            painter.translate(px(14), margin_t + plot_h / 2 + px(45))
            painter.rotate(-90)
            painter.drawText(0, 0, "True Positive Rate")
            painter.restore()

            # Diagonal reference line (random classifier)
            diag_pen = QtGui.QPen(QtGui.QColor("#94A3B8"))
            diag_pen.setWidthF(1.3)
            diag_pen.setStyle(_compat_enum(
                {"QtCore": QtCore}, "QtCore.Qt.PenStyle.DashLine", "QtCore.Qt.DashLine"
            ))
            painter.setPen(diag_pen)
            x0, y0 = to_xy(0, 0)
            x1, y1 = to_xy(1, 1)
            painter.drawLine(int(x0), int(y0), int(x1), int(y1))

            legend_lines = [("Random (AUC = 0.50)", "#94A3B8", True)]

            if self.fpr and self.tpr and len(self.fpr) > 1:
                curve_pen = QtGui.QPen(QtGui.QColor("#0891B2"))
                curve_pen.setWidthF(2.4)
                painter.setPen(curve_pen)
                points = [to_xy(x, y) for x, y in zip(self.fpr, self.tpr)]
                path = QtGui.QPainterPath()
                path.moveTo(*points[0])
                for pt in points[1:]:
                    path.lineTo(*pt)
                painter.drawPath(path)

                auc_txt = f"{self.auc:.2f}" if self.auc == self.auc else "n/a"
                legend_lines.insert(0, (f"Model (AUC = {auc_txt})", "#0891B2", False))

                # Legend box, upper-left inside the plot
                painter.setFont(tick_font)
                lx, ly = margin_l + px(10), margin_t + px(14)
                box_w, box_h = px(150), px(16 * len(legend_lines) + 8)
                painter.setBrush(QtGui.QColor(255, 255, 255, 230))
                painter.setPen(QtGui.QPen(QtGui.QColor("#CBD5E1")))
                painter.drawRect(QtCore.QRectF(lx, ly, box_w, box_h))
                for i, (label, color, dashed) in enumerate(legend_lines):
                    ly_i = ly + px(8) + i * px(16)
                    line_pen = QtGui.QPen(QtGui.QColor(color))
                    line_pen.setWidthF(2.0)
                    if dashed:
                        line_pen.setStyle(_compat_enum(
                            {"QtCore": QtCore}, "QtCore.Qt.PenStyle.DashLine", "QtCore.Qt.DashLine"
                        ))
                    painter.setPen(line_pen)
                    painter.drawLine(int(lx + px(6)), int(ly_i), int(lx + px(24)), int(ly_i))
                    painter.setPen(QtGui.QColor("#0F172A"))
                    painter.drawText(int(lx + px(30)), int(ly_i + px(3)), label)

                if self.n_actives is not None and self.n_decoys is not None:
                    painter.setPen(QtGui.QColor("#64748B"))
                    painter.drawText(
                        int(margin_l + plot_w - px(150)), int(margin_t + plot_h - px(8)),
                        f"n = {self.n_actives} actives, {self.n_decoys} decoys"
                    )
            else:
                painter.setFont(axis_font)
                painter.setPen(QtGui.QColor("#94A3B8"))
                painter.drawText(
                    int(margin_l + px(20)), int(margin_t + plot_h / 2),
                    "Run validation to see the ROC curve"
                )

        def paintEvent(self, event):
            painter = QtGui.QPainter(self)
            self._render(painter, self.width(), self.height(), scale=1.0)

        def export_png(self, path, width=1200, height=800):
            """Render at higher resolution than the on-screen widget, for a
            print-quality figure suitable for a paper, poster, or slide."""
            fmt = _compat_enum(
                {"QtGui": QtGui}, "QtGui.QImage.Format.Format_ARGB32", "QtGui.QImage.Format_ARGB32"
            )
            image = QtGui.QImage(width, height, fmt)
            image.fill(QtGui.QColor("white"))
            painter = QtGui.QPainter(image)
            scale = width / 420.0  # scale factor relative to a typical on-screen widget width
            self._render(painter, width, height, scale=scale)
            painter.end()
            return image.save(path)

    return RocCurveWidget()


def make_dialog():
    from pymol import cmd
    from pymol.Qt import QtWidgets, QtGui, QtCore

    dlg = QtWidgets.QDialog()
    dlg.setWindowTitle("PhArMol")
    dlg.resize(600, 720)

    qtmods = {"QtCore": QtCore, "QtWidgets": QtWidgets, "QtGui": QtGui}
    horizontal = _compat_enum(qtmods, "QtCore.Qt.Orientation.Horizontal", "QtCore.Qt.Horizontal")
    no_edit = _compat_enum(qtmods, "QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers",
                            "QtWidgets.QAbstractItemView.NoEditTriggers")

    outer_layout = QtWidgets.QVBoxLayout(dlg)
    tabs = QtWidgets.QTabWidget()
    outer_layout.addWidget(tabs)

    state = {"pharmacophore": None, "library_results": [], "library_path": None,
             "validation": None, "decoy_file_path": None, "active_entries": [],
             "stability_lookup": {}, "validation_used_external_actives": False,
             "ligands_visible": True}

    def scrollable_tab():
        """Each tab gets its own scroll area, so content never overlaps
        regardless of platform font metrics/DPI -- a real bug found on
        macOS during development, where a fixed-height layout caused
        widgets to visibly overlap."""
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        content = QtWidgets.QWidget()
        scroll.setWidget(content)
        lay = QtWidgets.QVBoxLayout(content)
        lay.setContentsMargins(14, 14, 14, 14)
        lay.setSpacing(8)
        return scroll, lay

    # ================================================================
    # TAB 1 — Build Model
    # ================================================================
    tab1_scroll, layout = scrollable_tab()
    tabs.addTab(tab1_scroll, "1. Build Model")

    layout.addWidget(_section_label("Known Actives"))
    input_edit = QtWidgets.QPlainTextEdit()
    input_edit.setPlaceholderText("One per line: name, SMILES")
    input_edit.setMinimumHeight(90)
    input_edit.setMaximumHeight(160)
    layout.addWidget(input_edit)

    btn_row = QtWidgets.QHBoxLayout()
    load_example_btn = QtWidgets.QPushButton("Load MAO-B Example")
    analyze_btn = QtWidgets.QPushButton("Analyze \u2192 Load into PyMOL")
    btn_row.addWidget(load_example_btn)
    btn_row.addWidget(analyze_btn)
    layout.addLayout(btn_row)

    adv_row = QtWidgets.QHBoxLayout()
    adv_row.addWidget(QtWidgets.QLabel("Conformers per molecule:"))
    confs_spin = QtWidgets.QSpinBox()
    confs_spin.setMinimum(10)
    confs_spin.setMaximum(200)
    confs_spin.setValue(DEFAULT_N_CONFS)
    confs_spin.setToolTip(
        "More conformers = better chance of finding the true bioactive shape, "
        "at the cost of more computation time. Default (50) matches common "
        "practice for ligand-based pharmacophore work."
    )
    adv_row.addWidget(confs_spin)
    protonate_check = QtWidgets.QCheckBox("Protonate at physiological pH (~7.4)")
    protonate_check.setChecked(True)
    protonate_check.setToolTip(
        "Protonates aliphatic amines and deprotonates carboxylic acids before "
        "3D embedding. A conservative, explicit rule set — not a full pKa "
        "predictor. See README for exact scope."
    )
    adv_row.addWidget(protonate_check)
    adv_row.addStretch()
    layout.addLayout(adv_row)

    tol_row = QtWidgets.QHBoxLayout()
    tol_row.addWidget(QtWidgets.QLabel("Match tolerance:"))
    tol_spin = QtWidgets.QDoubleSpinBox()
    tol_spin.setMinimum(0.5)
    tol_spin.setMaximum(4.0)
    tol_spin.setSingleStep(0.1)
    tol_spin.setValue(1.8)
    tol_spin.setSuffix(" \u00c5")
    tol_spin.setToolTip(
        "Maximum distance (\u00c5) between a candidate's feature and a consensus centroid to "
        "still count as a match, used everywhere a candidate is scored (single candidate, "
        "batch screening, validation). Tighter = stricter, fewer matches; looser = more "
        "forgiving. 1.8 \u00c5 is a reasonable default, not a rigorously fit value for every case."
    )
    tol_row.addWidget(tol_spin)
    tol_row.addStretch()
    layout.addLayout(tol_row)

    status_label = QtWidgets.QLabel("Ready. (Results will appear in your PyMOL viewer.)")
    status_label.setWordWrap(True)
    layout.addWidget(status_label)

    layout.addWidget(_section_label("Consensus Pharmacophore"))

    thresh_row = QtWidgets.QHBoxLayout()
    thresh_row.addWidget(QtWidgets.QLabel("Support threshold:"))
    thresh_slider = QtWidgets.QSlider(horizontal)
    thresh_slider.setMinimum(0)
    thresh_slider.setMaximum(100)
    thresh_slider.setValue(50)
    thresh_label = QtWidgets.QLabel("50%")
    thresh_row.addWidget(thresh_slider)
    thresh_row.addWidget(thresh_label)
    layout.addLayout(thresh_row)

    eps_row = QtWidgets.QHBoxLayout()
    eps_row.addWidget(QtWidgets.QLabel("Clustering radius (eps):"))
    eps_slider = QtWidgets.QSlider(horizontal)
    eps_slider.setMinimum(8)   # represents 0.8 A (slider works in tenths of an Angstrom)
    eps_slider.setMaximum(25)  # represents 2.5 A
    eps_slider.setValue(15)    # represents 1.5 A, the previous fixed default
    eps_label = QtWidgets.QLabel("1.5 \u00c5")
    eps_slider.setToolTip(
        "How close (in \u00c5) feature points from different actives must be to count as the same "
        "consensus feature. Tighter (lower) for a rigid, closely-matched active set; looser "
        "(higher) if actives are more conformationally flexible or diverse. 1.5 \u00c5 is a common "
        "default (matching e.g. LigandScout/MOE), not a rigorously-derived value for every case."
    )
    eps_row.addWidget(eps_slider)
    eps_row.addWidget(eps_label)
    layout.addLayout(eps_row)

    feature_table = QtWidgets.QTableWidget(0, 5)
    feature_table.setHorizontalHeaderLabels(["Feature", "Support", "Points", "Stability", "Consensus?"])
    feature_table.horizontalHeader().setStretchLastSection(True)
    feature_table.setEditTriggers(no_edit)
    feature_table.setMinimumHeight(140)
    layout.addWidget(feature_table)

    robustness_row = QtWidgets.QHBoxLayout()
    robustness_btn = QtWidgets.QPushButton("Run Leave-One-Out Robustness Check")
    robustness_btn.setEnabled(False)
    robustness_row.addWidget(robustness_btn)
    layout.addLayout(robustness_row)

    robustness_progress = QtWidgets.QProgressBar()
    robustness_progress.setVisible(False)
    layout.addWidget(robustness_progress)

    robustness_note = QtWidgets.QLabel(
        "Stability = fraction of leave-one-out reruns (each excluding one active) in which this "
        "feature still appears. Low stability means the feature depends on one or two specific "
        "compounds, not a pattern genuinely shared across the whole active set."
    )
    robustness_note.setWordWrap(True)
    robustness_note.setStyleSheet("color: #64748B; font-size: 11px;")
    layout.addWidget(robustness_note)

    view_row = QtWidgets.QHBoxLayout()
    zoom_btn = QtWidgets.QPushButton("Zoom to Pharmacophore")
    toggle_ligs_btn = QtWidgets.QPushButton("Toggle Ligands")
    view_row.addWidget(zoom_btn)
    view_row.addWidget(toggle_ligs_btn)
    layout.addLayout(view_row)

    model_io_row = QtWidgets.QHBoxLayout()
    save_model_btn = QtWidgets.QPushButton("Save Model...")
    load_model_btn = QtWidgets.QPushButton("Load Model...")
    compare_models_btn = QtWidgets.QPushButton("Compare Models...")
    save_model_btn.setEnabled(False)
    model_io_row.addWidget(save_model_btn)
    model_io_row.addWidget(load_model_btn)
    model_io_row.addWidget(compare_models_btn)
    layout.addLayout(model_io_row)

    layout.addWidget(_legend_widget())
    layout.addStretch()

    # ================================================================
    # TAB 2 — Screen
    # ================================================================
    tab2_scroll, layout2 = scrollable_tab()
    tabs.addTab(tab2_scroll, "2. Screen")

    layout2.addWidget(_section_label("Screen a Candidate"))
    cand_row = QtWidgets.QHBoxLayout()
    candidate_edit = QtWidgets.QLineEdit()
    candidate_edit.setPlaceholderText("Paste a SMILES to score + load into PyMOL")
    score_btn = QtWidgets.QPushButton("Score")
    cand_row.addWidget(candidate_edit)
    cand_row.addWidget(score_btn)
    layout2.addLayout(cand_row)

    verdict_label = QtWidgets.QLabel("")
    verdict_label.setWordWrap(True)
    layout2.addWidget(verdict_label)

    layout2.addWidget(_section_label("Batch Screen a Library (SDF)"))
    lib_row = QtWidgets.QHBoxLayout()
    lib_path_label = QtWidgets.QLabel("No file selected.")
    lib_path_label.setWordWrap(True)
    choose_lib_btn = QtWidgets.QPushButton("Choose SDF File...")
    lib_row.addWidget(choose_lib_btn)
    layout2.addLayout(lib_row)
    layout2.addWidget(lib_path_label)

    lib_action_row = QtWidgets.QHBoxLayout()
    screen_lib_btn = QtWidgets.QPushButton("Screen Library")
    export_csv_btn = QtWidgets.QPushButton("Export Results (CSV)")
    export_sdf_btn = QtWidgets.QPushButton("Export Aligned Hits (SDF)")
    export_csv_btn.setEnabled(False)
    export_sdf_btn.setEnabled(False)
    lib_action_row.addWidget(screen_lib_btn)
    lib_action_row.addWidget(export_csv_btn)
    lib_action_row.addWidget(export_sdf_btn)
    layout2.addLayout(lib_action_row)

    diverse_check = QtWidgets.QCheckBox("Diverse hits only (remove near-duplicate analogs)")
    diverse_check.setToolTip(
        "Keeps only the best-scoring representative of each cluster of structurally "
        "similar hits (Tanimoto > 0.6), instead of letting minor substituent variants "
        "on one scaffold fill the entire top of the ranked list."
    )
    layout2.addWidget(diverse_check)

    mw_filter_row = QtWidgets.QHBoxLayout()
    mw_filter_check = QtWidgets.QCheckBox("Restrict to MW range of training actives \u00b1")
    mw_tolerance_spin = QtWidgets.QSpinBox()
    mw_tolerance_spin.setMinimum(0)
    mw_tolerance_spin.setMaximum(500)
    mw_tolerance_spin.setValue(50)
    mw_tolerance_spin.setSuffix(" Da")
    mw_filter_check.setToolTip(
        "A candidate much larger than every training active likely cannot physically occupy "
        "the same binding pocket, regardless of pharmacophore fit \u2014 found, in real use, to be "
        "a more diagnostic red flag than chemical-plausibility similarity alone."
    )
    mw_filter_row.addWidget(mw_filter_check)
    mw_filter_row.addWidget(mw_tolerance_spin)
    mw_filter_row.addStretch()
    layout2.addLayout(mw_filter_row)

    prefilter_check = QtWidgets.QCheckBox("Fast Screen (skip full scoring for obvious non-matches)")
    prefilter_check.setToolTip(
        "For large libraries: a quick, alignment-free geometric check (~6x faster) rejects "
        "compounds that clearly cannot match before running full O3A-based scoring on the "
        "rest. Deliberately lenient \u2014 verified on real data to be noticeably more permissive "
        "than full scoring, so it only saves time on the clearly-hopeless majority, never "
        "replaces the full score for anything that survives it. Off by default; most useful "
        "for libraries of several thousand compounds or more."
    )
    layout2.addWidget(prefilter_check)

    ensemble_check = QtWidgets.QCheckBox("Ensemble tolerance matching (nearest raw point, not just centroid)")
    ensemble_check.setToolTip(
        "Instead of measuring a candidate feature's distance to the single averaged consensus "
        "position, measures distance to whichever individual contributing point (from the real "
        "aligned actives that formed that cluster) is closest. Tolerates real positional spread "
        "across a flexible or diverse active set rather than collapsing it to one fixed point. "
        "Off by default (uses the centroid, as before) for consistency with existing results; "
        "only available for features built with this version or later (falls back to centroid "
        "automatically for older saved models)."
    )
    layout2.addWidget(ensemble_check)

    lib_progress = QtWidgets.QProgressBar()
    lib_progress.setVisible(False)
    layout2.addWidget(lib_progress)

    lib_status_label = QtWidgets.QLabel("")
    lib_status_label.setWordWrap(True)
    layout2.addWidget(lib_status_label)

    lib_table = QtWidgets.QTableWidget(0, 9)
    lib_table.setHorizontalHeaderLabels(["Name", "Matched", "Fraction", "Mean Dist (\u00c5)", "Align RMSD (\u00c5)", "Plausibility", "Shape", "Scaffold", "Verdict"])
    lib_table.horizontalHeader().setStretchLastSection(True)
    lib_table.setEditTriggers(no_edit)
    lib_table.setMinimumHeight(160)
    lib_table.setToolTip(
        "Click a row to load that compound into PyMOL alongside the pharmacophore. "
        "Align RMSD = how well the whole candidate fit the reference during alignment (low = "
        "credible fit; a high fraction_matched with a high RMSD suggests the matched features "
        "landed in tolerance somewhat by chance). "
        "Plausibility = Tanimoto similarity to the nearest known active (a chemical-plausibility "
        "heuristic, not formal QSAR applicability-domain) \u2014 a candidate can match the 3D geometry "
        "by coincidence while being chemically unlike anything the model was actually built from.\n\n"
        "Shape = % of this candidate's own atoms that fall within the training actives' combined "
        "shape. A LIGAND-shape check only \u2014 not a receptor steric-clash check, since this tool has "
        "no receptor structure. Low shape consistency flags real substituent bulk the training "
        "actives never had, something the pharmacophore's isolated feature points can miss entirely.\n\n"
        "Scaffold = Bemis-Murcko comparison to the training actives' own core ring systems: Same "
        "(e.g. a prodrug or minor variant), Analog (related), or Novel (a genuine scaffold-hop).\n\n"
        "Verdict combines match + plausibility: a high match + low plausibility is labeled a "
        "'Scaffold-Hop' risk, not discarded \u2014 pharmacophore modeling exists specifically to find "
        "real hits with different scaffolds from the training actives. Low plausibility flags it as "
        "needing extra scrutiny (size, rotatable bonds, docking), never as an automatic rejection."
    )
    layout2.addWidget(lib_table)
    layout2.addStretch()

    # ================================================================
    # TAB 3 — Validate
    # ================================================================
    tab3_scroll, layout3 = scrollable_tab()
    tabs.addTab(tab3_scroll, "3. Validate")

    layout3.addWidget(_section_label("Validate This Model (EF / ROC-AUC / GH)"))
    val_note = QtWidgets.QLabel(
        "Tests whether this pharmacophore actually discriminates the known actives "
        "(Tab 1) from a decoy set you supply \u2014 a model is a hypothesis to be checked, "
        "not a final answer. Bring your own decoys: e.g. LIDeB Tools, DeepCoy, or the "
        "DUD-E decoy server all generate well-validated decoy sets for a given active "
        "set, and do that job better than an automatic generator built into this plugin "
        "could."
    )
    val_note.setWordWrap(True)
    layout3.addWidget(val_note)

    loocv_check = QtWidgets.QCheckBox("Use leave-one-out cross-validation for actives (recommended)")
    loocv_check.setChecked(True)
    loocv_check.setToolTip(
        "Scores each active only against a model rebuilt WITHOUT it, avoiding the "
        "methodological flaw of testing on the same molecules used to build the model. "
        "Automatically not used when external test actives (below) are supplied, since "
        "leave-one-out only makes sense for actives that were part of training."
    )
    layout3.addWidget(loocv_check)

    layout3.addWidget(_section_label("External test actives (optional)"))
    ext_actives_note = QtWidgets.QLabel(
        "Leave empty to validate the Tab 1 training actives via leave-one-out (the default, "
        "and correct choice whenever your actives were used to build this model). Paste "
        "genuinely different, held-out actives here \u2014 real compounds NOT used to build "
        "this model \u2014 to validate against them directly instead. When this box has content, "
        "leave-one-out is skipped automatically (there's nothing to leave out; these actives "
        "were never part of training)."
    )
    ext_actives_note.setWordWrap(True)
    ext_actives_note.setStyleSheet("color: #64748B; font-size: 11px;")
    layout3.addWidget(ext_actives_note)
    ext_actives_edit = QtWidgets.QPlainTextEdit()
    ext_actives_edit.setPlaceholderText("One per line: name, SMILES (only used if non-empty)")
    ext_actives_edit.setMinimumHeight(60)
    ext_actives_edit.setMaximumHeight(100)
    layout3.addWidget(ext_actives_edit)

    layout3.addWidget(_section_label("Decoy file"))
    decoy_file_row = QtWidgets.QHBoxLayout()
    decoy_file_label = QtWidgets.QLabel("No decoy file selected.")
    decoy_file_label.setWordWrap(True)
    choose_decoy_file_btn = QtWidgets.QPushButton("Choose Decoy File...")
    decoy_file_row.addWidget(choose_decoy_file_btn)
    layout3.addLayout(decoy_file_row)
    layout3.addWidget(decoy_file_label)

    val_btn_row = QtWidgets.QHBoxLayout()
    validate_btn = QtWidgets.QPushButton("Run Validation")
    val_btn_row.addWidget(validate_btn)
    layout3.addLayout(val_btn_row)

    val_progress = QtWidgets.QProgressBar()
    val_progress.setVisible(False)
    layout3.addWidget(val_progress)

    val_status_label = QtWidgets.QLabel("")
    val_status_label.setWordWrap(True)
    layout3.addWidget(val_status_label)

    val_results_label = QtWidgets.QLabel("")
    val_results_label.setWordWrap(True)
    val_results_label.setStyleSheet("font-family: monospace; font-size: 11px;")
    layout3.addWidget(val_results_label)

    layout3.addWidget(_section_label("ROC Curve"))
    roc_widget = _make_roc_widget(QtWidgets, QtGui, QtCore)
    layout3.addWidget(roc_widget)

    export_roc_row = QtWidgets.QHBoxLayout()
    export_roc_btn = QtWidgets.QPushButton("Export Plot (PNG)")
    export_report_btn = QtWidgets.QPushButton("Export Full Report (HTML)")
    export_roc_btn.setEnabled(False)
    export_report_btn.setEnabled(False)
    export_roc_row.addWidget(export_roc_btn)
    export_roc_row.addWidget(export_report_btn)
    export_roc_row.addStretch()
    layout3.addLayout(export_roc_row)
    layout3.addStretch()

    # ---------------------------------------------------- actions

    def do_load_example():
        input_edit.setPlainText(EXAMPLE_TEXT)

    def do_analyze():
        text = input_edit.toPlainText().strip()
        if not text:
            status_label.setText("Please enter at least 2 molecules.")
            return
        entries = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if "," in line:
                name, smi = line.split(",", 1)
            else:
                name, smi = f"mol_{len(entries)+1}", line
            entries.append((name.strip(), smi.strip()))

        status_label.setText("Generating conformers and aligning...")
        QtWidgets.QApplication.processEvents()

        n_confs = confs_spin.value()
        do_protonate = protonate_check.isChecked()
        ligands = [Ligand(name, smi, n_confs=n_confs, protonate=do_protonate) for name, smi in entries]
        valid = [l for l in ligands if l.is_valid]
        invalid = [l for l in ligands if not l.is_valid]
        if len(valid) < 2:
            status_label.setText(f"Only {len(valid)} valid molecule(s) parsed \u2014 need at least 2.")
            return

        pharm = LigandBasedPharmacophore(valid)
        pharm.align()
        pharm.extract_consensus(eps=eps_slider.value() / 10.0, min_support_frac=0.0)
        state["pharmacophore"] = pharm
        state["active_entries"] = [(l.name, l.smiles) for l in valid]
        state["stability_lookup"] = {}

        msg = f"Analyzed {len(valid)} ligand(s). Loaded into PyMOL viewer."
        if invalid:
            msg += f" ({len(invalid)} could not be parsed/embedded.)"
        status_label.setText(msg)

        _load_into_pymol(cmd, pharm)
        state["ligands_visible"] = True  # freshly-loaded objects are enabled by default
        _refresh_feature_table(feature_table, pharm, thresh_slider.value() / 100.0)
        _refresh_pharmacophore_spheres(cmd, pharm, thresh_slider.value() / 100.0)
        save_model_btn.setEnabled(True)
        robustness_btn.setEnabled(True)

    def do_threshold_changed():
        pct = thresh_slider.value()
        thresh_label.setText(f"{pct}%")
        pharm = state["pharmacophore"]
        if pharm is not None:
            _refresh_feature_table(feature_table, pharm, pct / 100.0, stability_lookup=state["stability_lookup"])
            _refresh_pharmacophore_spheres(cmd, pharm, pct / 100.0)

    def do_eps_changed():
        eps_val = eps_slider.value() / 10.0
        eps_label.setText(f"{eps_val:.1f} \u00c5")
        pharm = state["pharmacophore"]
        if pharm is None or not getattr(pharm, "ligands", None):
            return  # re-clustering needs the original aligned Ligand objects (not available
                     # after a Load Model, which only stores already-extracted consensus points)
        pharm.extract_consensus(eps=eps_val, min_support_frac=0.0)
        state["stability_lookup"] = {}  # stale after re-clustering with a different eps
        pct = thresh_slider.value() / 100.0
        _refresh_feature_table(feature_table, pharm, pct, stability_lookup={})
        _refresh_pharmacophore_spheres(cmd, pharm, pct)

    def do_zoom():
        cmd.zoom("ligand_* or pharm_*", buffer=2.0)

    def do_toggle_ligands():
        # cmd.toggle()'s first argument is a representation name (e.g. "sticks"),
        # not an object/selection pattern -- cmd.toggle("ligand_*") was raising
        # "unknown representation" and silently doing nothing useful. PyMOL has
        # no built-in "flip visibility" for a selection, so state is tracked
        # explicitly and cmd.enable/cmd.disable (which do support wildcards) are
        # used instead.
        state["ligands_visible"] = not state["ligands_visible"]
        if state["ligands_visible"]:
            cmd.enable("ligand_*")
        else:
            cmd.disable("ligand_*")

    def do_run_robustness():
        import numpy as np
        pharm = state["pharmacophore"]
        entries = state["active_entries"]
        if pharm is None or len(entries) < 3:
            status_label.setText("Need at least 3 known actives (Analyze first) to run a robustness check.")
            return

        robustness_btn.setEnabled(False)
        robustness_progress.setVisible(True)
        robustness_progress.setValue(0)
        status_label.setText("Running leave-one-out robustness check (rebuilding the model N times)...")
        QtWidgets.QApplication.processEvents()

        def progress_cb(i, n_total):
            if n_total > 0:
                robustness_progress.setMaximum(n_total)
                robustness_progress.setValue(i + 1)
            QtWidgets.QApplication.processEvents()

        result = leave_one_out_stability(
            entries, n_confs=confs_spin.value(), protonate=protonate_check.isChecked(),
            progress_callback=progress_cb,
        )
        robustness_progress.setVisible(False)
        robustness_btn.setEnabled(True)

        if "error" in result:
            status_label.setText(result["error"])
            return

        # Match each leave-one-out "full_points" result back to the currently
        # displayed pharm.points by family + spatial proximity (they're separately
        # computed objects, even though the underlying computation is deterministic).
        lookup = {}
        for loo_p, stab in zip(result["full_points"], result["stability"]):
            best, best_dist = None, 1e9
            for p in pharm.points:
                if p.family != loo_p.family:
                    continue
                dist = float(np.linalg.norm(np.array(p.centroid) - np.array(loo_p.centroid)))
                if dist < best_dist:
                    best, best_dist = p, dist
            if best is not None and best_dist <= 2.0:
                lookup[id(best)] = stab
        state["stability_lookup"] = lookup

        # Contextual guidance: specifically check the "core" (high-support) features --
        # the ones that SHOULD be robust, since nearly every active contributed to them.
        # If even those show low stability under leave-one-out, that's real evidence
        # the active set may not share one consistent binding mode (see README on why
        # automatic multi-hypothesis *detection* was tried and abandoned -- this check,
        # prompting the user to decide, is the safer alternative).
        core_points = [p for p in pharm.points if p.support_frac >= 0.8]
        unstable_core = [p for p in core_points if lookup.get(id(p), 1.0) < 0.6]
        guidance = ""
        if unstable_core:
            fams = ", ".join(sorted(set(p.family for p in unstable_core)))
            guidance = (
                f"\n\nRobustness check indicates instability: {len(unstable_core)} high-support "
                f"feature(s) ({fams}) dropped below 60% stability under leave-one-out testing. "
                f"Consider splitting your actives into separate, more chemically consistent "
                f"groups, or manually selecting a different alignment reference (Tab 1)."
            )

        status_label.setText(
            f"Robustness check complete ({result['n_runs']} leave-one-out reruns). "
            f"See the Stability column below.{guidance}"
        )
        _refresh_feature_table(feature_table, pharm, thresh_slider.value() / 100.0, stability_lookup=lookup)

    def do_save_model():
        pharm = state["pharmacophore"]
        if pharm is None:
            status_label.setText("Run Analyze first \u2014 nothing to save yet.")
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            dlg, "Save pharmacophore model", "pharmacophore_model.json", "JSON files (*.json)"
        )
        if not path:
            return
        metadata = {
            "n_ligands": len(pharm.ligands),
            "ligand_names": [l.name for l in pharm.ligands],
            "ligand_smiles": [l.smiles for l in pharm.ligands],
            "conformers_per_molecule": confs_spin.value(),
            "protonated": protonate_check.isChecked(),
        }
        try:
            export_pharmacophore_json(pharm.points, path, metadata=metadata, ligands=pharm.ligands)
            status_label.setText(f"Model saved to {path}")
        except Exception as e:
            status_label.setText(f"Failed to save model: {e}")

    def do_load_model():
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            dlg, "Load a saved pharmacophore model", "", "JSON files (*.json);;All files (*)"
        )
        if not path:
            return
        try:
            loaded = load_pharmacophore_json(path)
        except Exception as e:
            status_label.setText(f"Failed to load model: {e}")
            return

        state["pharmacophore"] = loaded
        n_ligs = loaded.metadata.get("n_ligands", len(loaded.ligand_mols))
        msg = f"Loaded model: {len(loaded.points)} feature cluster(s) from {n_ligs} original ligand(s)."
        if not loaded.ligand_mols:
            msg += " (No saved ligand structures in this file \u2014 spheres only.)"
        status_label.setText(msg)

        _load_into_pymol(cmd, loaded)
        state["ligands_visible"] = True  # freshly-loaded objects are enabled by default
        _refresh_feature_table(feature_table, loaded, thresh_slider.value() / 100.0)
        _refresh_pharmacophore_spheres(cmd, loaded, thresh_slider.value() / 100.0)
        save_model_btn.setEnabled(True)

    def do_compare_models():
        sub = QtWidgets.QDialog(dlg)
        sub.setWindowTitle("Compare Two Pharmacophore Models")
        sub.resize(560, 480)
        sub_layout = QtWidgets.QVBoxLayout(sub)

        sub_state = {"path_a": None, "path_b": None}

        row_a = QtWidgets.QHBoxLayout()
        label_a = QtWidgets.QLabel("Model A: (none selected)")
        label_a.setWordWrap(True)
        btn_a = QtWidgets.QPushButton("Choose Model A...")
        row_a.addWidget(btn_a)
        sub_layout.addLayout(row_a)
        sub_layout.addWidget(label_a)

        row_b = QtWidgets.QHBoxLayout()
        label_b = QtWidgets.QLabel("Model B: (none selected)")
        label_b.setWordWrap(True)
        btn_b = QtWidgets.QPushButton("Choose Model B...")
        row_b.addWidget(btn_b)
        sub_layout.addLayout(row_b)
        sub_layout.addWidget(label_b)

        compare_btn = QtWidgets.QPushButton("Compare")
        compare_btn.setEnabled(False)
        sub_layout.addWidget(compare_btn)

        sub_status = QtWidgets.QLabel("")
        sub_status.setWordWrap(True)
        sub_layout.addWidget(sub_status)

        compare_table = QtWidgets.QTableWidget(0, 5)
        compare_table.setHorizontalHeaderLabels(["Feature", "In A?", "Support A", "In B?", "Support B"])
        compare_table.horizontalHeader().setStretchLastSection(True)
        compare_table.setEditTriggers(no_edit)
        sub_layout.addWidget(compare_table)

        def check_ready():
            compare_btn.setEnabled(bool(sub_state["path_a"] and sub_state["path_b"]))

        def choose_a():
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                sub, "Choose Model A", "", "JSON files (*.json);;All files (*)"
            )
            if path:
                sub_state["path_a"] = path
                label_a.setText(f"Model A: {path}")
                check_ready()

        def choose_b():
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                sub, "Choose Model B", "", "JSON files (*.json);;All files (*)"
            )
            if path:
                sub_state["path_b"] = path
                label_b.setText(f"Model B: {path}")
                check_ready()

        def do_compare():
            try:
                model_a = load_pharmacophore_json(sub_state["path_a"])
                model_b = load_pharmacophore_json(sub_state["path_b"])
            except Exception as e:
                sub_status.setText(f"Failed to load one of the models: {e}")
                return
            rows = compare_pharmacophore_models(model_a, {}, model_b, {})
            n_shared = sum(1 for r in rows if r["shared"])
            sub_status.setText(
                f"{len(rows)} distinct feature(s) total across both models \u2014 "
                f"{n_shared} shared, {len(rows) - n_shared} unique to one model."
            )
            compare_table.setRowCount(len(rows))
            for i, r in enumerate(rows):
                compare_table.setItem(i, 0, QtWidgets.QTableWidgetItem(r["family"]))
                compare_table.setItem(i, 1, QtWidgets.QTableWidgetItem("\u2713" if r["in_a"] else "\u2014"))
                compare_table.setItem(i, 2, QtWidgets.QTableWidgetItem(
                    f"{r['support_a']*100:.0f}%" if r["support_a"] is not None else "\u2014"))
                compare_table.setItem(i, 3, QtWidgets.QTableWidgetItem("\u2713" if r["in_b"] else "\u2014"))
                compare_table.setItem(i, 4, QtWidgets.QTableWidgetItem(
                    f"{r['support_b']*100:.0f}%" if r["support_b"] is not None else "\u2014"))
                if r["shared"]:
                    for col in range(5):
                        item = compare_table.item(i, col)
                        item.setForeground(QtGui.QColor("#0F6E56"))

        btn_a.clicked.connect(choose_a)
        btn_b.clicked.connect(choose_b)
        compare_btn.clicked.connect(do_compare)

        sub.exec()

    def do_score_candidate():
        pharm = state["pharmacophore"]
        if pharm is None:
            verdict_label.setText("Run Analyze first.")
            return
        smi = candidate_edit.text().strip()
        if not smi:
            return
        pct = thresh_slider.value() / 100.0
        consensus_pts = pharm.consensus_at_threshold(pct)
        if not consensus_pts:
            verdict_label.setText("No features meet the current threshold.")
            return
        ref_mol, ref_cid = pharm.get_alignment_reference()
        n_matched, frac, mean_d, mol, cid, match_details, align_rmsd, shape_consistency = score_molecule(
            smi, consensus_pts, reference_mol=ref_mol, reference_conf_id=ref_cid,
            tol=tol_spin.value(), n_confs=confs_spin.value(), protonate=protonate_check.isChecked(),
            match_mode="ensemble" if ensemble_check.isChecked() else "centroid",
        )
        if mol is None:
            verdict_label.setText("Could not parse/embed that SMILES.")
            return

        active_smiles = _active_smiles_from(pharm)
        plausibility = None
        if active_smiles:
            max_sim, mean_sim = applicability_domain_similarity(smi, active_smiles)
            plausibility = max_sim  # None if the candidate SMILES itself failed to parse here too

        verdict_text, verdict_color, verdict_explanation = classify_verdict(frac, plausibility)
        verdict_label.setStyleSheet(f"color: {verdict_color}; font-weight: 600;")

        plaus_txt = f" | Chemical plausibility: {plausibility:.2f}" if plausibility is not None else ""
        rmsd_txt = f" | Alignment RMSD: {align_rmsd:.2f} \u00c5" if align_rmsd is not None else ""

        verdict_label.setText(
            f"[{verdict_text}]  Matched {n_matched}/{len(consensus_pts)} consensus features "
            f"({frac*100:.0f}%), mean fit distance {mean_d:.2f} \u00c5 \u2014 loaded into PyMOL as 'candidate'."
            f"{rmsd_txt}{plaus_txt}\n{verdict_explanation}"
        )
        _load_candidate_into_pymol(cmd, mol, cid, match_details=match_details)

    def _active_smiles_from(pharm):
        """Works for both a freshly-built LigandBasedPharmacophore (pharm.ligands)
        and a reloaded LoadedPharmacophore (SMILES saved in pharm.metadata)."""
        if getattr(pharm, "ligands", None):
            return [l.smiles for l in pharm.ligands]
        return pharm.metadata.get("ligand_smiles", [])

    def _active_entries_from(pharm):
        """Like _active_smiles_from, but returns (name, smiles) pairs -- needed
        for leave-one-out validation, which must rebuild the model from
        scratch and needs names for its per-active result records. Works for
        both a fresh model and one reloaded from a saved file (as long as
        that file's metadata included ligand_names alongside ligand_smiles,
        which export_pharmacophore_json always writes)."""
        if getattr(pharm, "ligands", None):
            return [(l.name, l.smiles) for l in pharm.ligands]
        names = pharm.metadata.get("ligand_names", [])
        smiles = pharm.metadata.get("ligand_smiles", [])
        if len(names) == len(smiles):
            return list(zip(names, smiles))
        return [(f"active_{i+1}", s) for i, s in enumerate(smiles)]

    def do_choose_decoy_file():
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            dlg, "Choose a decoy file", "", "SDF or SMILES files (*.sdf *.smi *.txt);;All files (*)"
        )
        if path:
            state["decoy_file_path"] = path
            decoy_file_label.setText(path)

    def _get_decoys(active_smiles):
        """Returns (decoy_smiles, source_description). Decoys always come
        from a user-supplied file — see the note in the Validate tab for why
        automatic decoy generation was deliberately removed in favor of
        purpose-built external tools (LIDeB Tools, DeepCoy, the DUD-E decoy
        server, etc.), which do this job better than a generic in-plugin
        generator could."""
        path = state["decoy_file_path"]
        if not path:
            return [], "No decoy file selected."
        return load_decoys_from_file(path)

    def _parse_actives_text(text):
        """Same 'name, SMILES per line' format as Tab 1's input box."""
        out = []
        for line in text.strip().splitlines():
            line = line.strip()
            if not line:
                continue
            if "," in line:
                name, smi = line.split(",", 1)
            else:
                name, smi = f"active_{len(out)+1}", line
            out.append((name.strip(), smi.strip()))
        return out

    def do_validate():
        pharm = state["pharmacophore"]
        if pharm is None:
            val_status_label.setText("Run Analyze (or Load Model) first, on Tab 1.")
            return

        ext_text = ext_actives_edit.toPlainText().strip()
        using_external = bool(ext_text)

        if using_external:
            ext_entries = _parse_actives_text(ext_text)
            active_smiles = [smi for _, smi in ext_entries]
            use_loocv = False  # leave-one-out doesn't apply -- these were never part of training
        else:
            active_smiles = _active_smiles_from(pharm)
            use_loocv = loocv_check.isChecked()

        if len(active_smiles) < 2:
            val_status_label.setText(
                "Need at least 2 active SMILES to validate against \u2014 "
                + ("check the external test actives box." if using_external
                   else "none were available from the current model.")
            )
            return

        pct = thresh_slider.value() / 100.0

        if not use_loocv:
            consensus_pts = pharm.consensus_at_threshold(pct)
            if not consensus_pts:
                val_status_label.setText("No features meet the current support threshold (Tab 1).")
                return

        validate_btn.setEnabled(False)
        val_progress.setVisible(True)
        val_progress.setValue(0)
        val_results_label.setText("")
        roc_widget.clear_curve()
        QtWidgets.QApplication.processEvents()

        decoys, decoy_source_desc = _get_decoys(active_smiles)
        if len(decoys) < 2:
            val_status_label.setText(
                f"Only {len(decoys)} decoy(s) available ({decoy_source_desc}) \u2014 not enough to validate against."
            )
            validate_btn.setEnabled(True)
            val_progress.setVisible(False)
            return

        val_status_label.setText(
            f"Scoring {len(active_smiles)} actives + {len(decoys)} decoys"
            f"{' (leave-one-out: rebuilding the model once per active)' if use_loocv else ''}..."
        )
        QtWidgets.QApplication.processEvents()

        def progress_cb(i, n_total):
            if n_total > 0:
                val_progress.setMaximum(n_total)
                val_progress.setValue(i + 1)
            QtWidgets.QApplication.processEvents()

        if use_loocv:
            active_entries = _active_entries_from(pharm)
            result = validate_pharmacophore_loocv(
                active_entries, decoys,
                eps=eps_slider.value() / 10.0, support_thresh=pct, tol=tol_spin.value(),
                n_confs=confs_spin.value(), protonate=protonate_check.isChecked(),
                progress_callback=progress_cb,
            )
        else:
            ref_mol, ref_cid = pharm.get_alignment_reference()
            result = validate_pharmacophore(
                consensus_pts, active_smiles, decoys,
                reference_mol=ref_mol, reference_conf_id=ref_cid, tol=tol_spin.value(),
                n_confs=confs_spin.value(), protonate=protonate_check.isChecked(),
                progress_callback=progress_cb,
            )
        state["validation"] = result
        state["validation_used_external_actives"] = using_external
        val_progress.setVisible(False)
        validate_btn.setEnabled(True)

        if "error" in result:
            val_status_label.setText(result["error"])
            return

        # Explicit, unambiguous wording: validation-set size vs. decoy-source size are
        # two different numbers and were previously easy to misread as summing together.
        if using_external:
            mode_note = "External test actives (genuinely held-out \u2014 not used to build this model)."
        elif result.get("loocv_actives"):
            mode_note = (
                f"Leave-one-out ({result.get('n_loocv_runs', '?')} rebuilds \u2014 each active scored "
                f"only against a model that never saw it)."
            )
        else:
            mode_note = "In-sample (actives were also used to build this model \u2014 see tooltip)."
        val_status_label.setText(
            f"Validation set: {result['n_actives']} known active(s) vs. {result['n_decoys']} decoy(s) "
            f"(decoys from: {decoy_source_desc}). {mode_note}"
        )
        auc_txt = f"{result['auc']:.2f}" if result['auc'] == result['auc'] else "n/a"  # NaN check
        if result['gh'] == result['gh']:
            gh_txt = f"{result['gh']:.3f}"
            gh_detail = (f"(hit = \u2265{int(result['hit_frac_cutoff']*100)}% features matched; "
                         f"{result['Ha']}/{result['Ht']} hits were true actives)")
        else:
            gh_txt = "undefined"
            gh_detail = "(no compound cleared the hit threshold at this support level \u2014 try lowering it)"
        p_txt = ""
        if result.get("p_value") == result.get("p_value"):  # not NaN
            p_val = result["p_value"]
            p_res = result.get("p_resolution")
            sig_note = "significant at p<0.05" if p_val < 0.05 else "not significant at p<0.05"
            res_note = f", smallest possible p-value given dataset size \u2248 {p_res:.4f}" if p_res else ""
            p_txt = f"\nPermutation p-value:  {p_val:.4f}  ({sig_note}{res_note})"
        val_results_label.setText(
            f"ROC-AUC:  {auc_txt}   (0.5 = random, 1.0 = perfect)\n"
            f"EF @ top {int(result['top_pct']*100)}%:  {result['ef']:.2f}x   (vs. random chance)\n"
            f"G\u00fcner-Henry (GH):  {gh_txt}   {gh_detail}"
            f"{p_txt}"
        )
        roc_widget.set_curve(result["fpr"], result["tpr"], result["auc"],
                              n_actives=result["n_actives"], n_decoys=result["n_decoys"])
        export_roc_btn.setEnabled(True)
        export_report_btn.setEnabled(True)

    def do_export_roc():
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            dlg, "Export ROC curve", "roc_curve.png", "PNG files (*.png)"
        )
        if not path:
            return
        ok = roc_widget.export_png(path, width=1200, height=800)
        val_status_label.setText(
            f"ROC curve exported to {path}" if ok else "Failed to export ROC curve image."
        )

    def do_export_report():
        result = state["validation"]
        pharm = state["pharmacophore"]
        if result is None or "error" in (result or {}) or pharm is None:
            val_status_label.setText("Run a successful validation first.")
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            dlg, "Export validation report", "pharmacophore_validation_report.html", "HTML files (*.html)"
        )
        if not path:
            return

        import base64, tempfile, datetime, html as _html

        tmp_png = tempfile.mktemp(suffix=".png")
        roc_widget.export_png(tmp_png, width=1100, height=750)
        with open(tmp_png, "rb") as f:
            img_b64 = base64.b64encode(f.read()).decode("ascii")
        os.unlink(tmp_png)

        active_smiles = _active_smiles_from(pharm)
        if getattr(pharm, "ligands", None):
            active_names = [l.name for l in pharm.ligands]
        else:
            active_names = pharm.metadata.get("ligand_names", [f"active_{i+1}" for i in range(len(active_smiles))])
        pct = thresh_slider.value()
        auc_txt = f"{result['auc']:.3f}" if result['auc'] == result['auc'] else "n/a"
        if result['gh'] == result['gh']:
            gh_txt = f"{result['gh']:.3f}"
            gh_row_detail = (f"(hit = \u2265{int(result['hit_frac_cutoff']*100)}% features matched; "
                             f"{result['Ha']}/{result['Ht']} hits were true actives)")
        else:
            gh_txt = "undefined"
            gh_row_detail = "(no compound cleared the hit threshold at this support level)"

        rows_html = ""
        for r in sorted(result["records"], key=lambda x: (-x["label"], -x["fraction_matched"])):
            label_txt = "Active" if r["label"] == 1 else "Decoy"
            rows_html += (
                f"<tr><td>{_html.escape(label_txt)}</td>"
                f"<td style='font-family:monospace;font-size:12px;'>{_html.escape(r['smiles'])}</td>"
                f"<td>{r['n_matched']}/{r['n_total']}</td>"
                f"<td>{r['fraction_matched']*100:.0f}%</td>"
                f"<td>{r['mean_dist']:.2f}</td></tr>"
            )

        if result.get("p_value") == result.get("p_value"):
            p_value_txt = f"{result['p_value']:.4f}"
            if result.get("p_resolution"):
                resolution_note = " (smallest possible value given dataset size \u2248 %.4f)" % result["p_resolution"]
            else:
                resolution_note = ""
            p_value_display = p_value_txt + resolution_note
        else:
            p_value_display = "n/a"

        if state.get("validation_used_external_actives"):
            validation_mode_display = (
                "External test actives \u2014 genuinely held-out, not used to build this model "
                "(no leave-one-out needed; no leakage concern)"
            )
        elif result.get("loocv_actives"):
            validation_mode_display = (
                "Leave-one-out cross-validation (%d rebuilds) \u2014 each active scored only against "
                "a model that never saw it during training" % result.get("n_loocv_runs", 0)
            )
        else:
            validation_mode_display = (
                "In-sample \u2014 actives were also used to build this model "
                "(see documentation for why leave-one-out is recommended)"
            )

        report_html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>PhArMol Validation Report</title>
<style>
  body {{ font-family: Arial, sans-serif; color: #0F172A; max-width: 900px; margin: 40px auto; padding: 0 20px; }}
  h1 {{ font-size: 22px; border-bottom: 2px solid #0891B2; padding-bottom: 8px; }}
  h2 {{ font-size: 16px; color: #0891B2; margin-top: 32px; }}
  table {{ border-collapse: collapse; width: 100%; margin-top: 8px; font-size: 13px; }}
  th, td {{ border: 1px solid #E2E8F0; padding: 6px 10px; text-align: left; }}
  th {{ background: #F1F5F9; }}
  .meta {{ color: #475569; font-size: 13px; }}
  .metrics td:first-child {{ font-weight: 600; width: 220px; }}
  img {{ max-width: 100%; border: 1px solid #E2E8F0; margin-top: 8px; }}
</style></head>
<body>
<h1>PhArMol Validation Report</h1>
<p class="meta">Generated {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')} \u2014
PhArMol (PyMOL plugin)</p>

<h2>Model</h2>
<table>
<tr><td style="font-weight:600;">Known actives</td><td>{len(active_smiles)}: {_html.escape(', '.join(active_names))}</td></tr>
<tr><td style="font-weight:600;">Conformers per molecule</td><td>{confs_spin.value()}</td></tr>
<tr><td style="font-weight:600;">Protonation applied</td><td>{"Yes" if protonate_check.isChecked() else "No"}</td></tr>
<tr><td style="font-weight:600;">Support threshold used</td><td>{pct}%</td></tr>
<tr><td style="font-weight:600;">Consensus features at this threshold</td><td>{len(pharm.consensus_at_threshold(pct/100.0))}</td></tr>
</table>

<h2>Validation Set</h2>
<table>
<tr><td style="font-weight:600;">Known actives</td><td>{result['n_actives']}</td></tr>
<tr><td style="font-weight:600;">Decoys</td><td>{result['n_decoys']}</td></tr>
<tr><td style="font-weight:600;">Validation mode</td><td>{validation_mode_display}</td></tr>
</table>

<h2>Metrics</h2>
<table class="metrics">
<tr><td>ROC-AUC</td><td>{auc_txt} (0.5 = random, 1.0 = perfect)</td></tr>
<tr><td>Enrichment Factor @ top {int(result['top_pct']*100)}%</td><td>{result['ef']:.2f}\u00d7 vs. random chance</td></tr>
<tr><td>G\u00fcner-Henry (GH) score</td><td>{gh_txt} {gh_row_detail}</td></tr>
<tr><td>Permutation test p-value</td><td>{p_value_display}</td></tr>
</table>

<h2>ROC Curve</h2>
<img src="data:image/png;base64,{img_b64}" alt="ROC Curve"/>

<h2>Per-Compound Results</h2>
<table>
<tr><th>Label</th><th>SMILES</th><th>Matched</th><th>Fraction</th><th>Mean Dist (\u00c5)</th></tr>
{rows_html}
</table>

</body></html>"""

        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(report_html)
            val_status_label.setText(f"Full validation report exported to {path}")
        except Exception as e:
            val_status_label.setText(f"Failed to export report: {e}")

    def do_choose_library():
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            dlg, "Choose a compound library (SDF)", "", "SDF files (*.sdf);;All files (*)"
        )
        if path:
            state["library_path"] = path
            lib_path_label.setText(path)

    def do_screen_library():
        pharm = state["pharmacophore"]
        if pharm is None:
            lib_status_label.setText("Run Analyze first (Tab 1) to build a consensus pharmacophore.")
            return
        path = state["library_path"]
        if not path:
            lib_status_label.setText("Choose an SDF file first.")
            return

        pct = thresh_slider.value() / 100.0
        consensus_pts = pharm.consensus_at_threshold(pct)
        if not consensus_pts:
            lib_status_label.setText("No features meet the current support threshold (Tab 1).")
            return

        lib_progress.setVisible(True)
        lib_progress.setValue(0)
        lib_status_label.setText("Screening library...")
        screen_lib_btn.setEnabled(False)
        QtWidgets.QApplication.processEvents()

        def progress_cb(i, n_total):
            if n_total > 0:
                lib_progress.setMaximum(n_total)
                lib_progress.setValue(i + 1)
            QtWidgets.QApplication.processEvents()

        ref_mol, ref_cid = pharm.get_alignment_reference()
        active_smiles_for_ad = _active_smiles_from(pharm)
        shape_envelope = pharm.get_shape_envelope()
        results, n_skipped, n_prefiltered = screen_library(
            path, consensus_pts,
            reference_mol=ref_mol, reference_conf_id=ref_cid, tol=tol_spin.value(),
            n_confs=confs_spin.value(), protonate=protonate_check.isChecked(),
            progress_callback=progress_cb,
            use_prefilter=prefilter_check.isChecked(),
            shape_envelope=shape_envelope, active_smiles_list=active_smiles_for_ad,
            match_mode="ensemble" if ensemble_check.isChecked() else "centroid",
        )
        n_before_filters = len(results)
        n_mw_removed = 0
        if mw_filter_check.isChecked() and active_smiles_for_ad:
            results, n_mw_removed = filter_by_mw_range(
                results, active_smiles_for_ad, tolerance_da=mw_tolerance_spin.value()
            )
        if diverse_check.isChecked():
            results = diverse_top_hits(results, max_similarity=0.6)
        state["library_results"] = results

        lib_progress.setVisible(False)
        screen_lib_btn.setEnabled(True)
        msg = f"Screened {n_before_filters} compound(s)."
        if n_prefiltered:
            msg += f" ({n_prefiltered} rejected by the fast pre-filter before full scoring.)"
        if n_skipped:
            msg += f" ({n_skipped} skipped \u2014 failed to parse or embed.)"
        if n_mw_removed:
            lo, hi = mw_range_from_actives(active_smiles_for_ad)
            msg += (f" {n_mw_removed} removed by MW filter (outside "
                    f"{lo-mw_tolerance_spin.value():.0f}\u2013{hi+mw_tolerance_spin.value():.0f} Da).")
        if diverse_check.isChecked():
            msg += f" Showing {len(results)} diverse hit(s) after removing near-duplicate analogs."
        lib_status_label.setText(msg)

        lib_table.setRowCount(len(results))
        for row, r in enumerate(results):
            lib_table.setItem(row, 0, QtWidgets.QTableWidgetItem(r["name"]))
            lib_table.setItem(row, 1, QtWidgets.QTableWidgetItem(f'{r["n_matched"]}/{r["n_total"]}'))
            lib_table.setItem(row, 2, QtWidgets.QTableWidgetItem(f'{r["fraction_matched"]*100:.0f}%'))
            lib_table.setItem(row, 3, QtWidgets.QTableWidgetItem(f'{r["mean_dist"]:.2f}'))
            rmsd_val = r.get("align_rmsd")
            lib_table.setItem(row, 4, QtWidgets.QTableWidgetItem(f'{rmsd_val:.2f}' if rmsd_val is not None else "\u2014"))
            plaus_txt = "\u2014"
            plausibility_val = None
            if active_smiles_for_ad:
                max_sim, _mean_sim = applicability_domain_similarity(r["smiles"], active_smiles_for_ad)
                if max_sim is not None:
                    plaus_txt = f"{max_sim:.2f}"
                    plausibility_val = max_sim
                    r["ad_max_similarity"] = max_sim
            plaus_item = QtWidgets.QTableWidgetItem(plaus_txt)
            if plaus_txt != "\u2014" and float(plaus_txt) < 0.4:
                plaus_item.setForeground(QtGui.QColor("#D97706"))
            lib_table.setItem(row, 5, plaus_item)

            shape_val = r.get("shape_consistency")
            shape_txt = f'{shape_val*100:.0f}%' if shape_val is not None else "\u2014"
            shape_item = QtWidgets.QTableWidgetItem(shape_txt)
            if shape_val is not None and shape_val < 0.7:
                shape_item.setForeground(QtGui.QColor("#D97706"))
            lib_table.setItem(row, 6, shape_item)

            scaffold_txt = r.get("scaffold_novelty") or "\u2014"
            lib_table.setItem(row, 7, QtWidgets.QTableWidgetItem(scaffold_txt))

            verdict_text, verdict_color, verdict_explanation = classify_verdict(
                r["fraction_matched"], plausibility_val
            )
            r["verdict"] = verdict_text
            verdict_item = QtWidgets.QTableWidgetItem(verdict_text)
            verdict_item.setForeground(QtGui.QColor(verdict_color))
            f_verdict = verdict_item.font()
            f_verdict.setBold(True)
            verdict_item.setFont(f_verdict)
            verdict_item.setToolTip(verdict_explanation)
            lib_table.setItem(row, 8, verdict_item)
        export_csv_btn.setEnabled(len(results) > 0)
        export_sdf_btn.setEnabled(len(results) > 0)

    def do_export_csv():
        results = state["library_results"]
        if not results:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            dlg, "Export screening results", "pharmacophore_screening_results.csv", "CSV files (*.csv)"
        )
        if not path:
            return
        import csv
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["name", "smiles", "n_matched", "n_total", "fraction_matched",
                              "mean_dist_angstrom", "align_rmsd_angstrom", "chemical_plausibility",
                              "shape_consistency", "scaffold_novelty", "verdict"])
            for r in results:
                rmsd_val = r.get("align_rmsd")
                shape_val = r.get("shape_consistency")
                writer.writerow([r["name"], r["smiles"], r["n_matched"], r["n_total"],
                                  f'{r["fraction_matched"]:.3f}', f'{r["mean_dist"]:.3f}',
                                  f'{rmsd_val:.3f}' if rmsd_val is not None else "",
                                  f'{r.get("ad_max_similarity", "")}',
                                  f'{shape_val:.3f}' if shape_val is not None else "",
                                  r.get("scaffold_novelty", ""), r.get("verdict", "")])
        lib_status_label.setText(f"Exported {len(results)} result(s) to {path}")

    def do_export_sdf():
        results = state["library_results"]
        if not results:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            dlg, "Export aligned hits", "pharmacophore_screening_aligned_hits.sdf", "SDF files (*.sdf)"
        )
        if not path:
            return
        from rdkit import Chem
        writer = Chem.SDWriter(path)
        n_written = 0
        for r in results:
            mol = r.get("mol")
            conf_id = r.get("conf_id")
            if mol is None or conf_id is None:
                continue
            mol_copy = Chem.Mol(mol)
            mol_copy.SetProp("_Name", r["name"])
            mol_copy.SetProp("fraction_matched", f'{r["fraction_matched"]:.3f}')
            mol_copy.SetProp("mean_dist_angstrom", f'{r["mean_dist"]:.3f}')
            if r.get("align_rmsd") is not None:
                mol_copy.SetProp("align_rmsd_angstrom", f'{r["align_rmsd"]:.3f}')
            if r.get("ad_max_similarity") is not None:
                mol_copy.SetProp("chemical_plausibility", f'{r["ad_max_similarity"]:.3f}')
            writer.write(mol_copy, confId=conf_id)
            n_written += 1
        writer.close()
        lib_status_label.setText(
            f"Exported {n_written} aligned hit structure(s) to {path} \u2014 open directly in PyMOL "
            f"to see exactly how each one aligned to the reference."
        )

    def do_lib_row_clicked(row, _col):
        results = state["library_results"]
        if row < 0 or row >= len(results):
            return
        r = results[row]
        _load_candidate_into_pymol(cmd, r["mol"], r["conf_id"], match_details=r.get("match_details"))
        verdict_label.setStyleSheet("color: #06B6D4; font-weight: 600;")
        verdict_label.setText(
            f'Loaded \u201c{r["name"]}\u201d from library into PyMOL as \u2018candidate\u2019 '
            f'({r["n_matched"]}/{r["n_total"]} features, {r["fraction_matched"]*100:.0f}%).'
        )

    load_example_btn.clicked.connect(do_load_example)
    analyze_btn.clicked.connect(do_analyze)
    thresh_slider.valueChanged.connect(do_threshold_changed)
    eps_slider.valueChanged.connect(do_eps_changed)
    zoom_btn.clicked.connect(do_zoom)
    toggle_ligs_btn.clicked.connect(do_toggle_ligands)
    robustness_btn.clicked.connect(do_run_robustness)
    save_model_btn.clicked.connect(do_save_model)
    load_model_btn.clicked.connect(do_load_model)
    compare_models_btn.clicked.connect(do_compare_models)
    validate_btn.clicked.connect(do_validate)

    def do_ext_actives_changed():
        has_content = bool(ext_actives_edit.toPlainText().strip())
        loocv_check.setEnabled(not has_content)
        if has_content:
            loocv_check.setToolTip(
                "Disabled: leave-one-out doesn't apply to external test actives, since they "
                "were never part of training in the first place \u2014 there's nothing to leave out."
            )
        else:
            loocv_check.setToolTip(
                "Scores each active only against a model rebuilt WITHOUT it, avoiding the "
                "methodological flaw of testing on the same molecules used to build the model. "
                "Automatically not used when external test actives (above) are supplied, since "
                "leave-one-out only makes sense for actives that were part of training."
            )

    ext_actives_edit.textChanged.connect(do_ext_actives_changed)
    export_roc_btn.clicked.connect(do_export_roc)
    export_report_btn.clicked.connect(do_export_report)
    choose_decoy_file_btn.clicked.connect(do_choose_decoy_file)
    score_btn.clicked.connect(do_score_candidate)
    choose_lib_btn.clicked.connect(do_choose_library)
    screen_lib_btn.clicked.connect(do_screen_library)
    export_csv_btn.clicked.connect(do_export_csv)
    export_sdf_btn.clicked.connect(do_export_sdf)
    lib_table.cellClicked.connect(do_lib_row_clicked)

    return dlg


# -------------------------------------------------------------- helpers

def _section_label(text):
    from pymol.Qt import QtWidgets
    lbl = QtWidgets.QLabel(text)
    f = lbl.font()
    f.setBold(True)
    f.setPointSize(f.pointSize() + 1)
    lbl.setFont(f)
    return lbl


def _legend_widget():
    from pymol.Qt import QtWidgets, QtCore
    lbl = QtWidgets.QLabel()
    qtmods = {"QtCore": QtCore, "QtWidgets": QtWidgets}
    rich_text = _compat_enum(qtmods, "QtCore.Qt.TextFormat.RichText", "QtCore.Qt.RichText")
    lbl.setTextFormat(rich_text)
    parts = []
    for fam, (r, g, b) in FAMILY_COLORS.items():
        if fam == "LumpedHydrophobe":
            continue
        hexcol = "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))
        parts.append(f'<span style="color:{hexcol};">\u25CF</span> {fam}')
    lbl.setText("&nbsp;&nbsp;".join(parts))
    return lbl


def _refresh_feature_table(table, pharm, thresh, stability_lookup=None):
    """stability_lookup, if given, maps id(point) -> stability fraction from
    a leave-one-out robustness check; shown as '—' until that's been run."""
    from pymol.Qt import QtWidgets, QtGui
    pts = pharm.points
    table.setRowCount(len(pts))
    for row, p in enumerate(pts):
        is_consensus = p.support_frac >= thresh
        r, g, b = p.color
        fam_item = QtWidgets.QTableWidgetItem(p.family)
        fam_item.setForeground(QtGui.QColor(int(r * 255), int(g * 255), int(b * 255)))
        table.setItem(row, 0, fam_item)
        table.setItem(row, 1, QtWidgets.QTableWidgetItem(f"{p.support_frac*100:.0f}%"))
        table.setItem(row, 2, QtWidgets.QTableWidgetItem(str(p.n_points)))
        stability_txt = "\u2014"
        if stability_lookup is not None and id(p) in stability_lookup:
            stability_txt = f"{stability_lookup[id(p)]*100:.0f}%"
        table.setItem(row, 3, QtWidgets.QTableWidgetItem(stability_txt))
        v_item = QtWidgets.QTableWidgetItem("\u2713" if is_consensus else "\u2014")
        table.setItem(row, 4, v_item)


def _load_into_pymol(cmd, pharm):
    """Load the aligned ligands into the user's live, already-open PyMOL session.
    Works for both a freshly-built LigandBasedPharmacophore (pharm.ligands, a list
    of Ligand objects) and a reloaded LoadedPharmacophore (pharm.ligand_mols, a
    list of (name, Mol) tuples with the conformer already at index 0)."""
    cmd.delete("ligand_* or pharm_* or candidate")
    import tempfile
    from rdkit import Chem
    colors = ["marine", "yelloworange", "hotpink", "aquamarine", "violet",
              "orange", "cyan", "salmon", "limegreen", "lightblue"]

    if getattr(pharm, "ligands", None):
        entries = [(lig.name, lig.mol, lig.conf_id) for lig in pharm.ligands]
    else:
        entries = [(name, mol, 0) for name, mol in getattr(pharm, "ligand_mols", [])]

    for i, (name, mol, conf_id) in enumerate(entries):
        with tempfile.NamedTemporaryFile(suffix=".mol", delete=False, mode="w") as tf:
            path = tf.name
        Chem.MolToMolFile(mol, path, confId=conf_id)
        objname = f"ligand_{i}_{_safe(name)}"
        cmd.load(path, objname)
        cmd.show_as("sticks", objname)
        cmd.color(colors[i % len(colors)], objname)
        os.unlink(path)
    if entries:
        cmd.zoom("ligand_*", buffer=2.0)


def _refresh_pharmacophore_spheres(cmd, pharm, thresh):
    cmd.delete("pharm_*")
    pts = pharm.consensus_at_threshold(thresh)
    for i, p in enumerate(pts):
        name = f"pharm_{i}_{p.family}"
        x, y, z = [float(v) for v in p.centroid]
        cmd.pseudoatom(name, pos=[x, y, z])
        cmd.show("spheres", name)
        cmd.set("sphere_scale", 0.55 + 0.35 * p.support_frac, name)
        cmd.set("sphere_transparency", 0.25, name)
        r, g, b = p.color
        cmd.set_color(f"col_{name}", [r, g, b])
        cmd.color(f"col_{name}", name)


def _load_candidate_into_pymol(cmd, mol, conf_id, match_details=None):
    """Loads a scored candidate into PyMOL as sticks, colored white/CNC by
    default, then highlights the specific atoms that matched each consensus
    feature in that feature's own color -- using match_details from
    score_molecule() (family, color, atom_ids per matched feature). PyMOL's
    atom `index` selector is 1-indexed and matches RDKit's 0-indexed atom
    order exactly (verified directly before this was implemented: loading a
    molecule via MolToMolFile + cmd.load preserves atom order element-for-
    element, PyMOL index N == RDKit atom index N-1, with no exceptions
    across a real test molecule)."""
    import tempfile
    from rdkit import Chem
    cmd.delete("candidate")
    with tempfile.NamedTemporaryFile(suffix=".mol", delete=False, mode="w") as tf:
        path = tf.name
    Chem.MolToMolFile(mol, path, confId=conf_id)
    cmd.load(path, "candidate")
    cmd.show_as("sticks", "candidate")
    cmd.color("white", "candidate")
    cmd.util.cnc("candidate")
    os.unlink(path)

    if match_details:
        for i, detail in enumerate(match_details):
            atom_ids = detail.get("atom_ids", ())
            color = detail.get("color")
            if not atom_ids or not color:
                continue
            color_name = f"cand_match_{i}_{_safe(detail.get('family', 'feat'))}"
            r, g, b = color
            cmd.set_color(color_name, [r, g, b])
            pymol_indices = [str(a + 1) for a in atom_ids]  # RDKit 0-indexed -> PyMOL 1-indexed
            sel = f"candidate and index {'+'.join(pymol_indices)}"
            cmd.color(color_name, sel)
            cmd.show("spheres", sel)
            cmd.set("sphere_scale", 0.25, sel)

    cmd.zoom("candidate or pharm_*", buffer=2.0)


def _safe(name):
    return "".join(c if c.isalnum() else "_" for c in name)[:20]
