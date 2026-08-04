# Manuscript Figures

This folder contains the manuscript's data-driven figures (Figures 2-4) and the script that generates them.

## Contents

- `figure2_coverage.png` — Coverage response of HSCN-calling accuracy (SEACON, Alleloscope, CNNaive)
- `figure3_phaseswitch.png` — Sensitivity to phase-switch error at fixed coverage
- `figure4_phylogeny.png` — CNA phylogeny accuracy, clone recovery, and tree completeness (SCICoNE, CONET)
- `generate_figures.py` — Standalone script that regenerates all three figures from `results/master_summary.csv`

## Regenerating the figures

From the project root:

```bash
pip3 install pandas matplotlib
python3 plots/generate_figures.py
```

This reads `results/master_summary.csv` (all parameters and metrics across all completed benchmark runs) and writes updated PNGs back into this folder.

## Editing the style

Each method has a fixed color/marker/linestyle defined in the `STYLE` dictionary near the top of `generate_figures.py`:

- SEACON: blue circle, solid line
- Alleloscope: orange square, dashed line
- CNRein (displayed as "CNNaive"): green triangle, dotted line
- SCICoNE: red diamond
- CONET: purple plus

To change which metrics appear in each panel, or to add/remove panels, edit the `make_figure2()`, `make_figure3()`, and `make_figure4()` functions — each one builds a 1x4 grid of subplots using the shared `panel()` (grouped by x-axis condition) or `two_tool_panel()` (simple two-tool comparison) helper functions.

Note: Figure 4's Panels C and D (tree topology metrics) show "SCICoNE: N/A" because `scripts/evaluate_scicone.py` does not currently compute normalized RF distance or tree coverage — see the Known Limitations section of `PROJECT_FILES.md`.
