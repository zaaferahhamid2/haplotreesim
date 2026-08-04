"""
Generate manuscript Figures 2-4 (coverage response, phase-switch sensitivity,
and CNA phylogeny accuracy) from results/master_summary.csv.

Usage: python3 plots/generate_figures.py
Run from the project root (haplotreesim/), or update the CSV_PATH below.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CSV_PATH = 'results/master_summary.csv'
OUTPUT_DIR = 'plots'

df = pd.read_csv(CSV_PATH)
for col in ['hscn_error', 'tcn_mse', 'clone_ari', 'loh_f1', 'nrf_distance',
            'tree_coverage', 'mean_allelic_coverage', 'phase_switch']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Fixed style per method: color, marker, linestyle -- keep consistent across all figures
STYLE = {
    'SEACON':      dict(color='#1f77b4', marker='o', ls='-'),
    'Alleloscope': dict(color='#ff7f0e', marker='s', ls='--'),
    'CNRein':      dict(color='#2ca02c', marker='^', ls=':'),
    'SCICoNE':     dict(color='#d62728', marker='D'),
    'CONET':       dict(color='#9467bd', marker='P'),
}
LABELS = {'CNRein': 'CNNaive'}  # display label override; underlying data column stays 'CNRein'
rng = np.random.default_rng(0)


def panel(ax, data, x_col, x_order, y_col, methods, na_methods=None, letter=''):
    """Grouped boxplot with jittered points, one box cluster per x value."""
    n = len(methods)
    width = 0.7
    for xi, xval in enumerate(x_order):
        for mi, method in enumerate(methods):
            sub = data[(data[x_col] == xval) & (data['method'] == method)][y_col].dropna()
            pos = xi + (mi - (n - 1) / 2) * (width / n)
            st = STYLE[method]
            if len(sub) == 0:
                continue
            bp = ax.boxplot([sub.values], positions=[pos], widths=width / n * 0.8,
                             patch_artist=True, showfliers=False, zorder=2)
            bp['boxes'][0].set_facecolor('white')
            bp['boxes'][0].set_edgecolor('black')
            bp['medians'][0].set_color(st['color'])
            bp['medians'][0].set_linewidth(2)
            jitter = rng.uniform(-width / n * 0.15, width / n * 0.15, size=len(sub))
            ax.scatter(np.full(len(sub), pos) + jitter, sub.values,
                       color=st['color'], marker=st['marker'], s=35,
                       alpha=0.85, edgecolors='white', linewidths=0.5, zorder=3)
    if na_methods:
        for i, method in enumerate(na_methods):
            ax.text(0.95, 0.05 + i * 0.09, f"{LABELS.get(method, method)}: N/A",
                    transform=ax.transAxes, ha='right', va='bottom', fontsize=8, style='italic')
    ax.set_xticks(range(len(x_order)))
    ax.set_xticklabels([str(x) for x in x_order])
    ax.set_ylim(bottom=0)
    if letter:
        ax.text(-0.12, 1.05, letter, transform=ax.transAxes, fontsize=16, fontweight='bold')


def two_tool_panel(ax, data, y_col, methods, na_methods=None, letter=''):
    """Simple two-box comparison (used for the SCICoNE/CONET phylogeny figure)."""
    for i, method in enumerate(methods):
        sub = data[data['method'] == method][y_col].dropna()
        st = STYLE[method]
        if len(sub) == 0:
            continue
        bp = ax.boxplot([sub.values], positions=[i], widths=0.5,
                         patch_artist=True, showfliers=False, zorder=2)
        bp['boxes'][0].set_facecolor('white')
        bp['boxes'][0].set_edgecolor('black')
        bp['medians'][0].set_color(st['color'])
        bp['medians'][0].set_linewidth(2)
        jitter = rng.uniform(-0.1, 0.1, size=len(sub))
        ax.scatter(np.full(len(sub), i) + jitter, sub.values,
                   color=st['color'], marker=st['marker'], s=40,
                   alpha=0.85, edgecolors='white', linewidths=0.5, zorder=3)
    if na_methods:
        for i, method in enumerate(na_methods):
            ax.text(0.95, 0.05 + i * 0.09, f"{method}: N/A", transform=ax.transAxes,
                    ha='right', va='bottom', fontsize=8, style='italic')
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods)
    ax.set_ylim(bottom=0)
    ax.text(-0.15, 1.05, letter, transform=ax.transAxes, fontsize=16, fontweight='bold')


def make_legend(fig, methods):
    handles = [plt.Line2D([0], [0], color=STYLE[m]['color'], marker=STYLE[m]['marker'],
                           linestyle=STYLE[m].get('ls', 'None'), markersize=8,
                           label=LABELS.get(m, m)) for m in methods]
    fig.legend(handles=handles, loc='upper center', ncol=len(methods),
               bbox_to_anchor=(0.5, 1.02), frameon=True)


def make_figure2():
    cov = df[df['dataset'].str.startswith('coverage_', na=False)].copy()
    cov_order = [0.005, 0.02, 0.05, 0.1]
    methods = ['SEACON', 'Alleloscope', 'CNRein']

    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
    panel(axes[0], cov, 'mean_allelic_coverage', cov_order, 'hscn_error', methods, letter='A')
    axes[0].set_title('Haplotype-specific state error'); axes[0].set_ylabel('HSCN error')
    panel(axes[1], cov, 'mean_allelic_coverage', cov_order, 'loh_f1', methods, letter='B')
    axes[1].set_title('LOH recovery'); axes[1].set_ylabel('LOH F1')
    panel(axes[2], cov, 'mean_allelic_coverage', cov_order, 'clone_ari',
          ['SEACON', 'Alleloscope'], na_methods=['CNRein'], letter='C')
    axes[2].set_title('Clone assignment'); axes[2].set_ylabel('Clone ARI')
    panel(axes[3], cov, 'mean_allelic_coverage', cov_order, 'tcn_mse',
          ['CNRein'], na_methods=['SEACON', 'Alleloscope'], letter='D')
    axes[3].set_title('Total-copy-number error'); axes[3].set_ylabel('TCN MSE')
    for ax in axes:
        ax.set_xlabel('Mean per-cell coverage (x)')
    make_legend(fig, methods)
    fig.suptitle('Coverage response of HSCN-calling accuracy', y=1.12, fontsize=15)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/figure2_coverage.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {OUTPUT_DIR}/figure2_coverage.png")


def make_figure3():
    ps = df[df['dataset'].str.startswith('phaseswitch_', na=False)].copy()
    ps_order = [0.0, 0.05]
    methods = ['SEACON', 'Alleloscope', 'CNRein']

    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
    panel(axes[0], ps, 'phase_switch', ps_order, 'hscn_error', methods, letter='A')
    axes[0].set_title('HSCN state error'); axes[0].set_ylabel('HSCN error')
    panel(axes[1], ps, 'phase_switch', ps_order, 'loh_f1', methods, letter='B')
    axes[1].set_title('LOH recovery'); axes[1].set_ylabel('LOH F1')
    panel(axes[2], ps, 'phase_switch', ps_order, 'clone_ari',
          ['SEACON', 'Alleloscope'], na_methods=['CNRein'], letter='C')
    axes[2].set_title('Clone assignment'); axes[2].set_ylabel('Clone ARI')
    panel(axes[3], ps, 'phase_switch', ps_order, 'tcn_mse',
          ['CNRein'], na_methods=['SEACON', 'Alleloscope'], letter='D')
    axes[3].set_title('Total-copy-number error'); axes[3].set_ylabel('TCN MSE')
    for ax in axes:
        ax.set_xlabel('Phase-switch probability')
    make_legend(fig, methods)
    fig.suptitle('Sensitivity to phase-switch error at fixed coverage', y=1.12, fontsize=15)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/figure3_phaseswitch.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {OUTPUT_DIR}/figure3_phaseswitch.png")


def make_figure4():
    core = df[~df['dataset'].str.contains('coverage_|phaseswitch_|wgd_|doublet_', na=False)].copy()
    methods = ['SCICoNE', 'CONET']

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    two_tool_panel(axes[0], core, 'tcn_mse', methods, letter='A')
    axes[0].set_title('Total-CN accuracy'); axes[0].set_ylabel('TCN MSE')
    two_tool_panel(axes[1], core, 'clone_ari', methods, letter='B')
    axes[1].set_title('Clone assignment'); axes[1].set_ylabel('Clone ARI')
    two_tool_panel(axes[2], core, 'nrf_distance', methods, na_methods=['SCICoNE'], letter='C')
    axes[2].set_title('Tree topology'); axes[2].set_ylabel('Normalized RF distance')
    two_tool_panel(axes[3], core, 'tree_coverage', methods, na_methods=['SCICoNE'], letter='D')
    axes[3].set_title('Tree coverage'); axes[3].set_ylabel('Tree coverage')
    fig.suptitle('CNA phylogeny accuracy, clone recovery, and tree completeness', y=1.08, fontsize=15)
    fig.text(0.5, -0.03,
             'Boxes show the median and interquartile range; points show individual replicate instances.',
             ha='center', fontsize=9, style='italic')
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/figure4_phylogeny.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {OUTPUT_DIR}/figure4_phylogeny.png")


if __name__ == '__main__':
    make_figure2()
    make_figure3()
    make_figure4()
    print("\nAll figures regenerated.")
