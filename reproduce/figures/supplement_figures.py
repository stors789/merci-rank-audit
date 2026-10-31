"""Redraw supplementary figures S1-S4 from fixed aggregate source tables."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter

from figure_style import TEAL, TERRA, SLATE, MUTED, GRID, panel, grid, save


def read(data, name):
    path = data / name
    if not path.exists() and (data.parent / 'bio06' / name).exists():
        path = data.parent / 'bio06' / name
    return pd.read_csv(path, float_precision='round_trip')


def gate(data, out):
    d = read(data, 'FigureS1_gate.csv')
    cases = [('bench1_split21_negative', 'B1 split 21 negative'), ('bench2_split1_negative', 'B2 split 1 negative'),
             ('P2', 'P2'), ('GBML181', 'GBML181')]
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.05), sharex=True, sharey=True)
    for ax, (case, label), letter in zip(axes.flat, cases, 'ABCD'):
        z = d[d.case == case]
        for method, color, ls, marker in [('uniform', TERRA, '-', 'o'),
                                           ('conditional_weighted', TEAL, '--', 's')]:
            q = z[z.method == method].sort_values('target_seconds')
            x = q.target_seconds.to_numpy(float)
            y, lo, hi = (q[c].to_numpy(float) for c in ('gate_median', 'gate_q25', 'gate_q75'))
            ax.fill_between(x, np.where(lo > 0, lo, np.nan), np.where(hi > 0, hi, np.nan),
                            color=color, alpha=.10, linewidth=0)
            ax.plot(x, np.where(y > 0, y, np.nan), color=color, ls=ls,
                    marker=marker, ms=3.3, lw=1.15)
        half = float(z.reference_gate_95_halfwidth.iloc[0])
        if half > 0:
            ax.axhline(half, color=SLATE, lw=.9, ls=':')
        p = float(z.reference_P_G.iloc[0])
        panel(ax, letter, f'{label} · P(G) ≈ {p:.4g}')
        ax.set(xscale='log', yscale='log', xlim=(.075, 13), ylim=(1e-8, .03))
        ax.set_xticks([.1, .5, 2, 10], ['0.1', '0.5', '2', '10'])
        ax.set_yticks([1e-7, 1e-5, 1e-3]); grid(ax)
    axes[0, 0].set_ylabel('Absolute gate error')
    axes[1, 0].set_ylabel('Absolute gate error')
    axes[1, 0].set_xlabel('Sampling time (s)'); axes[1, 1].set_xlabel('Sampling time (s)')
    fig.legend(handles=[Line2D([], [], color=TERRA, marker='o', ms=3.5, lw=1.15, label='Uniform'),
                        Line2D([], [], color=TEAL, marker='s', ls='--', ms=3.5, lw=1.15, label='Weighted'),
                        Line2D([], [], color=SLATE, ls=':', lw=1, label='Reference 95% MC half-width')],
               loc='upper center', ncol=3, frameon=False, bbox_to_anchor=(.51, .99),
               handlelength=1.6, columnspacing=1.0)
    fig.text(.09, .055, 'Certain controls', fontsize=8, weight='bold')
    fig.text(.27, .055, 'B1 split 21 mixed  ·  B2 split 1 mixed  ·  MC only', fontsize=7.8)
    fig.text(.27, .015, 'P(G) = 1; zero error at all four times, 20/20 in each.', fontsize=7.5, color=MUTED)
    fig.subplots_adjust(left=.12, right=.98, bottom=.24, top=.82, hspace=.64, wspace=.23)
    save(fig, out, 'FigureS1_gate')


def precision(data, out):
    summary = read(data, 'FigureS3_precision.csv')
    reference = read(data, 'FigureS3_reference_quantiles.csv')
    series = [('original', 'HALLMARK_MTORC1_SIGNALING'), ('original', 'HALLMARK_PROTEIN_SECRETION'),
              ('alternative', 'HALLMARK_MTORC1_SIGNALING'), ('alternative', 'HALLMARK_PROTEIN_SECRETION')]
    colors = {.025: SLATE, .5: TEAL, .975: TERRA}
    markers = {.025: 'o', .5: 's', .975: '^'}
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.30), sharey=True, sharex=True)
    for ax, (rep, feature), letter in zip(axes.flat, series, 'ABCD'):
        part = summary[(summary.representation == rep) & (summary.feature == feature)]
        ref = reference[(reference.representation == rep) & (reference.feature == feature)]
        for quantile in [.025, .5, .975]:
            z = part[part['quantile'] == quantile].set_index('budget').reindex([100, 1000, 5000, 20000])
            y = 100*z.difference_over_old_width_median.to_numpy(float)
            lo = 100*z.difference_over_old_width_q25.to_numpy(float)
            hi = 100*z.difference_over_old_width_q75.to_numpy(float)
            color = colors[quantile]
            ax.fill_between(range(4), lo, hi, color=color, alpha=.08, linewidth=0)
            ax.plot(range(4), y, color=color, marker=markers[quantile], ms=3.7, lw=1.1)
            rr = ref[ref['quantile'] == quantile].iloc[0]
            ax.axhline(100*rr.reference_max_half_width/rr.old_interval_width, color=color, ls=':', lw=.8)
        ax.axhline(2.5, color=MUTED, ls='--', lw=.75)
        title = ('mTORC1' if 'MTORC1' in feature else 'Secretion') + ' · ' + rep
        panel(ax, letter, title); ax.set_ylim(0, 3.4); ax.set_yticks([0, 1, 2, 3])
        ax.set_xlim(-.15, 3.15); ax.set_xticks(range(4), ['100', '1,000', '5,000', '20,000'])
        grid(ax)
    axes[0, 0].set_ylabel('Difference / old width (%)')
    axes[1, 0].set_ylabel('Difference / old width (%)')
    axes[1, 0].set_xlabel('Proposal budget'); axes[1, 1].set_xlabel('Proposal budget')
    fig.legend(handles=[Line2D([], [], color=colors[q], marker=markers[q], ms=3.5, lw=1.1, label=label)
                        for q, label in [(.025, '2.5%'), (.5, 'Median'), (.975, '97.5%')]] +
                       [Line2D([], [], color=MUTED, ls='--', lw=.8, label='2.5% scale')],
               loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(.52, .985),
               handlelength=1.4, columnspacing=1.0)
    fig.text(.5, .025, '100 proposals: 0/20 defined  ·  1,000–20,000: 20/20 defined  ·  dotted: reference MC error',
             ha='center', fontsize=7.5, color=MUTED)
    fig.subplots_adjust(left=.12, right=.98, bottom=.18, top=.82, hspace=.62, wspace=.18)
    save(fig, out, 'FigureS2_precision')


def union(data, out):
    d = read(data, 'FigureS4_union.csv')
    budgets = [100, 1000, 5000, 20000]
    fig = plt.figure(figsize=(7.0, 4.15))
    ax = fig.add_axes([.11, .57, .83, .29])
    table = fig.add_axes([.11, .12, .83, .29]); table.axis('off')
    for method, color, marker, ls, face in [('uniform', TERRA, 'o', '-', TERRA),
                                             ('conditional', TEAL, 's', '--', 'white')]:
        q = d[(d.case == 'P2') & (d.method == method)].set_index('budget').reindex(budgets)
        y = 100*q.union_unseen.to_numpy(float)/q.exact_passing_union.to_numpy(float)
        ax.plot(range(4), y, color=color, ls=ls, marker=marker, markerfacecolor=face,
                markeredgecolor=color, ms=4.2, lw=1.2)
    ax.set(ylim=(-4, 104), xlim=(-.1, 3.1)); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(range(4), ['100', '1k', '5k', '20k'])
    ax.set_ylabel('Exact union unseen (%)'); ax.set_xlabel('Proposal budget'); grid(ax)
    panel(ax, 'A', 'P2 · rare gate')
    ax.legend(handles=[Line2D([], [], color=TERRA, marker='o', ms=4, lw=1.1, label='Uniform'),
                       Line2D([], [], color=TEAL, marker='s', markerfacecolor='white', ms=4, ls='--', lw=1.1, label='Conditional')],
              loc='lower right', bbox_to_anchor=(1, 1.19), ncol=2, frameon=False,
              handlelength=1.4, handletextpad=.4, columnspacing=.8)
    table.text(0, 1.11, 'B   Accepted orders', transform=table.transAxes, fontsize=8.6, weight='bold')
    table.text(.37, 1.11, 'uniform / conditional at each budget', transform=table.transAxes, fontsize=7.5, color=MUTED)
    xs = [.40, .59, .78, .98]
    for x, label in zip(xs, ['100', '1k', '5k', '20k']):
        table.text(x, .79, label, transform=table.transAxes, ha='right', fontsize=7.5, weight='bold')
    rows = [('all_tied_certain', 'N8 certain'), ('all_tied_single', 'N8 single'),
            ('all_tied_overlapping', 'N8 overlap'), ('P2', 'P2')]
    for i, (case, label) in enumerate(rows):
        y = .65-i*.16
        table.text(0, y, label, transform=table.transAxes, fontsize=7.5, va='center')
        a = d[(d.case == case) & (d.method == 'uniform')].set_index('budget').reindex(budgets)
        b = d[(d.case == case) & (d.method == 'conditional')].set_index('budget').reindex(budgets)
        for x, u, c in zip(xs, a.accepted, b.accepted):
            table.text(x, y, f'{int(u):,}/{int(c):,}', transform=table.transAxes,
                       ha='right', fontsize=7.5, va='center')
        table.plot([0, 1], [y-.08, y-.08], color=GRID, lw=.5, transform=table.transAxes, clip_on=False)
    fig.text(.11, .045, 'N8: 0% unseen for both methods. P2 uniform at 100: zero accepted; frequencies undefined.',
             fontsize=7.5, color=MUTED)
    save(fig, out, 'FigureS3_union')


def sensitivity(data, out):
    d = read(data, 'rarity_comparison.csv')
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.4), sharey=True)
    for ax, case, letter in zip(axes, ['P2', 'GBML181'], 'AB'):
        z = d[d.case == case].set_index('null_id').loc[['original'] + [f'null_{i:02d}' for i in range(1, 11)]]
        x = np.arange(11); y = z.P_G.to_numpy(float); err = z.MC_95_halfwidth.to_numpy(float)
        ax.axhspan(y[0]-err[0], y[0]+err[0], color=TEAL, alpha=.12, linewidth=0)
        ax.axhline(y[0], color=TEAL, lw=.9, ls='--')
        ax.errorbar(x[1:], y[1:], yerr=err[1:], fmt='o', color=TERRA, ms=3.7, capsize=1.5, lw=.9)
        ax.errorbar([0], [y[0]], yerr=[err[0]], fmt='s', color=TEAL, ms=4, capsize=1.5, lw=.9)
        ax.set_yscale('log'); ax.set_ylim(.0002, .002); ax.set_xlim(-.65, 10.65)
        ax.set_xticks(x, ['Old'] + [str(i) for i in range(1, 11)])
        ax.set_yticks([.0002, .0005, .001, .002], labels=['0.0002', '0.0005', '0.001', '0.002'])
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel('Fixed threshold seed'); grid(ax); panel(ax, letter, case)
    axes[0].set_ylabel('P(gate passes)')
    fig.legend(handles=[Line2D([], [], color=TEAL, marker='s', ls='--', ms=4, lw=1, label='Original threshold'),
                        Line2D([], [], color=TERRA, marker='o', ls='', ms=4, label='New null seeds')],
               loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(.53, .99),
               handlelength=1.5, columnspacing=1)
    fig.subplots_adjust(left=.12, right=.985, bottom=.24, top=.73, wspace=.18)
    save(fig, out, 'FigureS4_bio06')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--figures', nargs='+', choices=['S1', 'S2', 'S3', 'S4'],
                   default=['S1', 'S2', 'S3', 'S4'])
    a = p.parse_args()
    for n in a.figures:
        {'S1': gate, 'S2': precision, 'S3': union, 'S4': sensitivity}[n](a.data, a.output)
    print('Redrew supplement figures:', ', '.join(a.figures))


if __name__ == '__main__':
    main()
