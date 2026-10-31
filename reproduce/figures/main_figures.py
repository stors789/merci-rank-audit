"""Redraw main figures 1-5 from the fixed aggregate source tables."""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from figure_style import INK, TEAL, TERRA, SLATE, GRID, MUTED, panel, grid, save


def read(data, name):
    return pd.read_csv(data / name, float_precision='round_trip')


def workflow(out):
    fig, ax = plt.subplots(figsize=(7.0, 2.65))
    ax.set(xlim=(0, 100), ylim=(0, 36)); ax.axis('off')
    positions = [2, 27, 52, 77]
    for x, letter, title in zip(positions, 'ABCD',
            ['Fixed inputs', 'Orders + gate', 'Label families', 'Audit outputs']):
        ax.text(x, 33.4, letter, weight='bold', fontsize=9.6, va='center')
        ax.text(x + 4, 33.4, title, weight='bold', fontsize=8.5, va='center')
    for x in (25, 50, 75):
        ax.plot([x, x], [5.5, 30], color=GRID, lw=.7)

    def chips(x, y, letters, selected=(), width=2.28):
        for i, letter in enumerate(letters):
            xx=x+i*width+(0.65 if i>=4 else 0)
            ax.add_patch(Rectangle((xx, y), width-.14, 3.05, linewidth=.55,
                                   edgecolor=GRID, facecolor='#E8F3F2' if letter in selected else 'white'))
            ax.text(xx+(width-.14)/2, y+1.52, letter, ha='center', va='center',
                    color=TEAL if letter in selected else INK, fontsize=7.5)
        ax.plot([x+4*width+.25, x+4*width+.25], [y-.45, y+3.5], color=TEAL, lw=.9)

    rna=set('aceg')
    ax.text(2, 28, 'RNA rank · top 4', fontsize=8, va='center')
    chips(2, 23.2, 'acegbdfh', rna)
    ax.text(2, 19.3, 'Neff tie blocks', fontsize=8, va='center')
    for x,w,label in [(2,4.1,'a b'),(6.5,9.1,'c d e f'),(16,4.1,'g h')]:
        ax.add_patch(Rectangle((x,14.5),w,3.4,linewidth=.7,edgecolor=SLATE,facecolor='white'))
        ax.text(x+w/2,16.2,label,ha='center',va='center',fontsize=7.5,color=INK)
    ax.text(2, 10.5, 'Gate maximum = 2', fontsize=7.8, va='center')

    ax.text(27, 28, 'One gate: overlap > 2', fontsize=7.8, va='center')
    ax.text(27, 24.8, 'I', fontsize=8, weight='bold', va='center')
    chips(29, 23.1, 'abcdefgh', set('ac'))
    ax.text(29, 20.1, '2 selected  ·  fail', fontsize=7.5, color=TERRA, va='center')
    ax.text(27, 16.8, 'II', fontsize=8, weight='bold', va='center')
    chips(29, 15.1, 'baecdfhg', set('ace'))
    ax.text(29, 12.1, '3 selected  ·  pass', fontsize=7.5, color=TEAL, va='center')
    ax.text(29, 8.4, 'DNA top 4  |  remainder', fontsize=7.5, color=MUTED, va='center')

    ax.text(65.5, 27.8, 'I', fontsize=7.5, weight='bold', ha='center', va='center')
    ax.text(72.0, 27.8, 'II', fontsize=7.5, weight='bold', ha='center', va='center')
    for y,color,label,left,right in [(23,SLATE,'Raw','{a,c}','{a,c,e}'),
                                    (17,TERRA,'Default','∅','{a,c,e}'),
                                    (11,TEAL,'Passing-\nonly','—','{a,c,e}')]:
        ax.plot([52,53.3],[y,y],color=color,lw=2,solid_capstyle='round')
        ax.text(54, y, label, fontsize=7.8, weight='bold', va='center')
        ax.text(65.5,y,left,fontsize=7.5,ha='center',va='center')
        ax.text(72.0,y,right,fontsize=7.5,ha='center',va='center')
    ax.text(52, 7.0, '— excluded from passing', fontsize=7.5, color=MUTED, va='center')

    ax.text(77, 27, 'Exact · all orders', fontsize=8, weight='bold', color=TEAL, va='center')
    ax.text(77, 22.8, '∩ core   ∪ union   ✓ witness', fontsize=7.5, va='center')
    ax.text(77, 17.8, 'Analytic marginals', fontsize=7.8, weight='bold', va='center')
    ax.text(77, 14.3, 'certain / single gate', fontsize=7.5, color=MUTED, va='center')
    ax.text(77, 9.9, 'Sampled orders  →  nuclear', fontsize=7.5, va='center')
    for x in (25, 50, 75):
        ax.annotate('', xy=(x+1.3, 5.2), xytext=(x-1.3, 5.2),
                    arrowprops={'arrowstyle': '-|>', 'color': SLATE, 'lw': .7, 'mutation_scale': 7})
    ax.text(2, 2.2, 'Single-gate schematic · K = 4 · two orders shown, not exhaustive',
            fontsize=7.5, color=MUTED, va='center')
    fig.subplots_adjust(left=.01, right=.99, bottom=.04, top=.98)
    save(fig, out, 'Figure1_workflow')


def sampling(data, out, aspect=1.4):
    d = read(data, 'figure2_equal_time.csv')
    cases = [('P2', 'P2 · rare'), ('GBML181', 'GBML181 · rare'),
             ('bench2_split1_negative', 'B2 negative · common')]
    styles = [('uniform', TERRA, '-', 'o', TERRA),
              ('conditional_weighted', TEAL, '-', 's', TEAL),
              ('conditional_accepted', TEAL, 'None', 'D', 'white')]
    left, right, gap = .10, .985, .035
    width = (right-left-2*gap)/3
    plot_height = aspect*width*7
    height = plot_height+1.20
    fig = plt.figure(figsize=(7, height))
    axes = [fig.add_axes([left+i*(width+gap), .52/height, width, plot_height/height]) for i in range(3)]
    for ax, (case, title), letter in zip(axes, cases, 'ABC'):
        sub = d[d.case == case]
        for method, color, ls, marker, face in (styles[0], styles[2], styles[1]):
            z = sub[sub.method == method].sort_values('target_seconds')
            x = z.target_seconds.to_numpy(float)
            y = z.member_median.to_numpy(float)
            lo, hi = z.member_q25.to_numpy(float), z.member_q75.to_numpy(float)
            ax.fill_between(x, lo, hi, color=color, alpha=.08, linewidth=0)
            size = 4.8 if method == 'conditional_accepted' else 3.1
            ax.plot(x, y, color=color, ls=ls, marker=marker, lw=1.05, ms=size,
                    markerfacecolor=face, markeredgecolor=color, markeredgewidth=.8)
            if method == 'uniform':
                for _, r in z[z.member_evaluable < 20].iterrows():
                    ax.annotate(f'{int(r.member_evaluable)}/20', (r.target_seconds, r.member_median),
                                xytext=(2, 4), textcoords='offset points', color=TERRA, fontsize=7.5)
        ax.set(xscale='log', yscale='log', ylim=(.0004, 1), xlim=(.075, 13))
        ax.set_xticks([.1, .5, 2, 10], labels=['0.1', '0.5', '2', '10'])
        ax.set_yticks([.001, .01, .1, 1])
        ax.grid(False); ax.tick_params(which='minor', length=0); panel(ax, letter, title)
    axes[0].set_ylabel('Member RMSE')
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    fig.text(.52, .018, 'Sampling time (s)', ha='center', fontsize=8.2)
    handles = [Line2D([], [], color=c, ls=ls, marker=m, markerfacecolor=face,
                      markeredgecolor=c, markeredgewidth=.8, lw=1.05, ms=4, label=label)
               for (_, c, ls, m, face), label in zip(styles, ['Uniform', 'Weighted', 'Accepted'])]
    fig.legend(handles=handles, loc='upper center', ncol=3, frameon=False,
               bbox_to_anchor=(.52, .995), handlelength=1.5, columnspacing=1.25)
    save(fig, out, 'Figure2_sampling')


def annotation(data, out):
    errors = read(data, 'figure3_error_partitions.csv')
    changes = read(data, 'figure3_paired_changes.csv')
    stages = ['frozen_evaluation', 'retrospective_second_experiment']
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.25), sharey=True)
    for i, stage in enumerate(stages):
        for ax, table, specs in [
            (axes[0], errors, [('FP_core_fraction', TEAL, 'o', -.13),
                               ('FN_outside_fraction', TERRA, 's', .13)]),
            (axes[1], changes, [('precision_delta', TEAL, 'o', -.13),
                                ('sensitivity_delta', TERRA, 's', .13)])]:
            for metric, color, marker, offset in specs:
                r = table[(table.stage == stage) & (table.metric == metric)].iloc[0]
                med, lo, hi = (100*float(r[x]) for x in ('median', 'q25', 'q75'))
                ax.errorbar(med, i+offset, xerr=[[med-lo], [hi-med]], fmt=marker,
                            color=color, lw=1.2, ms=4.2, capsize=2)
    for ax in axes:
        ax.set_ylim(1.55, -.55)
        ax.set_yticks([0, 1], ['B1 frozen (n = 80)', 'B2 retro (n = 20)'])
        grid(ax, 'x')
    axes[0].set_xlim(70, 101); axes[0].set_xticks([70, 80, 90, 100])
    axes[0].set_xlabel('Fraction of original errors (%)')
    axes[1].set_xlim(-4.4, 3); axes[1].set_xticks([-4, -2, 0, 2])
    axes[1].axvline(0, color=SLATE, lw=.8, ls=':')
    axes[1].set_xlabel('Paired change (percentage points)')
    panel(axes[0], 'A', 'Stable error partition')
    panel(axes[1], 'B', 'Core versus original call')
    left = [Line2D([], [], color=TEAL, marker='o', ls='', ms=4.2, label='FP in core'),
            Line2D([], [], color=TERRA, marker='s', ls='', ms=4.2, label='FN outside')]
    right = [Line2D([], [], color=TEAL, marker='o', ls='', ms=4.2, label='Precision'),
             Line2D([], [], color=TERRA, marker='s', ls='', ms=4.2, label='Sensitivity')]
    axes[0].legend(handles=left, loc='lower left', bbox_to_anchor=(-.02, 1.22), ncol=2,
                   frameon=False, handletextpad=.4, columnspacing=.9)
    axes[1].legend(handles=right, loc='lower left', bbox_to_anchor=(-.02, 1.22), ncol=2,
                   frameon=False, handletextpad=.4, columnspacing=.9)
    fig.subplots_adjust(left=.19, right=.985, bottom=.25, top=.68, wspace=.19)
    save(fig, out, 'Figure3_annotation')


def composition(data, out, aspect=1.4):
    curves = read(data, 'figure4_core_curves.csv')
    transitions = read(data, 'figure4_endpoint_transitions.csv')
    stages = [('frozen_evaluation', TEAL, 'o', 'B1 frozen'),
              ('retrospective_second_experiment', TERRA, 's', 'B2 retrospective')]
    left_width = .29
    plot_height = aspect*left_width*7
    height = plot_height+1.20
    fig = plt.figure(figsize=(7, height))
    bottom = .52
    a = fig.add_axes([.12, bottom/height, left_width, plot_height/height])
    forest_height = 1.22
    forest_bottom = bottom+(plot_height-forest_height)/2
    b = fig.add_axes([.60, forest_bottom/height, .36, forest_height/height])
    for stage, color, marker, label in stages:
        z = curves[curves.stage == stage].sort_values('level')
        x = z.level.to_numpy(float)
        y, lo, hi = (100*z[c].to_numpy(float) for c in ('median', 'q25', 'q75'))
        a.fill_between(x, lo, hi, color=color, alpha=.12, linewidth=0)
        a.plot(x, y, color=color, marker=marker, ms=3.8, lw=1.2, label=label)
    a.set_xlim(-5, 105); a.set_ylim(0, 30)
    a.set_yticks([0, 10, 20, 30]); a.set_xticks([0, 25, 50, 75, 100])
    a.set_xlabel('Background slots replaced by CC (%)'); a.set_ylabel('Anchor in core (%)')
    a.grid(False)
    a.text(0, 1.045, 'A', transform=a.transAxes, fontsize=9.6, fontweight='bold', color=INK, va='bottom')
    a.text(.105, 1.045, 'Fixed-N trend', transform=a.transAxes, fontsize=8.5,
           fontweight='bold', color=INK, va='bottom')
    handles, labels = a.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(.12, .995), ncol=2,
               frameon=False, handlelength=1.4, handletextpad=.4, columnspacing=.8)
    for i, (stage, color, _, _) in enumerate(stages):
        for design, offset, marker, face in [('fixedN', -.13, 'o', color), ('nested', .13, 's', 'white')]:
            r = transitions[(transitions.stage == stage) & (transitions.design == design)].iloc[0]
            med, lo, hi = (100*float(r[x]) for x in ('median', 'q25', 'q75'))
            b.errorbar(med, i+offset, xerr=[[med-lo], [hi-med]], fmt=marker, color=color,
                       markerfacecolor=face, markeredgecolor=color, markeredgewidth=.9,
                       lw=1.2, ms=4.2, capsize=2)
    b.set_ylim(1.55, -.55); b.set_yticks([0, 1], ['B1 frozen', 'B2 retro'])
    b.set_xlim(0, 20.5); b.set_xticks([0, 5, 10, 15, 20])
    b.set_xlabel('Core → outside (% of reference)'); b.grid(False)
    panel(b, 'B', 'Endpoint change')
    b.legend(handles=[Line2D([], [], color=INK, marker='o', ls='', ms=4, label='Fixed N'),
                      Line2D([], [], color=INK, marker='s', markerfacecolor='white', ls='', ms=4, label='Nested')],
             loc='lower left', bbox_to_anchor=(-.02, 1.20), ncol=2, frameon=False,
             handletextpad=.4, columnspacing=.8)
    save(fig, out, 'Figure4_composition')


def nuclear(data, out):
    ref = read(data, 'figure5_reference_quantiles.csv')
    fixed = pd.concat([read(data, 'figure5_fixed_original.csv'), read(data, 'figure5_fixed_alternative.csv')])
    witness = read(data, 'figure5_selected_witness.csv')
    features = ['HALLMARK_MTORC1_SIGNALING', 'HALLMARK_PROTEIN_SECRETION']
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.36), sharey=True)
    for ax, rep, letter in zip(axes, ['original', 'alternative'], 'AB'):
        for y, feature in enumerate(features):
            z = ref[(ref.representation == rep) & (ref.feature == feature)].sort_values('quantile')
            q = z.reference_estimate.to_numpy(float)
            raw = fixed[(fixed.representation == rep) & (fixed.nuclear_feature == feature)].beta.iloc[0]
            selected = witness[(witness.representation == rep) & (witness.nuclear_feature == feature)].beta.iloc[0]
            ax.plot([q[0], q[2]], [y, y], color=TEAL, lw=2.1, solid_capstyle='round')
            ax.plot(q[1], y, 'o', color=TEAL, ms=3.8, markeredgecolor='white', markeredgewidth=.5)
            ax.plot(raw, y, 'D', color=INK, ms=3.7)
            ax.plot(selected, y, 'x', color=TERRA, ms=4.4, mew=1.05)
        ax.axvline(0, color=SLATE, lw=.8, ls=':')
        ax.set_yticks([0, 1], ['mTORC1', 'Protein secretion']); ax.set_ylim(1.5, -.5)
        ax.set_xlabel('β coefficient'); grid(ax, 'x'); panel(ax, letter, rep.capitalize())
    axes[0].set_xlim(-.012, .035); axes[0].set_xticks([-.01, 0, .01, .02, .03])
    axes[1].set_xlim(-.0022, .0092); axes[1].set_xticks([-.002, 0, .002, .004, .006, .008])
    handles = [Line2D([], [], color=TEAL, lw=2.1, marker='o', ms=3.8, markeredgecolor='white',
                      label='Passing: median and 2.5–97.5%'),
               Line2D([], [], color=INK, marker='D', ls='', ms=4, label='Fixed raw'),
               Line2D([], [], color=TERRA, marker='x', ls='', ms=5, label='Selected witness')]
    fig.legend(handles=handles, loc='upper center', ncol=3, frameon=False,
               bbox_to_anchor=(.54, .995), handlelength=1.6, columnspacing=1.15)
    fig.subplots_adjust(left=.16, right=.985, bottom=.25, top=.73, wspace=.28)
    save(fig, out, 'Figure5_nuclear')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--figures', nargs='+', choices=['1', '2', '3', '4', '5'], default=['1', '2', '3', '4', '5'])
    p.add_argument('--aspect', type=float, choices=[1.4, 1.6], default=1.4,
                   help='Height/width of each curve panel')
    a = p.parse_args()
    for number in a.figures:
        {'1': lambda: workflow(a.output), '2': lambda: sampling(a.data, a.output, a.aspect),
         '3': lambda: annotation(a.data, a.output), '4': lambda: composition(a.data, a.output, a.aspect),
         '5': lambda: nuclear(a.data, a.output)}[number]()
    print('Redrew main figures:', ', '.join(a.figures))


if __name__ == '__main__':
    main()
