#!/usr/bin/env python3

import json
import numpy
from datetime import date
from math import floor, ceil, sqrt
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.pylab as pylab

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
vlines_kwargs = dict(colors='lightgray', linewidth=1, linestyles='dashed', zorder=-1)


# Note, taken from: https://www.statsmodels.org/dev/_modules/statsmodels/graphics/boxplots.html
# And slightly changed
def _single_violin(ax, pos, pos_data, width, side, plot_opts):
    # bw_factor = plot_opts.get('bw_factor', None)

    def _violin_range(pos_data, plot_opts):
        """Return array with correct range, with which violins can be plotted."""
        # cutoff = plot_opts.get('cutoff', False)
        # cutoff_type = plot_opts.get('cutoff_type', 'std')
        # cutoff_val = plot_opts.get('cutoff_val', 1.5)

        # s = 0.0
        # if not cutoff:
        #     if cutoff_type == 'std':
        #         s = cutoff_val * numpy.std(pos_data)
        #     else:
        #         s = cutoff_val

        _range = kde.dataset.max() - kde.dataset.min()
        # if _range < 10000:
        #     cutoff_val = 1.5
        if _range < 20000:
            cutoff_val = 1
        elif _range < 30000:
            cutoff_val = 0.4
        else:
            cutoff_val = 0.15

        s = cutoff_val * numpy.std(pos_data)

        x_lower = kde.dataset.min() - s
        x_upper = kde.dataset.max() + s
        return numpy.linspace(x_lower, x_upper, 100)

    _ran = (max(pos_data) - min(pos_data))
    if _ran < 10000:
        bw_factor = 0.2
    elif _ran < 20000:
        bw_factor = 0.1
    else:
        bw_factor = 0.05

    pos_data = numpy.asarray(pos_data)
    # Kernel density estimate for data at this position.
    kde = gaussian_kde(pos_data, bw_method=bw_factor)

    # Create violin for pos, scaled to the available space.
    xvals = _violin_range(pos_data, plot_opts)
    violin = kde.evaluate(xvals)
    violin = width * violin / violin.max()

    if side == 'both':
        envelope_l, envelope_r = (-violin + pos, violin + pos)
    elif side == 'right':
        envelope_l, envelope_r = (pos, violin + pos)
    elif side == 'left':
        envelope_l, envelope_r = (-violin + pos, pos)
    else:
        msg = "`side` parameter should be one of {'left', 'right', 'both'}."
        raise ValueError(msg)

    # Draw the violin.
    ax.fill_between(xvals, envelope_l, envelope_r,
                     facecolor=plot_opts.get('violin_fc', '#66c2a5'),
                     edgecolor=plot_opts.get('violin_ec', 'k'),
                     lw=plot_opts.get('violin_lw', 1),
                     alpha=plot_opts.get('violin_alpha', 0.5))

    return xvals, violin


def _show(plt, filename_base, params):
    if params.svg:
        plt.savefig(f'{filename_base}.svg')
    plt.savefig(f'{filename_base}.png', dpi=params.dpi)
    if not params.no_preview:
        plt.show()

#
def graph_data(data):
    num_families = len(data)
    orders = []  # List of tuples (clade name, order name, from_idx, list of family names)
    clades = []  # List of tuples (clade name, from_idx)
    for f_data in data:
        clade, order = f_data['clade'], f_data['order']
        if not orders:
            orders.append((clade, order, 0, []))
        elif orders[-1][1] != order:
            orders.append((clade, order, orders[-1][-2] + len(orders[-1][-1]), []))
        orders[-1][-1].append(f_data['family'])
        #
        if not clades or clades[-1][0] != clade:
            clades.append((clade, orders[-1][2]))

    margine = 0.5
    y_pos = numpy.arange(num_families)
    y_min = y_pos[0] - margine
    y_max = y_pos[-1] + margine

    #
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(9, 12), sharey=True)
    fig.subplots_adjust(wspace=0., top=0.96, bottom=0.05, right=0.98, left=0.27)

    # Draw figures
    _num_sequences(axs[0], y_pos, data)
    _seq_dates(axs[2], y_pos, data)
    _seq_lengths(axs[1], y_pos, data)

    # Set axis
    for axis in axs:
        axis.set_ylim([y_min, y_max])
        axis.invert_yaxis()
        axis.spines['top'].set_visible(False)
        # Background
        for _, _, f, fams in orders[1::2]:
            axis.axhspan(f - 0.5, f + len(fams) - 0.5, facecolor='#F1F1F1', zorder=-2)

    for axis in axs[1:]:
        axis.set_yticklabels([])
        axis.tick_params(axis='y', direction='in', length=0, pad=0)

    # Family names
    axis = axs[0]
    axis.axes.set_yticks(y_pos, labels=[f['family'] for f in data])
    axis.tick_params(axis='y', direction='in', length=0, pad=3)

    # Order names
    l_x = 0
    for clade_name, order, from_idx, families in orders:
        l_x = min(l_x, _set_order(fig, axs[0], order, from_idx, from_idx + len(families), families, clade_name))

    # Clade names
    y = clades[1][1] - 0.5
    # Grrr: length last
    # l_x = l_x / 125_000  # Depends on subplots_adjust values. I don't know how :-/
    # Grrr: length in the middle
    l_x = -5.25
    axis.add_artist(lines.Line2D([l_x, l_x + 0.6], [y, y], color='gray', linewidth=2, clip_on=False))
    axis.text(l_x, y - 0.5, clades[0][0], ha='left', va='bottom', color='blue', rotation='vertical', clip_on=False)
    axis.text(l_x, y + 0.5, clades[1][0], ha='left', va='top', color='green', rotation='vertical', clip_on=False)

    # Annotations
    y_h = -margine - 0.5
    txts = ['No. of species\nTotal no. of species', 'Sequence lengths', 'Publication year']
    for axis, c, txt in zip(axs, 'abcd', txts):
        x_min, x_max = axis.get_xlim()
        axis.text(x_min, y_h, c, ha='left', va='baseline', color='black', weight='bold', clip_on=False)
        axis.text((x_min + x_max) / 2, y_h, txt, ha='center', va='baseline', color='black', clip_on=False)

    return plt


def _num_sequences(axis, y_pos, data):
    bar_width = 0.5
    overlap = 0.3
    delta = (bar_width - overlap) / 2
    # n<50, 50<=n<100, 100<=n<200, 200<=n
    # Filipove crvena
    cs = ('#ffffb2', '#fecc5c', '#fd8d3c', '#e31a1c')

    colors = []
    for f_data in data:
        n = f_data['num_sequences']
        if n < 50:
            colors.append(cs[0])
        elif n < 100:
            colors.append(cs[1])
        elif n < 200:
            colors.append(cs[2])
        else:
            colors.append(cs[3])

    op = numpy.log10
    kwargs = dict(height=bar_width)
    axis.barh(y_pos + delta, op([f['num_species'] for f in data]), color='lightgray', **kwargs)
    axis.barh(y_pos - delta, op([f['num_sequences'] for f in data]), color=colors, **kwargs)

    # ToDo: find max value
    values = [1, 2, 3, 4]
    axis.set_xlim([0, 5])
    axis.set_xticks(values)
    # axis.xaxis.set_major_formatter(lambda x, pos: f'$10^{x}$')
    axis.set_xlabel('Logarithm scale')
    axis.vlines(values, *axis.get_ylim(), **vlines_kwargs)


def _seq_dates(axis, y_pos, data):
    fromisoformat = date.fromisoformat

    def _d2n(d):
        d = fromisoformat(d)
        # return d.year + (d.month - 1) / 12
        return d.year

    f_data = [[_d2n(d) for d in f['dates']] for f in data]
    mins = [int(min(ds)) for ds in f_data]
    maxs = [int(max(ds)) for ds in f_data]
    min_d = min(mins)
    max_d = max(maxs)
    parts = axis.violinplot(f_data, positions=y_pos, vert=False, showextrema=False)
    for pc in parts['bodies']:
        pc.set_alpha(0.6)

    #
    axis.set_xticks([min_d, (min_d + max_d) // 2, max_d])
    axis.set_xlabel('Year')


def _seq_lengths(axis, y_pos, data):
    boxes = axis.boxplot([f['lengths'] for f in data],
                         positions=y_pos, vert=False,
                         showcaps=False,
                         medianprops=dict(color='black'),
                         flierprops=dict(marker='o', markersize=3, markeredgecolor='gray'),  # markerfacecolor='black',
                         # flierprops=dict(marker='o', markersize=6,
                         #                 alpha=0.5, markerfacecolor='gray',
                         #                 markeredgecolor='black'),
                         )

    values = [100000, 150000, 200000]
    axis.set_xticks(values)
    # axis.tick_params(axis='x', labelrotation=90)
    axis.xaxis.set_major_formatter(lambda x, pos: f'{x // 1000}')
    axis.set_xlabel('Basepairs (kb)')

    axis.vlines(values, *axis.get_ylim(), **vlines_kwargs)


def _set_order(fig, axis, order, from_idx, to_idx, families, clade_name):
    max_l = 0
    for f in families:
        t = pylab.text(-10, -10, f, color='red')  # Out of canvas
        fig.canvas.draw()
        bbox = pylab.gca().transData.inverted().transform_bbox(t.get_window_extent())
        max_l = max(max_l, abs(bbox.x1 - bbox.x0))

    # Grrr: lengths at the end
    # x = -max_l / 45_000 - 0.2 # Depends on subplots_adjust values. I don't know how :-/
    # Grrr: lengths in the middle
    x = -max_l / 3 - 0.35
    d = 0.2
    cl = 'blue' if clade_name == 'asterids' else 'green'
    axis.add_artist(lines.Line2D([x, x], [from_idx - 0.5 + d, to_idx - 0.5 - d], color='black', linewidth=2, clip_on=False))
    t = axis.text(x - 0.2, (from_idx + to_idx - 1) / 2, order, ha='right', va='center', color=cl, weight='bold', clip_on=False)
    #
    fig.canvas.draw()
    bbox = pylab.gca().transData.inverted().transform_bbox(t.get_window_extent())
    return bbox.x0


#
def graph_iqr_median(iqr_median, histogram_max, histogram_bins):
    fig, axis = plt.subplots(nrows=1, ncols=1, figsize=(8, 4))
    fig.subplots_adjust(bottom=0.12, left=0.07, right=0.96, top=0.96)
    step = 1 / histogram_bins
    v_bins = histogram_max * histogram_bins

    def _p(values):
        ret = [0] * (v_bins + 1)
        for x in values:
            ret[min(v_bins, int(floor(x / step)))] += 1
        return ret

    delta = 0.1
    width = 0.7
    kwargs = dict(width=width, alpha=0.8)
    y_pos = numpy.linspace(0, v_bins, num=v_bins + 1)
    y_pos += delta + width / 2
    axis.bar(y_pos - delta, _p(iqr_median['genera']), color='blue', label='Genera', **kwargs)
    axis.bar(y_pos + delta, _p(iqr_median['families']), color='red', label='Families', **kwargs)

    for side in ('right', 'top'):
        axis.spines[side].set_visible(False)

    axis.set_xlabel('IRQ/median ratio')
    axis.set_xticks([p * histogram_bins for p in range(histogram_max)] + [v_bins + width / 2])
    axis.set_xticklabels([f'{p}%' for p in range(histogram_max)] + [f'>{histogram_max}%'])

    axis.set_ylabel('No. of sequences')
    axis.vlines([histogram_bins], *axis.get_ylim(), colors='lightgray', linewidth=2, linestyles='dashed', zorder=-1)

    # axis.legend()
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(reversed(handles), reversed(labels))  # , title='Line', loc='upper left')

    return plt


#
def _g_axis(axs, idx, nrows, ncols):
    if nrows == 1:
        return axs
    if ncols == 1:
        return axs[idx]
    return axs[idx // ncols, idx % ncols]

def _w_axis(axis, num_x):
    for side in ('left', 'right', 'top', 'bottom'):
        axis.spines[side].set_visible(False)

    # Y axis
    axis.set_ylim([0, 1.2])
    axis.set_yticklabels([])
    axis.tick_params(axis='y', direction='in', length=0, pad=0)

    # X axis
    axis.xaxis.set_major_locator(plt.MaxNLocator(num_x))
    axis.xaxis.set_major_formatter(lambda x, pos: f'{int(x // 1000)}kb')


def _empty_graph(num_graph, nrows, ncols, axs):
    for idx in range(num_graph, nrows * ncols):
        axis = _g_axis(axs, idx, nrows, ncols)
        for side in ('left', 'right', 'top', 'bottom'):
            axis.spines[side].set_visible(False)
        axis.set_yticklabels([])
        axis.tick_params(axis='y', direction='in', length=0, pad=0)
        axis.set_xticklabels([])
        axis.tick_params(axis='x', direction='in', length=0, pad=0)


class _NameColor:
    def __init__(self, clade_names, colors=('blue', 'green', 'red')):
        assert len(colors) >= len(clade_names), f'Not enough colors!!! {len(colors)} {len(clade_names)}'
        self.name_2_idx = dict()
        self.colors = colors
        for idx, names in enumerate(clade_names.values()):
            self.name_2_idx.update((n, idx) for n in names)


    def get_color(self, name):
        return self.colors[self.name_2_idx[name]]


def graph_wide_families(data, clade_names, ncols):
    colors = ['#ff0000', '#cc0033', '#990066', '#660099', '#3300cc', '#0000ff']
    clade_cols = _NameColor(clade_names)
    nrows = int(ceil(len(data) / ncols))
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 7))
    fig.subplots_adjust(hspace=1.2, left=0.02, right=0.98, bottom=0.06, top=0.9)  # wspace=0.2,

    def _one_family(axis, f_data):
        family_col = clade_cols.get_color(f_data['family'])
        num_genera = 0
        num_in_genera = 0
        max_in_gen = len(f_data['genera'][0]['lengths'])
        for g_idx, (g_data, col) in enumerate(zip(f_data['genera'], colors)):
            if len(g_data['lengths']) == 1:
                break
            num_genera += 1
            num_in_genera += len(g_data['lengths'])
            pos = 0.2 + g_idx * 0.15
            # height = 0.3
            height = 0.4 * sqrt(len(g_data['lengths']) / max_in_gen)
            _single_violin(axis, pos, g_data['lengths'], height, 'right', dict(violin_fc=col, violin_ec=col))

        axis.boxplot(f_data['lengths'],
                     vert=False, widths=0.1,
                     positions=[0.1],
                     showfliers=True,
                     boxprops=dict(color=family_col),
                     whiskerprops=dict(color=family_col),
                     flierprops=dict(marker='.'))

        # General
        axis.set_title(f"{f_data['family']} ($n_g$ = {num_genera}, $n_s$ = {num_in_genera})\nIQR/med = {f_data['iqr_median']}%",
                       color=family_col)
        _w_axis(axis, 3)

    for idx, f_data in enumerate(data):
        _one_family(_g_axis(axs, idx, nrows, ncols), f_data)
    _empty_graph(len(data), nrows, ncols, axs)
    return plt


#
def graph_wide_genera(data, clade_names, ncols):
    nrows = int(ceil(len(data) / ncols))
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14, 10))  # , sharey=True)
    fig.subplots_adjust(wspace=0.2, hspace=1.8, left=0.02, right=0.98, bottom=0.05, top=0.95)

    v_color = 'r'
    c_kwargs = dict(violin_fc=v_color, violin_ec=v_color, violin_alpha=0.5)

    def _one_genus(axis, g_data):
        gen_color = _NameColor(clade_names).get_color(g_data['genus'])
        _single_violin(axis, 0.2, g_data['lengths'], 0.8, 'right', c_kwargs)
        axis.boxplot(g_data['lengths'],
                     vert=False, widths=0.1,
                     positions=[0.1],
                     showfliers=True,
                     boxprops=dict(color=gen_color),
                     whiskerprops=dict(color=gen_color),
                     flierprops=dict(marker='.'))

        # General
        axis.set_title(f"${g_data['genus']}$, {g_data['family']}\nn = {len(g_data['lengths'])}, IQR/med = {g_data['iqr_median']}%",
                       color=gen_color)
        _w_axis(axis, 2)


    for idx, g_data in enumerate(data):
        _one_genus(_g_axis(axs, idx, nrows, ncols), g_data)
    _empty_graph(len(data), nrows, ncols, axs)
    return plt


#
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Show and save figures from results')
    parser.add_argument('result_data', help='JSON file created by main script')
    parser.add_argument('figure', choices=('data', 'im', 'iqr_median', 'wf', 'wg'), help='Type of figure to create')

    parser.add_argument('--histogram-max', type=int, default=10, help='IQR/median histogram maximal percentage')
    parser.add_argument('--histogram-bins', type=int, default=10, help='IQR/median histogram number of bins in 1%')
    parser.add_argument('--graph-columns', type=int, help='Number of columns used for wide distribution graphs.')

    parser.add_argument('--dpi', type=int, default=150, help='DPI value use for saving images')
    parser.add_argument('--svg', action='store_true', help='Store image in SVG format')
    parser.add_argument('--no-preview', action='store_true', help='Do not show image')

    params = parser.parse_args()

    with open(params.result_data, 'r', encoding='utf-8') as _in:
        results = json.load(_in)

    base_name = results['base_name']
    if params.figure == 'data':
        _show(graph_data(results['data']), f'{base_name}_data', params)
    elif params.figure in ('im', 'iqr_median'):
        _show(graph_iqr_median(results['iqr_median'], params.histogram_max, params.histogram_bins), f'{base_name}_iqr_median', params)
    elif params.figure == 'wf':
        _show(graph_wide_families(results['wf'], results['clade_names'], params.graph_columns or 2), f'{base_name}_wide_families', params)
    elif params.figure == 'wg':
        _show(graph_wide_genera(results['wg'], results['clade_names'], params.graph_columns or 3), f'{base_name}_wide_genera', params)
    else:
        print(f'Nothing to create for figure {params.type}?!')
