import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple
import numpy as np
import subprocess
import json

def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

# solver key, legend label, color (Okabe-Ito colorblind-safe palette), marker
SOLVERS = [
    ('admm_fp',  r'ADMM-FP',  '#0072B2', 'x'),  # blue
    ('ofp',      r'OFP',      '#CC79A7', '^'),  # reddish purple
    ('admm',     r'ADMM',     '#009E73', 's'),  # bluish green
    ('admm_rnd', r'ADMM-RND', '#E69F00', '*'),  # orange
    ('gurobi',   r'Gurobi',   '#D55E00', 'o'),  # vermillion
    ('scip',     r'SCIP',     '#56B4E9', 'd'),  # sky blue
]

# distinct dash pattern per solver for band_plot: when two series have nearly
# identical values (e.g. ADMM-FP vs ADMM memory use), alpha alone can't
# separate coincident lines, since the topmost one still dominates the shared
# pixels. A differing dash phase lets the other color show through the gaps,
# without displacing either curve from its true value.
LINESTYLES = {
    'admm_fp':  '-',
    'ofp':      '--',
    'admm':     '-.',
    'admm_rnd': ':',
    'gurobi':   (0, (3, 1, 1, 1)),  # densely dash-dotted
    'scip':     (0, (5, 1)),        # densely dashed
}

parser = argparse.ArgumentParser(description='Plot reach-avoid heuristic comparison data.')
parser.add_argument('results_file', nargs='?', default='reach_avoid_results.json',
                     help='Path to the results JSON file (default: reach_avoid_results.json).')
parser.add_argument('--style', choices=['box', 'band'], default='box',
                     help='"band" (default) draws median + IQR bands; "box" draws grouped boxplots.')
args = parser.parse_args()

# load json data
with open(args.results_file, 'r') as f:
    data = json.load(f)

sample_factor_arr = data['sample_factor_arr']
n_seeds = data['n_seeds']

# only keep solvers that are actually present in this results file, so older
# result files that predate a given solver (e.g. ADMM-RND) still plot.
active_solvers = []
for key, label, color, marker in SOLVERS:
    if f'n_feas_{key}_arr' not in data or f't_{key}_arr' not in data:
        print(f'No data for {label} in {args.results_file}; omitting from plots.')
        continue
    active_solvers.append((key, label, color, marker))
n_groups = len(active_solvers)

n_feas_arr = {key: data[f'n_feas_{key}_arr'] for key, *_ in active_solvers}
t_arr = {key: data[f't_{key}_arr'] for key, *_ in active_solvers}
m_arr = {key: data.get(f'm_{key}_arr') for key, *_ in active_solvers}

have_mem = all(m_arr[key] is not None for key, *_ in active_solvers) and \
    any(np.isfinite(mi) for key, *_ in active_solvers for m_sample_factor in m_arr[key] for mi in m_sample_factor)
if not have_mem:
    print(f'No peak memory data in {args.results_file}; plotting solution time and solution rate only.')


# max improvement in convergence rate for ADMM-FP compared to ADMM
d_conv_arr = [n_feas_arr['admm_fp'][i]/n_feas_arr['admm'][i] if n_feas_arr['admm'][i] > 0 else np.inf
              for i in range(len(sample_factor_arr))]
print(f'Max improvement in convergence rate for ADMM-FP compared to ADMM: {np.max(d_conv_arr)}')



### statistics on ADMM-FP vs ADMM solves times ###
sol_time_ratios = np.array([])
for i, sample_factor in enumerate(sample_factor_arr):
    t_admm_fp = np.array(t_arr['admm_fp'][i])
    t_admm = np.array(t_arr['admm'][i])

    # filter out failed solves (represented as NaN/inf)
    valid_indices = np.isfinite(t_admm_fp) & np.isfinite(t_admm)
    t_admm_fp_valid = t_admm_fp[valid_indices]
    t_admm_valid = t_admm[valid_indices]

    if len(t_admm_fp_valid) == 0:
        print(f"Sample factor {sample_factor}: No valid solves for ADMM-FP and ADMM at sample factor {sample_factor}.")
        continue

    sol_time_ratios = np.concatenate((sol_time_ratios, t_admm_valid / t_admm_fp_valid))

print('ADMM-FP sol time / ADMM sol time for cases where both produce feasible solutions:')
print(f'median = {np.median(sol_time_ratios)}, min = {np.min(sol_time_ratios)}, max = {np.max(sol_time_ratios)}')

### statistics on peak memory use ###
if have_mem:
    for key, label, *_ in active_solvers:
        if key in ('admm_fp', 'ofp'):  # not compared against themselves / not included in this stat
            continue
        mem_ratios = np.array([])
        for i, sample_factor in enumerate(sample_factor_arr):
            m_admm_fp = np.array(m_arr['admm_fp'][i])
            m_ref = np.array(m_arr[key][i])

            # only compare cases where both solvers produced a feasible solution
            valid_indices = np.isfinite(m_admm_fp) & np.isfinite(m_ref)
            if not np.any(valid_indices):
                continue

            mem_ratios = np.concatenate((mem_ratios, m_ref[valid_indices] / m_admm_fp[valid_indices]))

        if len(mem_ratios) == 0:
            print(f'No cases where both ADMM-FP and {label} produce feasible solutions.')
            continue

        print(f'{label} peak memory / ADMM-FP peak memory for cases where both produce feasible solutions:')
        print(f'median = {np.median(mem_ratios)}, min = {np.min(mem_ratios)}, max = {np.max(mem_ratios)}')

### plots ###

# generate plots
textwidth_pt = 10.
if is_latex_installed():
    rc_context = {
        "text.usetex": True,
        "font.size": textwidth_pt,
        "font.family": "serif",  # Choose a serif font like 'Times New Roman' or 'Computer Modern'
        "pgf.texsystem": "pdflatex",
        "pgf.rcfonts": False,
    }
else:
    print("LaTeX not installed, using default font.")
    rc_context = {
        "font.size": textwidth_pt,
    }


inches_per_pt = 1 / 72.27

with plt.rc_context(rc_context):

    # solution time / peak memory / rate of solution comparison plot
    figwidth_pt = 505.89
    fig_height_frac = 0.75 if have_mem else 0.5
    figsize = (figwidth_pt * inches_per_pt, fig_height_frac*figwidth_pt * inches_per_pt)  # Convert pt to inches
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    gs = fig.add_gridspec(nrows=3 if have_mem else 2, figure=fig,
                          height_ratios=[2,1,1] if have_mem else [2,1])

    # time steps
    N_arr = [sample_factor * 10 for sample_factor in sample_factor_arr]

    # percent feasible solutions
    r_feas_arr = {key: [float(n_feas)/n_seeds*100. for n_feas in n_feas_arr[key]] for key, *_ in active_solvers}

    # remove non-finite entries from time and memory data before plotting
    for key, *_ in active_solvers:
        t_arr[key] = [[t for t in sample if np.isfinite(t)] for sample in t_arr[key]]
        if have_mem:
            m_arr[key] = [[m for m in sample if np.isfinite(m)] for sample in m_arr[key]]

    # geometry for the grouped boxplot: box_width is derived from n_groups so it
    # adapts automatically as solvers are added/removed, rather than being
    # hardcoded independently of n_groups.
    GROUP_FILL = 0.82       # fraction of each unit-wide horizon group covered by boxes
    BOX_GAP_FRAC = 0.12     # gap between adjacent boxes, as a fraction of box_width
    box_width = GROUP_FILL / (n_groups + (n_groups - 1) * BOX_GAP_FRAC)
    box_pitch = box_width * (1 + BOX_GAP_FRAC)

    group_center_positions = np.arange(len(N_arr))
    positions = []
    for center in group_center_positions:
        start_pos = center - GROUP_FILL / 2 + box_width / 2
        positions.extend(start_pos + i * box_pitch for i in range(n_groups))

    BOX_FILL_SATURATION = 0.65  # 1.0 = full color, 0.0 = white

    def fill_color(color, saturation=BOX_FILL_SATURATION):
        """Blend a color toward white for use as a solid box fill."""
        rgb = np.array(to_rgb(color))
        return tuple(1 - saturation * (1 - rgb))

    def box_face_kwargs(color):
        """Appearance shared by the plotted boxes and their legend swatch,
        so the two cannot visually drift apart."""
        return {'facecolor': fill_color(color), 'edgecolor': color, 'linewidth': 0.6}

    def grouped_boxplot(ax, arrs_by_key, title, ylabel):
        """One boxplot group per horizon, one box per solver within a group."""

        all_data = []
        for i in range(len(N_arr)):
            for key, *_ in active_solvers:
                all_data.append(arrs_by_key[key][i])

        bp = ax.boxplot(
            all_data,
            positions=positions,
            widths=box_width,
            patch_artist=True,
            boxprops={'linewidth': 0.6},
            whiskerprops={'linewidth': 0.6},
            capprops={'linewidth': 0.6},
            capwidths=0.6 * box_width,
            medianprops={'color': 'black', 'linewidth': 0.9},
            whis=(0., 100.)  # whiskers cover all data
        )

        for i, box in enumerate(bp['boxes']):
            _, _, color, _ = active_solvers[i % n_groups]
            box.set(**box_face_kwargs(color))

        # shade alternating horizon groups so each reads as one visual unit
        for i, center in enumerate(group_center_positions):
            if i % 2 == 1:
                ax.axvspan(center - 0.5, center + 0.5, color='0.5', alpha=0.06, zorder=0, lw=0)

        ax.set_title(title, fontsize=textwidth_pt)
        ax.set_ylabel(ylabel, fontsize=textwidth_pt)
        ax.set_yscale('log')
        ax.grid(alpha=0.2)

        ax.set_xticks(group_center_positions)
        ax.set_xticklabels([str(N) for N in N_arr])
        ax.tick_params(axis='x', which='both', length=0)

        return bp

    def band_plot(ax, arrs_by_key, title, ylabel):
        """Median line with a 25th-75th percentile band, one per solver.

        Avoids the narrow-box readability problem entirely, and stays
        readable regardless of how many horizons/solvers are plotted.
        """

        for key, _, color, marker in active_solvers:
            arr = arrs_by_key[key]
            med = [np.median(a) if len(a) else np.nan for a in arr]
            lo = [np.percentile(a, 25) if len(a) else np.nan for a in arr]
            hi = [np.percentile(a, 75) if len(a) else np.nan for a in arr]

            ax.fill_between(N_arr, lo, hi, color=color, alpha=0.15, lw=0)
            # hollow markers: some solvers have only 1-2 valid seeds per horizon
            # (e.g. ADMM often has just one), so their "median" is an isolated
            # point that can land almost exactly on another solver's line/marker.
            # A filled marker would simply hide underneath; an open marker's
            # outline stays visible no matter which one is drawn on top.
            ax.plot(N_arr, med, color=color, marker=marker, linestyle=LINESTYLES[key],
                    markersize=4.5, markerfacecolor='none', markeredgewidth=1.1,
                    linewidth=1.0, alpha=0.5)

        ax.set_title(title, fontsize=textwidth_pt)
        ax.set_ylabel(ylabel, fontsize=textwidth_pt)
        ax.set_yscale('log')
        ax.grid(alpha=0.2)

        ax.set_xticks(N_arr)
        ax.set_xticklabels([str(N) for N in N_arr])

        return ax

    plot_panel = grouped_boxplot if args.style == 'box' else band_plot

    ax = fig.add_subplot(gs[0])
    plot_panel(ax, t_arr, r'Time to find a feasible solution', r'[sec]')

    # single figure-level legend, shared by all panels: color identifies the
    # solver and marker is a redundant, colorblind- and grayscale-safe cue.
    # The legend handle must depict the encoding actually drawn in each style,
    # so it is built per style rather than always using a plain Line2D.
    if args.style == 'box':
        # composite handle: a swatch matching the box exactly (same helper as
        # the boxes themselves, so they cannot visually drift apart) plus the
        # dashed/hollow marker used by the rate panel below.
        handles = [(Patch(**box_face_kwargs(color)),
                    Line2D([0], [0], color=color, marker=marker, linestyle=LINESTYLES[key],
                           markersize=4, markerfacecolor='none', markeredgewidth=1.1))
                   for key, _, color, marker in active_solvers]
        handler_map = {tuple: HandlerTuple(ndivide=None, pad=0.4)}
    else:
        handles = [Line2D([0], [0], color=color, marker=marker, linestyle=LINESTYLES[key],
                           markersize=5, markerfacecolor='none', markeredgewidth=1.1)
                   for key, _, color, marker in active_solvers]
        handler_map = None
    labels = [label for _, label, _, _ in active_solvers]
    fig.legend(handles, labels, loc='outside upper center', ncol=n_groups,
               fontsize=textwidth_pt, handler_map=handler_map)

    # plot peak memory per solve vs sample factor for each solver
    if have_mem:
        ax = fig.add_subplot(gs[1])
        plot_panel(ax, m_arr, r'Peak memory use', r'[MB]')

    # plot rate of solution vs sample factor for each solver
    ax = fig.add_subplot(gs[2] if have_mem else gs[1])
    for key, _, color, marker in active_solvers:
        # dashed lines + hollow markers, matching band_plot: several solvers
        # sit at exactly 100% for low N, so filled markers/solid lines would
        # otherwise stack and hide one another here too.
        ax.plot(N_arr, r_feas_arr[key], color=color, marker=marker, linestyle=LINESTYLES[key],
                markersize=4.5, markerfacecolor='none', markeredgewidth=1.1,
                linewidth=1.0, alpha=0.5)
    ax.set_title(r'Percentage of cases where solution is found', fontsize=textwidth_pt)
    ax.set_ylabel(r'[\%]', fontsize=textwidth_pt)
    ax.grid(which='major', alpha=0.2)
    ax.grid(which='minor', axis='y', alpha=0.1)
    ax.minorticks_on()
    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())

    ax.set_xticks([N for N in N_arr])
    ax.set_xticklabels([str(N) for N in N_arr])
    ax.tick_params(axis='x', which='minor', bottom=False)

    # labels
    fig.supxlabel(r'Planning horizon $N$', fontsize=textwidth_pt)

    # save
    if is_latex_installed():
        plt.savefig(f'motion_planning_fp_heuristic_comp_{args.style}.pgf')

    plt.show()
