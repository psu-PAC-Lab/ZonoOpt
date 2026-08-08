import argparse
import matplotlib.pyplot as plt
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
    ('gurobi',   r'Gurobi',   '#D55E00', 'o'),  # vermillion
    ('scip',     r'SCIP',     '#56B4E9', 'd'),  # sky blue
]

parser = argparse.ArgumentParser(description='Plot random hybrid zonotope heuristic comparison data.')
parser.add_argument('results_file', nargs='?', default='heuristic_comparison_results.json',
                     help='Path to the results JSON file (default: heuristic_comparison_results.json).')
args = parser.parse_args()

# load json data
with open(args.results_file, 'r') as f:
    data = json.load(f)

n_trials = data['n_trials']

# only keep solvers that are actually present in this results file, so older
# result files that predate a given solver still plot.
active_solvers = []
for key, label, color, marker in SOLVERS:
    if key not in data.get('t_arr', {}):
        print(f'No data for {label} in {args.results_file}; omitting from plots.')
        continue
    active_solvers.append((key, label, color, marker))
n_groups = len(active_solvers)

t_arr = {key: data['t_arr'][key] for key, *_ in active_solvers}
rel_obj_arr = {key: data['rel_obj_arr'][key] for key, *_ in active_solvers}
iter_arr = {key: data['iter_arr'][key] for key, *_ in active_solvers}

# print number of cases where solution is found
print(f'Out of {n_trials} trials:')
for key, label, *_ in active_solvers:
    print(f'  {label} found feasible solution in {len(t_arr[key])} cases.')

# plot results

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

    colors = [color for _, _, color, _ in active_solvers]
    data_labels = [label for _, label, _, _ in active_solvers]

    # solution time / rate of solution comparison plot
    figwidth_pt = 245.71
    figsize = (figwidth_pt * inches_per_pt, 1.3*figwidth_pt * inches_per_pt)  # Convert pt to inches
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    gs = fig.add_gridspec(ncols=1, nrows=4, height_ratios=[1, 2, 2, 2])

    ax = fig.add_subplot(gs[0,0])
    ax.bar(range(n_groups),
           [len(t_arr[key])*(100./n_trials) for key, *_ in active_solvers],
           color=colors,
           alpha=0.5)
    ax.set_xticklabels([])
    ax.set_ylabel(r'[\%]', fontsize=textwidth_pt)
    ax.set_title(r'Percent trials with feasible solution found', fontsize=textwidth_pt)
    ax.grid(axis='y', which='major', alpha=0.2)


    # suboptimality of found solution
    ax = fig.add_subplot(gs[1,0])

    all_data = [rel_obj_arr[key] for key, *_ in active_solvers]
    bp = ax.boxplot(
        all_data,
        patch_artist=True,
        medianprops={'color': 'black'},
        showfliers=True,
        sym='k.',
        whis=(0., 100.) # whiskers cover all data
    )

    for i, box in enumerate(bp['boxes']):
        box.set_facecolor(colors[i])
        box.set_alpha(0.5)

    ax.set_title(r'Suboptimality of found solution', fontsize=textwidth_pt)
    ax.set_ylabel(r'[\%]', fontsize=textwidth_pt)

    ax.grid(axis='y', which='major', alpha=0.2)

    ax.set_xticklabels([])


    # time to find a feasible solution
    ax = fig.add_subplot(gs[2,0])

    all_data = [t_arr[key] for key, *_ in active_solvers]
    bp = ax.boxplot(
        all_data,
        patch_artist=True,
        medianprops={'color': 'black'},
        showfliers=True,
        sym='k.',
        whis=(0., 100.) # whiskers cover all data
    )

    for i, box in enumerate(bp['boxes']):
        box.set_facecolor(colors[i])
        box.set_alpha(0.5)

    ax.set_title(r'Time to find a feasible solution', fontsize=textwidth_pt)
    ax.set_ylabel(r'[sec]', fontsize=textwidth_pt)
    ax.set_yscale('log')

    ax.grid(axis='y', which='major', alpha=0.2)

    ax.set_xticklabels([])


    # number of iterations to find a feasible solution
    ax = fig.add_subplot(gs[3,0])

    all_data = [iter_arr[key] for key, *_ in active_solvers]
    bp = ax.boxplot(
        all_data,
        patch_artist=True,
        medianprops={'color': 'black'},
        showfliers=True,
        sym='k.',
        whis=(0., 100.), # whiskers cover all data
        tick_labels=data_labels
    )

    for i, box in enumerate(bp['boxes']):
        box.set_facecolor(colors[i])
        box.set_alpha(0.5)

    ax.set_title(r'Number of iterations', fontsize=textwidth_pt)
    ax.set_ylabel(r'[num]', fontsize=textwidth_pt)
    ax.set_yscale('log')

    ax.grid(axis='y', which='major', alpha=0.2)

    ax.tick_params(axis='x', which='minor', bottom=False)

    if is_latex_installed():
        plt.savefig('heuristic_comparison.pgf')

    plt.show()
