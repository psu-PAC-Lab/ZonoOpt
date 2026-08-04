import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import numpy as np
import subprocess
import ast
import json
import sys

def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

# load json data
results_file = sys.argv[1] if len(sys.argv) > 1 else 'reach_avoid_results.json'
with open(results_file, 'r') as f:
    data = json.load(f)

    sample_factor_arr = data['sample_factor_arr']
    n_seeds = data['n_seeds']
    n_feas_admm_fp_arr = data['n_feas_admm_fp_arr']
    n_feas_admm_arr = data['n_feas_admm_arr']
    n_feas_admm_rnd_arr = data['n_feas_admm_rnd_arr']
    n_feas_gurobi_arr = data['n_feas_gurobi_arr']
    n_feas_scip_arr = data['n_feas_scip_arr']
    n_feas_ofp_arr = data['n_feas_ofp_arr']
    t_admm_fp_arr = data['t_admm_fp_arr']
    t_admm_arr = data['t_admm_arr']
    t_admm_rnd_arr = data['t_admm_rnd_arr']
    t_gurobi_arr = data['t_gurobi_arr']
    t_scip_arr = data['t_scip_arr']
    t_ofp_arr = data['t_ofp_arr']
    m_admm_fp_arr = data.get('m_admm_fp_arr')
    m_admm_arr = data.get('m_admm_arr')
    m_admm_rnd_arr = data.get('m_admm_rnd_arr')
    m_gurobi_arr = data.get('m_gurobi_arr')
    m_scip_arr = data.get('m_scip_arr')
    m_ofp_arr = data.get('m_ofp_arr')

    m_arrs = [m_admm_fp_arr, m_admm_arr, m_admm_rnd_arr, m_gurobi_arr, m_scip_arr, m_ofp_arr]
    have_mem = all(m is not None for m in m_arrs) and \
        any(np.isfinite(mi) for m in m_arrs for m_sample_factor in m for mi in m_sample_factor)
    if not have_mem:
        print(f'No peak memory data in {results_file}; plotting solution time and solution rate only.')


# max improvement in convergence rate for ADMM-FP compared to ADMM
d_conv_arr = [n_feas_admm_fp_arr[i]/n_feas_admm_arr[i] if n_feas_admm_arr[i] > 0 else np.inf for i in range(len(sample_factor_arr))]
print(f'Max improvement in convergence rate for ADMM-FP compared to ADMM: {np.max(d_conv_arr)}')



### statistics on ADMM-FP vs ADMM solves times ###
sol_time_ratios = np.array([])
for i, sample_factor in enumerate(sample_factor_arr):
    t_admm_fp = np.array(t_admm_fp_arr[i])
    t_admm = np.array(t_admm_arr[i])

    # filter out failed solves (represented as np.inf)
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
    for label, m_ref_arr in [('ADMM', m_admm_arr), ('ADMM-RND', m_admm_rnd_arr), ('Gurobi', m_gurobi_arr), ('SCIP', m_scip_arr)]:
        mem_ratios = np.array([])
        for i, sample_factor in enumerate(sample_factor_arr):
            m_admm_fp = np.array(m_admm_fp_arr[i])
            m_ref = np.array(m_ref_arr[i])

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
    r_feas_admm_fp_arr = [float(n_feas)/n_seeds*100. for n_feas in n_feas_admm_fp_arr]
    r_feas_admm_arr = [float(n_feas)/n_seeds*100. for n_feas in n_feas_admm_arr]
    r_feas_admm_rnd_arr = [float(n_feas)/n_seeds*100. for n_feas in n_feas_admm_rnd_arr]
    r_feas_gurobi_arr = [float(n_feas)/n_seeds*100. for n_feas in n_feas_gurobi_arr]
    r_feas_scip_arr = [float(n_feas)/n_seeds*100. for n_feas in n_feas_scip_arr]
    r_feas_ofp_arr = [float(n_feas)/n_seeds*100. for n_feas in n_feas_ofp_arr]

    # plot solution time vs sample factor for each solver
    box_width = 0.25
    gap_between_boxes = 0.05
    group_spacing = 0.8
    n_groups = 6
    group_width = n_groups * box_width + (n_groups - 1) * gap_between_boxes

    group_center_positions = np.arange(len(N_arr)) * ((n_groups * box_width) + ((n_groups-1) * gap_between_boxes) + group_spacing)
    positions = []
    for center in group_center_positions:
        start_pos = center - (group_width / 2) + box_width / 2
        group_positions = [start_pos + i * (box_width + gap_between_boxes) for i in range(n_groups)]
        positions.extend(group_positions)

    colors = ['b', 'm', 'g', 'orange', 'r', 'c']
    data_labels = [r'ADMM-FP', r'OFP', r'ADMM', r'ADMM-RND', r'Gurobi', r'SCIP']

    # remove infinities from time and memory data for boxplot
    for i in range(len(N_arr)):
        t_admm_fp_arr[i] = [t for t in t_admm_fp_arr[i] if np.isfinite(t)]
        t_ofp_arr[i] = [t for t in t_ofp_arr[i] if np.isfinite(t)]
        t_admm_arr[i] = [t for t in t_admm_arr[i] if np.isfinite(t)]
        t_admm_rnd_arr[i] = [t for t in t_admm_rnd_arr[i] if np.isfinite(t)]
        t_gurobi_arr[i] = [t for t in t_gurobi_arr[i] if np.isfinite(t)]
        t_scip_arr[i] = [t for t in t_scip_arr[i] if np.isfinite(t)]

        if have_mem:
            m_admm_fp_arr[i] = [m for m in m_admm_fp_arr[i] if np.isfinite(m)]
            m_ofp_arr[i] = [m for m in m_ofp_arr[i] if np.isfinite(m)]
            m_admm_arr[i] = [m for m in m_admm_arr[i] if np.isfinite(m)]
            m_admm_rnd_arr[i] = [m for m in m_admm_rnd_arr[i] if np.isfinite(m)]
            m_gurobi_arr[i] = [m for m in m_gurobi_arr[i] if np.isfinite(m)]
            m_scip_arr[i] = [m for m in m_scip_arr[i] if np.isfinite(m)]

    def grouped_boxplot(ax, arrs, title, ylabel):
        """One boxplot group per horizon, one box per solver within a group."""

        all_data = []
        for i in range(len(N_arr)):
            for arr in arrs:
                all_data.append(arr[i])

        bp = ax.boxplot(
            all_data,
            positions=positions,
            widths=box_width,
            patch_artist=True,
            medianprops={'color': 'black'},
            showfliers=True,
            sym='k.',
            whis=(0., 100.) # whiskers cover all data
        )

        for i, box in enumerate(bp['boxes']):
            box.set_facecolor(colors[i % n_groups])
            box.set_alpha(0.5)

        ax.set_title(title, fontsize=textwidth_pt)
        ax.set_ylabel(ylabel, fontsize=textwidth_pt)
        ax.set_yscale('log')
        ax.grid(alpha=0.2)

        ax.set_xticks(group_center_positions)
        ax.set_xticklabels([str(N) for N in N_arr])

        return bp

    ax = fig.add_subplot(gs[0])
    grouped_boxplot(ax, [t_admm_fp_arr, t_ofp_arr, t_admm_arr, t_admm_rnd_arr, t_gurobi_arr, t_scip_arr],
                    r'Time to find a feasible solution', r'[sec]')

    legend_patches = [plt.Rectangle((0, 0), 1, 1, alpha=0.5, fc=color) for color in colors]

    ax.legend(legend_patches, data_labels, bbox_to_anchor=(0., 1.25, 1., .102), loc=3,
        ncol=n_groups, mode="expand", borderaxespad=0., fontsize=textwidth_pt)

    # plot peak memory per solve vs sample factor for each solver
    if have_mem:
        ax = fig.add_subplot(gs[1])
        grouped_boxplot(ax, [m_admm_fp_arr, m_ofp_arr, m_admm_arr, m_admm_rnd_arr, m_gurobi_arr, m_scip_arr],
                        r'Peak memory use', r'[MB]')

    # plot rate of solution vs sample factor for each solver
    ax = fig.add_subplot(gs[2] if have_mem else gs[1])
    ax.plot(N_arr, r_feas_admm_fp_arr, color='b', marker='x', linestyle='-', alpha=0.5)
    ax.plot(N_arr, r_feas_ofp_arr, color='m', marker='^', linestyle='-', alpha=0.5)
    ax.plot(N_arr, r_feas_admm_arr, color='g', marker='s', linestyle='-', alpha=0.5)
    ax.plot(N_arr, r_feas_admm_rnd_arr, color='orange', marker='*', linestyle='-', alpha=0.5)
    ax.plot(N_arr, r_feas_gurobi_arr, color='r', marker='o', linestyle='-', alpha=0.5)
    ax.plot(N_arr, r_feas_scip_arr, color='c', marker='d', linestyle='-', alpha=0.5)
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
        plt.savefig('motion_planning_fp_heuristic_comp.pgf')

    plt.show()
