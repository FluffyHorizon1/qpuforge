#!/usr/bin/env python3
"""Fill designs/beta1/BETA-1.template.md from designs/beta1/beta1_review.json -> designs/beta1/BETA-1.md"""
import json, sys
r = json.load(open('designs/beta1/beta1_review.json'))
t = open('designs/beta1/BETA-1.template.md').read()
pc = lambda v: f'{100*v:.4f} %'
d, s, ref, lr = r['device'], r['distance_scaling'], r['alpha1_v11_reference'], r['long_run_d5']
fy, w, lay, ec = r['fab_yield_d5'], r['wiring'], r.get('layout', {}), r['exit_criteria']
names = {'lambda': 'Λ(3→5) ≥ 2', 'eps_d5': 'logical error per round at d = 5 ≤ 0.1 %', 'beats_alpha': 'd = 5 at least 3× below Alpha-1 v1.1 at d = 3',
         'lifetime': 'd = 5 logical lifetime > 5× best physical T1', 'yield': 'trimmed fabrication yield ≥ 80 % collision-free',
         'leakage': 'chosen leakage scheme at least as good as the alternative over 50 rounds', 'drc': 'first-cut layout passes placement and island-spacing checks (mask width/space DRC not gated, see §6)'}
gates = '\n'.join(f'| {names.get(k, k)} | {"✔ pass" if v else "✘ fail"} |' for k, v in ec.items())
gates = '| Criterion | Result |\n| --- | --- |\n' + gates
npass = sum(1 for v in ec.values() if v)
bud = r['budget_d3']; sens = r['sensitivity_d3']
budget = ', '.join(f'{k} {100*v:+.0f} %' for k, v in sorted(bud.items(), key=lambda kv: kv[1]))
sn = {'resid_2p5e-4': 'CZ residual 5e-4 → 2.5e-4', 'Q_x2': 'dielectric Q 6e6 → 1.2e7', 'readout_200': 'readout 300 → 200 ns', 'sep_0p25': 'readout separation 0.5 % → 0.25 %'}
sensl = ', '.join(f'{sn.get(k, k)} {100*v:+.0f} %' for k, v in sorted(sens.items(), key=lambda kv: kv[1]))
lrr, lrd = lr['lrr'], lr['lrr_dqlr']
better = lrr['P_L_50'] <= lrd['P_L_50']
dec = ('**Decision confirmed:** the leakage-removing reset alone gives the lower 50-round failure. DQLR clears data-qubit leakage, but its time and gate cost outweigh the benefit at this leakage rate, so it stays an option for hardware that turns out to leak more.'
       if better else '**Decision revisited:** DQLR gives the lower 50-round failure in this run; adopt it as the baseline.')
def drc_text(lay):
    d = lay.get('drc')
    if not d: return 'clean' if lay.get('drc_clean') else 'VIOLATIONS (see beta1_review.txt)'
    pl = 'clean' if d['placement_clean'] else f"{d['placement_problems']} problems"
    ws = 'clean' if d['etch_width_violations'] + d['etch_space_violations'] == 0 else \
         f"**not clean: {d['etch_width_violations']} width and {d['etch_space_violations']} space violations**; the mask is not fabrication-ready"
    return (f"placement and island-spacing checks {pl}. Minimum width/space on the merged etch layer "
            f"({d['etch_min_width_um']:.0f} µm / {d['etch_min_space_um']:.0f} µm, projection metric): {ws}")
y = lambda k: f'{fy[k]["collision_free"]}/{fy[k]["samples"]} ({100*fy[k]["collision_free"]/fy[k]["samples"]:.0f} %)'
vals = {
    'GATE_SUMMARY': f'{npass} of {len(ec)} simulation-gate criteria pass', 'RUNTIME': f'{r["runtime_s"]:.0f}', 'SHOTS': f'{int(r["shots"]):,}',
    'T1_MEAN': f'{d["mean_T1_us"]:.0f}', 'T1_BEST': f'{d["best_T1_us"]:.0f}', 'T2_MEAN': f'{d["mean_T2_us"]:.0f}', 'RO_ERR': f'{100*d["mean_readout_err"]:.2f}',
    'N_CZ': f'{d["couplers"]}', 'TCZ': f'{d["t_cz_mean_ns"]:.1f}', 'CZ_ERR': f'{d["cz_err_mean"]:.2e}', 'CZ_LEAK': f'{d["cz_leak_mean"]:.1e}',
    'EPS3': pc(s['eps_d3']), 'EPS5': pc(s['eps_d5']), 'EPS7': pc(s['eps_d7']), 'A_EPS3': pc(ref['eps_d3']), 'A_EPS5': pc(ref['eps_d5']),
    'N7': f'{int(s.get("shots_d7", 100000)):,}', 'N5_9': f'{int(s.get("shots_d5_9rounds", 50000)):,}', 'EPS5_9': pc(s['eps_d5_9rounds']),
    'L35': f'{s["Lambda_3_5"]:.2f} ± {s["Lambda_3_5"]*((s["eps_err_d3"]/s["eps_d3"])**2+(s["eps_err_d5"]/s["eps_d5"])**2)**0.5:.2f}',
    'L57': f'{s["Lambda_5_7"]:.2f} ± {s["Lambda_5_7"]*((s.get("eps_err_d5_9rounds",0)/s["eps_d5_9rounds"])**2+(s["eps_err_d7"]/max(s["eps_d7"],1e-30))**2)**0.5:.2f}', 'ROUND5': f'{s["round_ns_d5"]:.0f}', 'LIFE5': f'{s["logical_lifetime_d5_ms"]:.2f}',
    'LIFE_X': f'{s["logical_lifetime_d5_ms"]*1e3/d["best_T1_us"]:.1f}', 'BEAT': f'{ref["eps_d3"]/s["eps_d5"]:.1f}',
    'LR_PL': f'{lrr["P_L_50"]:.4f}', 'LR_LD': f'{lrr["leaked_data"]:.1e}', 'LR_LM': f'{lrr["leaked_meas"]:.1e}', 'LR_RN': f'{lrr["round_ns"]:.0f}',
    'LD_PL': f'{lrd["P_L_50"]:.4f}', 'LD_LD': f'{lrd["leaked_data"]:.1e}', 'LD_LM': f'{lrd["leaked_meas"]:.1e}', 'LD_RN': f'{lrd["round_ns"]:.0f}',
    'DQLR_DECISION': dec, 'BUDGET': budget, 'SENS': sensl,
    'Y20': y('spread_2.0'), 'Y10': y('spread_1.0'), 'Y05': y('spread_0.5'),
    'S20': f'{fy["spread_2.0"]["f01_sigma_MHz"]:.0f}', 'S10': f'{fy["spread_1.0"]["f01_sigma_MHz"]:.0f}', 'S05': f'{fy["spread_0.5"]["f01_sigma_MHz"]:.0f}',
    'W_XY': w['xy'], 'W_ZQ': w['flux_qubit'], 'W_ZC': w['flux_coupler'], 'W_RI': w['readout_in'], 'W_RO': w['readout_out'], 'W_TOT': w['total_lines'],
    'FEED': w['feedlines'], 'QPF': f'{49/w["feedlines"]:.0f}',
    'CHIP_W': f'{lay.get("width_mm", 0):.1f}', 'CHIP_H': f'{lay.get("height_mm", 0):.1f}', 'ETCH': f'{lay.get("etch_polygons", 0):,}', 'BUMPS': lay.get('bumps', 0),
    'DRC': drc_text(lay), 'GATES': gates,
}
for k, v in vals.items(): t = t.replace('{{' + k + '}}', str(v))
left = [x for x in t.split('{{')[1:]]
if left: print('unfilled placeholders:', [x.split('}}')[0] for x in left]); sys.exit(1)
open('designs/beta1/BETA-1.md', 'w').write(t)
print('wrote designs/beta1/BETA-1.md;', vals['GATE_SUMMARY'])
