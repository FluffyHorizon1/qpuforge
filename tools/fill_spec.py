#!/usr/bin/env python3
"""Fill a design spec from its review data.

    python3 tools/fill_spec.py beta  [DIR]   designs/beta1/BETA-1.template.md   + beta1_review.json   -> BETA-1.md
    python3 tools/fill_spec.py gamma [DIR]   designs/gamma1/GAMMA-1.template.md + gamma1_review.json  -> GAMMA-1.md

Every number in the spec comes from the review JSON; the script stops if a placeholder is left unfilled.
Run from the repository root. With DIR, the review JSON is read from DIR and the spec is written there (the directory
given to a review program with --out); the committed files stay untouched.
"""
import json
import os
import sys


def load(path):
    """Review JSON; a non-finite number (written as null, e.g. a ratio with no failures in a short run) reads as NaN."""
    def fix(v):
        if isinstance(v, dict):
            return {k: fix(x) for k, x in v.items()}
        if isinstance(v, list):
            return [fix(x) for x in v]
        return float('nan') if v is None else v
    with open(path, encoding='utf-8') as f:
        return fix(json.load(f))


def div(a, b):
    return a / b if b else float('nan')


def paths(design_dir, review, spec, scratch):
    """(review JSON, output spec): the committed locations, or both inside the scratch directory."""
    base = scratch if scratch else design_dir
    return os.path.join(base, review), os.path.join(base, spec)


def pc(v, digits=4):
    return f'{100 * v:.{digits}f} %'


def gate_table(ec, names):
    rows = '\n'.join(f'| {names.get(k, k)} | {"✔ pass" if v else "✘ fail"} |' for k, v in ec.items())
    return '| Criterion | Result |\n| --- | --- |\n' + rows, sum(1 for v in ec.values() if v)


def fill(template, out, vals):
    t = open(template, encoding='utf-8').read()
    for k, v in vals.items():
        t = t.replace('{{' + k + '}}', str(v))
    left = [x.split('}}')[0] for x in t.split('{{')[1:]]
    if left:
        print('unfilled placeholders:', left)
        sys.exit(1)
    open(out, 'w', encoding='utf-8').write(t)
    print(f'wrote {out}; {vals["GATE_SUMMARY"]}')


def rel_err(a, b):
    """Relative one-sigma error of a ratio from (value, error) pairs."""
    return ((a[1] / max(a[0], 1e-30)) ** 2 + (b[1] / max(b[0], 1e-30)) ** 2) ** 0.5


def drc_text(lay):
    d = lay.get('drc')
    if not d:
        return 'clean' if lay.get('drc_clean') else 'VIOLATIONS (see the review output)'
    pl = 'clean' if d['placement_clean'] else f"{d['placement_problems']} problems"
    ws = ('clean' if d['etch_width_violations'] + d['etch_space_violations'] == 0 else
          f"**not clean: {d['etch_width_violations']:,} width and {d['etch_space_violations']:,} space violations**; the mask is not fabrication-ready")
    return (f"placement and island-spacing checks {pl}. Minimum width/space on the merged etch layer "
            f"({d['etch_min_width_um']:.0f} µm / {d['etch_min_space_um']:.0f} µm, projection metric): {ws}")


def mask_clean(lay):
    d = lay.get('drc', {})
    return bool(d) and d.get('placement_clean') and d['etch_width_violations'] + d['etch_space_violations'] == 0


def mask_note(lay, sec):
    d = lay.get('drc', {})
    if mask_clean(lay):
        return (f"**Mask DRC.** The layout passes the placement and island-spacing checks and the mask-level check "
                f"(minimum width and space of {d['etch_min_width_um']:.0f} µm on the merged etch layer, projection metric) "
                f"with no violations, after {d.get('ground_sliver_fills', 0):,} ground slivers narrower than the minimum were merged "
                f"into the neighbouring cut-outs (`layout.hpp`, `fill_ground_slivers`; the extraction cells get the same clean-up). "
                f"Corner-to-corner (Euclidean) spacing is not checked: the foundry's DRC deck must still be run before tape-out.")
    return (f"**The mask is not fabrication-ready.** The mask-level check (minimum width and space on the merged etch layer) "
            f"reports {d.get('etch_width_violations', 0):,} width and {d.get('etch_space_violations', 0):,} space violations (§{sec}).")


def leak_note(d):
    mx, mean = d['cz_leak_max'], d['cz_leak_mean']
    if mx > 2e-4:
        return (f'The largest per-pair leakage, {mx:.1e} per gate, is above the 2e-4 limit of the Alpha simulation gate '
                f'(ROADMAP.md); the mean is {mean:.1e}.')
    return f'The largest per-pair leakage, {mx:.1e} per gate, is within the 2e-4 limit of the Alpha simulation gate (ROADMAP.md).'


def layout_vals(lay):
    d = lay.get('drc', {})
    return {
        'CHIP_W': f'{lay.get("width_mm", 0):.1f}', 'CHIP_H': f'{lay.get("height_mm", 0):.1f}', 'ETCH': f'{lay.get("etch_polygons", 0):,}',
        'JJ': lay.get('junction_regions', 0), 'BUMPS': lay.get('bumps', 0),
        'DRC': drc_text(lay), 'MASK_STATUS': 'passes the mask width/space DRC' if mask_clean(lay) else 'first cut, not mask-ready', 'DRC_W': f"{d.get('etch_width_violations', 0):,}", 'DRC_S': f"{d.get('etch_space_violations', 0):,}",
        'EC_ERR': f'{100 * lay.get("max_EC_err", 0):.1f} %', 'G_ERR': f'{100 * lay.get("max_coupling_err", 0):.1f} %',
        'RO_C_ERR': f'{100 * lay.get("readout_C_err", 0):+.1f} %', 'RO_C': f'{lay.get("readout_C_fF", 0):.2f}',
        'RO_G': f'{lay.get("readout_g_MHz", 0):.0f}', 'RO_G_T': f'{lay.get("readout_g_target_MHz", 0):.0f}',
        'PARK': f'{lay.get("park_fc_GHz", 0):.3f}', 'PARK_ZZ': f'{lay.get("park_zz_kHz", 0):+.2f}',
    }


def beta(scratch=None):
    src, dst = paths('designs/beta1', 'beta1_review.json', 'BETA-1.md', scratch)
    r = load(src)
    d, s, ref, lr, st = r['device'], r['distance_scaling'], r['alpha1_v11_reference'], r['long_run_d5'], r['stability_100_rounds_d5']
    fy, w, lay, ec = r['fab_yield_d5'], r['wiring'], r.get('layout', {}), r['exit_criteria']
    names = {'lambda': 'Λ(3→5) ≥ 2 in both bases', 'eps_d5': 'logical error per round at d = 5 ≤ 0.1 % in both bases',
             'beats_alpha': 'd = 5 at least 3× below Alpha-1 v1.1 at d = 3', 'lifetime': 'd = 5 logical lifetime > 5× best physical T1',
             'yield': 'trimmed fabrication yield ≥ 80 % collision-free',
             'leakage': 'chosen leakage scheme at least as good as the alternative over 50 rounds',
             'stability': 'detection rate flat over 100 rounds with the leakage-removing reset (drift < 5 %)',
             'drc': 'layout passes placement, island-spacing and mask width/space DRC (4 µm / 4 µm on the merged etch layer, see §6)',
             'layout': 'extracted E_C within 3 %, couplings and readout capacitance within 6 %; parked exchange re-nullable'}
    gates, npass = gate_table(ec, names)
    bud, sens = r['budget_d3'], r['sensitivity_d3']
    budget = ', '.join(f'{k} {100 * v:+.0f} %' for k, v in sorted(bud.items(), key=lambda kv: kv[1]))
    sn = {'resid_2p5e-4': 'CZ residual 5e-4 → 2.5e-4', 'Q_x2': 'dielectric Q 6e6 → 1.2e7', 'readout_200': 'readout 300 → 200 ns',
          'sep_0p25': 'readout separation 0.5 % → 0.25 %'}
    sensl = ', '.join(f'{sn.get(k, k)} {100 * v:+.0f} %' for k, v in sorted(sens.items(), key=lambda kv: kv[1]))
    lrr, lrd = lr['lrr'], lr['lrr_dqlr']
    n_long = r['shots'] / 2                                    # the review runs the 50-round experiments at half the shots
    sem = lambda q: (max(q * (1 - q), 1e-12) / n_long) ** 0.5
    zscore = abs(lrd['P_L_50'] - lrr['P_L_50']) / (sem(lrr['P_L_50']) ** 2 + sem(lrd['P_L_50']) ** 2) ** 0.5
    if lrr['P_L_50'] <= lrd['P_L_50'] * 1.05:
        dec = ('**Decision confirmed:** the leakage-removing reset alone does at least as well over 50 rounds '
               f'(the two differ by {zscore:.1f} standard errors). DQLR clears data-qubit leakage, but its time and gate '
               'cost are not repaid at this leakage rate, so it stays an option.')
    else:
        dec = (f'**Decision revisited:** DQLR gives the lower 50-round failure in this run, by {zscore:.1f} standard errors. '
               'Re-test before adopting it as the baseline.')
    y = lambda k: f'{fy[k]["collision_free"]}/{fy[k]["samples"]} ({100 * fy[k]["collision_free"] / fy[k]["samples"]:.0f} %)'
    stl, stp = st['lrr'], st['plain_reset']
    dc = r['decoding_d5']
    # the "Alpha-1 v1.1 sensitivity" column of the design table quotes the committed Alpha review
    asens = load('designs/alpha1-v11/alpha1_v11_review.json')['sensitivity']
    sg = lambda v: f'{"+" if v >= 0 else "−"}{abs(100 * v):.0f} %'
    vals = {
        'GATE_SUMMARY': f'{npass} of {len(ec)} simulation-gate criteria pass', 'RUNTIME': f'{r["runtime_s"]:.0f}', 'SHOTS': f'{int(r["shots"]):,}',
        'T1_MEAN': f'{d["mean_T1_us"]:.0f}', 'T1_BEST': f'{d["best_T1_us"]:.0f}', 'T2_MEAN': f'{d["mean_T2_us"]:.0f}', 'RO_ERR': f'{100 * d["mean_readout_err"]:.2f}',
        'N_CZ': f'{d["couplers"]}', 'TCZ': f'{d["t_cz_mean_ns"]:.1f}', 'CZ_ERR': f'{d["cz_err_mean"]:.2e}', 'CZ_LEAK': f'{d["cz_leak_mean"]:.1e}', 'CZ_LEAK_MAX': f'{d["cz_leak_max"]:.1e}', 'CZ_LEAK_NOTE': leak_note(d),
        'EPS3': pc(s['eps_d3']), 'EPS5': pc(s['eps_d5']), 'EPS7': pc(s['eps_d7']), 'EPS3X': pc(s['eps_d3_X']), 'EPS5X': pc(s['eps_d5_X']),
        'A_EPS3': pc(ref['eps_d3']), 'A_EPS5': pc(ref['eps_d5']), 'N7': f'{int(s["shots_d7"]):,}', 'F7': f'{int(s["d7_failures"])}',
        'N5_9': f'{int(s["shots_d5_9rounds"]):,}', 'EPS5_9': pc(s['eps_d5_9rounds']),
        'L35': f'{s["Lambda_3_5"]:.2f} ± {s["Lambda_3_5"] * rel_err((s["eps_d3"], s["eps_err_d3"]), (s["eps_d5"], s["eps_err_d5"])):.2f}',
        'L35X': f'{s["Lambda_3_5_X"]:.2f} ± {s["Lambda_3_5_X"] * rel_err((s["eps_d3_X"], s["eps_err_d3_X"]), (s["eps_d5_X"], s["eps_err_d5_X"])):.2f}',
        'L57': f'{s["Lambda_5_7"]:.2f} ± {s["Lambda_5_7"] * rel_err((s["eps_d5_9rounds"], s["eps_err_d5_9rounds"]), (s["eps_d7"], s["eps_err_d7"])):.2f}',
        'ROUND5': f'{s["round_ns_d5"]:.0f}',
        'LIFE5': f'{s["logical_lifetime_d5_ms"]:.2f}', 'LIFE_X': f'{s["logical_lifetime_d5_ms"] * 1e3 / d["best_T1_us"]:.1f}',
        'BEAT': f'{div(ref["eps_d3"], s["eps_d5"]):.1f}',
        'LR_PL': f'{lrr["P_L_50"]:.4f}', 'LR_LD': f'{lrr["leaked_data"]:.1e}', 'LR_LM': f'{lrr["leaked_meas"]:.1e}', 'LR_RN': f'{lrr["round_ns"]:.0f}',
        'LD_PL': f'{lrd["P_L_50"]:.4f}', 'LD_LD': f'{lrd["leaked_data"]:.1e}', 'LD_LM': f'{lrd["leaked_meas"]:.1e}', 'LD_RN': f'{lrd["round_ns"]:.0f}',
        'DQLR_DECISION': dec, 'BUDGET': budget, 'SENS': sensl,
        'A_SENS_Q': sg(asens['Q_x2']), 'A_SENS_T': sg(asens['Teff_35']), 'A_SENS_RO': sg(asens['readout_300']), 'A_SENS_CZ': sg(asens['resid_0p5e-3']),
        'DEC_R': int(dc['rounds']), 'DEC_N': f'{int(dc["shots"]):,}', 'DEC_DEF': f'{dc["mean_defects"]:.1f}', 'DEC_OVER': f'{100 * dc["frac_over_22"]:.0f} %',
        'DEC_PL_EX': f'{100 * dc["P_L_exact"]:.1f} %', 'DEC_PL_LG': f'{100 * dc["P_L_legacy"]:.1f} %',
        'DEC_F_EX': int(dc['failures_exact']), 'DEC_F_LG': int(dc['failures_legacy']),
        'ST_L_DRIFT': f'{100 * stl["drift"]:+.1f} %', 'ST_L_PL': f'{stl["P_L_100"]:.4f}', 'ST_L_DEF': f'{stl["mean_defects"]:.1f} ({int(stl["max_defects"])})',
        'ST_P_DRIFT': f'{100 * stp["drift"]:+.1f} %', 'ST_P_PL': f'{stp["P_L_100"]:.4f}', 'ST_P_DEF': f'{stp["mean_defects"]:.1f} ({int(stp["max_defects"])})',
        'Y20': y('spread_2.0'), 'Y10': y('spread_1.0'), 'Y05': y('spread_0.5'),
        'S20': f'{fy["spread_2.0"]["f01_sigma_MHz"]:.0f}', 'S10': f'{fy["spread_1.0"]["f01_sigma_MHz"]:.0f}', 'S05': f'{fy["spread_0.5"]["f01_sigma_MHz"]:.0f}',
        'W_XY': w['xy'], 'W_ZQ': w['flux_qubit'], 'W_ZC': w['flux_coupler'], 'W_RI': w['readout_in'], 'W_RO': w['readout_out'], 'W_TOT': w['total_lines'],
        'FEED': w['feedlines'], 'QPF': f'{w["qubits_per_feedline"]:.1f}', 'RO_LO': f'{w["readout_band_lo_GHz"]:.2f}', 'RO_HI': f'{w["readout_band_hi_GHz"]:.2f}',
        'EXTRA_LINES': 2 * (w['feedlines'] - 3), 'PARK_BELOW': f'{(8.5 - lay.get("park_fc_GHz", 8.5)) * 1e3:.0f}',
        'GATES': gates,
    }
    vals.update(layout_vals(lay)); vals['MASK_NOTE'] = mask_note(lay, '6')
    fill('designs/beta1/BETA-1.template.md', dst, vals)


def gamma(scratch=None):
    src, dst = paths('designs/gamma1', 'gamma1_review.json', 'GAMMA-1.md', scratch)
    r = load(src)
    d, s, dec, imp, lr = r['device'], r['distance_scaling'], r['decoding'], r['impacts'], r['long_run_d7']
    fy, w, lay, ec = r['fab_yield_d7'], r['wiring'], r.get('layout', {}), r['exit_criteria']
    names = {'lambda': 'Λ(5→7) ≥ 2 in both bases', 'eps_d7': 'distance-7 logical error per round ≤ 0.01 % in both bases',
             'lifetime': 'distance-7 logical lifetime > 50× best physical T1',
             'impacts': 'logical-error floor from impacts below 10 % of the distance-7 error (both scenarios)',
             'stability': 'detection rate flat over 40 rounds with the leakage-removing reset (drift < 5 %)',
             'yield': 'trimmed fabrication yield ≥ 80 % collision-free',
             'drc': 'layout passes placement, island-spacing and mask width/space DRC (4 µm / 4 µm on the merged etch layer, see §7)',
             'layout': 'extracted E_C within 3 %, couplings and readout capacitance within 6 %; parked exchange re-nullable'}
    gates, npass = gate_table(ec, names)
    e = lambda dd, b: pc(s[f'd{dd}'][f'eps_{b}'], 5 if dd == 7 else 4)
    f = lambda dd: (f'{int(s[f"d{dd}"]["failures_Z"])} / {int(s[f"d{dd}"]["shots_Z"]):,}; '
                    f'{int(s[f"d{dd}"]["failures_X"])} / {int(s[f"d{dd}"]["shots_X"]):,}')
    y = lambda k: f'{fy[k]["collision_free"]}/{fy[k]["samples"]} ({100 * fy[k]["collision_free"] / fy[k]["samples"]:.0f} %)'
    d7 = s['d7']
    rel = max(d7['eps_Z_err'] / max(d7['eps_Z'], 1e-30), d7['eps_X_err'] / max(d7['eps_X'], 1e-30))
    lam = lambda a, b, bs: (f'{s[f"Lambda_{a}_{b}_{bs}"]:.2f} ± '
                            f'{s[f"Lambda_{a}_{b}_{bs}"] * rel_err((s[f"d{a}"][f"eps_{bs}"], s[f"d{a}"][f"eps_{bs}_err"]), (s[f"d{b}"][f"eps_{bs}"], s[f"d{b}"][f"eps_{bs}_err"])):.2f}')
    b1, b2 = imp['unmitigated'], imp['gap_engineered']
    lg = dec['long_run']
    lrr, pr = lr['lrr'], lr['plain_reset']
    share = lambda b: (f'{100 * b["floor_over_eps_d7"]:.2f} %')
    vals = {
        'GATE_SUMMARY': f'{npass} of {len(ec)} memory-gate criteria pass', 'RUNTIME': f'{r["runtime_s"]:.0f}', 'SHOTS': f'{int(r["shots"]):,}',
        'T1_MEAN': f'{d["mean_T1_us"]:.0f}', 'T1_BEST': f'{d["best_T1_us"]:.0f}', 'T2_MEAN': f'{d["mean_T2_us"]:.0f}', 'RO_ERR': f'{100 * d["mean_readout_err"]:.2f}',
        'N_CZ': f'{d["couplers"]}', 'TCZ': f'{d["t_cz_mean_ns"]:.1f}', 'CZ_ERR': f'{d["cz_err_mean"]:.2e}', 'CZ_LEAK': f'{d["cz_leak_mean"]:.1e}',
        'CZ_LEAK_MAX': f'{d["cz_leak_max"]:.1e}', 'CZ_LEAK_NOTE': leak_note(d), 'CZ_VERSION': r['provenance']['cz_model_version'],
        'COLLISIONS': int(d['collisions']), 'ROUNDS': int(s['rounds']),
        'EZ3': e(3, 'Z'), 'EZ5': e(5, 'Z'), 'EZ7': e(7, 'Z'), 'EX3': e(3, 'X'), 'EX5': e(5, 'X'), 'EX7': e(7, 'X'),
        'F3': f(3), 'F5': f(5), 'F7': f(7),
        'L57Z': lam(5, 7, 'Z'), 'L57X': lam(5, 7, 'X'), 'L35Z': lam(3, 5, 'Z'), 'L35X': lam(3, 5, 'X'),
        'ROUND7': f'{d7["round_ns"]:.0f}', 'ROUND7_US': f'{d7["round_ns"] / 1e3:.2f}', 'LIFE7': f'{s["logical_lifetime_d7_ms"]:.1f}',
        'LIFE_X': f'{s["lifetime_over_best_T1"]:.0f}', 'E7ERR': f'{100 * rel:.0f} %',
        'DEF_MEAN': f'{dec["mean_defects_d7"]:.1f}', 'DEF_MAX': int(dec['max_defects_d7']), 'DEF_OVER': f'{100 * dec["frac_over_22_d7"]:.2f} %',
        'DEC_US': f'{dec["decode_us_per_round_d7"]:.2f}',
        'LONG_R': int(lg['rounds']), 'LONG_N': f'{int(lg["shots"]):,}', 'LONG_DEF': f'{lg["mean_defects"]:.1f}', 'LONG_OVER': f'{100 * lg["frac_over_22"]:.0f} %',
        'LONG_PL_EX': f'{lg["P_L_exact"]:.4f}', 'LONG_PL_LG': f'{lg["P_L_legacy"]:.4f}', 'LONG_F_EX': int(lg['failures_exact']), 'LONG_F_LG': int(lg['failures_legacy']),
        'LONG_US_EX': f'{lg["decode_us_per_round_exact"]:.1f}', 'LONG_US_LG': f'{lg["decode_us_per_round_legacy"]:.0f}',
        'B1_RATE': f'{b1["rate_hz"]:g}', 'B1_TAU': f'{b1["tau_ms"]:g}', 'B1_T1': f'{b1["t1_peak_us"]:g}', 'B1_PF': f'{b1["p_fail"]:.2f}',
        'B1_DEAD': f'{b1["dead_time_ms"]:.0f}', 'B1_VIS': f'{b1["detectable_ms"]:.0f}', 'B1_FLOOR': f'{b1["floor_per_round"]:.1e}', 'B1_SHARE': share(b1),
        'B1_DSTAR': f'{b1["distance_where_impacts_dominate"]:.0f}', 'B1_LOST': f'{100 * b1["lost_time_fraction"]:.1f} %',
        'B2_RATE': f'{b2["rate_hz"]:g}', 'B2_TAU': f'{b2["tau_ms"]:g}', 'B2_TP': f'{b2["tphi_peak_us"]:g}', 'B2_PF': f'{b2["p_fail"]:.2f}',
        'B2_DEAD': f'{b2["dead_time_ms"]:.1f}', 'B2_VIS': f'{b2["detectable_ms"]:.1f}', 'B2_FLOOR': f'{b2["floor_per_round"]:.1e}', 'B2_SHARE': share(b2),
        'LR_F': f'{int(lrr["failures"])} / {int(lrr["shots"]):,}', 'PR_F': f'{int(pr["failures"])} / {int(pr["shots"]):,}',
        'PR_RATIO': f'{div(pr["P_L_40"], lrr["P_L_40"]):.0f}',
        'LR_PL': f'{lrr["P_L_40"]:.5f}', 'LR_D0': f'{lrr["detection_by_round"][1]:.4f}', 'LR_D1': f'{lrr["detection_by_round"][-2]:.4f}',
        'LR_DRIFT': f'{100 * lrr["drift"]:+.1f} %', 'LR_LD': f'{lrr["leaked_data"]:.1e}', 'LR_LM': f'{lrr["leaked_meas"]:.1e}',
        'PR_PL': f'{pr["P_L_40"]:.5f}', 'PR_D0': f'{pr["detection_by_round"][1]:.4f}', 'PR_D1': f'{pr["detection_by_round"][-2]:.4f}',
        'PR_DRIFT': f'{100 * pr["drift"]:+.1f} %', 'PR_LD': f'{pr["leaked_data"]:.1e}', 'PR_LM': f'{pr["leaked_meas"]:.1e}',
        'Y10': y('spread_1.0'), 'Y05': y('spread_0.5'),
        'S10': f'{fy["spread_1.0"]["f01_sigma_MHz"]:.0f}', 'S05': f'{fy["spread_0.5"]["f01_sigma_MHz"]:.0f}',
        'W_XY': w['xy'], 'W_ZQ': w['flux_qubit'], 'W_ZC': w['flux_coupler'], 'W_RI': w['readout_in'], 'W_RO': w['readout_out'], 'W_TOT': w['total_lines'],
        'FEED': w['feedlines'], 'QPF': f'{w["qubits_per_feedline"]:.1f}', 'RO_LO': f'{w["readout_band_lo_GHz"]:.2f}', 'RO_HI': f'{w["readout_band_hi_GHz"]:.2f}',
        'BETA_LINES': w['beta1_lines'], 'ALPHA_LINES': w['alpha1_lines'],
        'GATES': gates,
    }
    vals.update(layout_vals(lay)); vals['MASK_NOTE'] = mask_note(lay, '7')
    fill('designs/gamma1/GAMMA-1.template.md', dst, vals)


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else ''
    scratch = sys.argv[2] if len(sys.argv) > 2 else None
    if which == 'beta':
        beta(scratch)
    elif which == 'gamma':
        gamma(scratch)
    else:
        print(__doc__)
        sys.exit(2)
