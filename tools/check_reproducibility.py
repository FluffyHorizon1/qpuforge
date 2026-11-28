#!/usr/bin/env python3
"""Results-policy checks for committed outputs (docs/RESULTS_POLICY.md).

    python3 tools/check_reproducibility.py policy                       size limits + provenance coverage (no build needed)
    python3 tools/check_reproducibility.py regen --tier fast --build B  regenerate the cheap deterministic outputs, diff
    python3 tools/check_reproducibility.py regen --tier full --build B  also rerun the three review programs (~30 min)
    python3 tools/check_reproducibility.py compare ROOT_A ROOT_B        diff the committed outputs of two checkouts

Run from the repository root. B is a built tree (Release, any compiler, any OMP_NUM_THREADS). Comparison is bit-for-bit
after masking only what is measured rather than computed:
  * JSON:  `runtime_s`, `date_utc`, `decode_us_per_round*`, and the provenance block except its seeds and CZ model version.
  * text:  the "review complete in N s" line and decode-time figures; lines are compared as a sorted multiset because one
           banner line is written on stderr and its position relative to stdout varies with buffering.
  * GDS:   the BGNLIB / BGNSTR date fields (set from SOURCE_DATE_EPOCH or the clock), nothing else.
Everything else must match byte for byte. Exit status 1 on any mismatch; a sha256 of the masked content is printed for
both sides so a log shows exactly what was compared.
"""
import argparse
import fnmatch
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import time

POLICY = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_policy.json')

# ---------------------------------------------------------------- masking

VOLATILE_JSON = re.compile(r'^(runtime_s|date_utc|decode_us_per_round.*)$')
KEEP_PROVENANCE = ('seeds', 'cz_model_version')


def _mask_json(o, top=True):
    if isinstance(o, dict):
        out = {}
        for k, v in o.items():
            if VOLATILE_JSON.match(k):
                continue
            if k == 'provenance' and isinstance(v, dict):
                out[k] = {x: v[x] for x in KEEP_PROVENANCE if x in v}
                continue
            out[k] = _mask_json(v, False)
        return out
    if isinstance(o, list):
        return [_mask_json(x, False) for x in o]
    return o


def norm_json(data):
    return json.dumps(_mask_json(json.loads(data)), sort_keys=True, separators=(',', ':')).encode()


def norm_text(data):
    lines = []
    for ln in data.decode('utf-8', 'replace').splitlines():
        ln = re.sub(r'review complete in \d+ s', 'review complete in N s', ln)
        if 'decode' in ln:
            ln = re.sub(r'[0-9.]+ us per round', 'N us per round', ln)
        lines.append(ln)
    return '\n'.join(sorted(lines)).encode()


def norm_gds(data):
    out, i = bytearray(data), 0
    while i + 4 <= len(data):
        n, rec = struct.unpack('>HB', data[i:i + 3])
        if n < 4:
            break
        if rec in (0x01, 0x05):                       # BGNLIB, BGNSTR: modification + access dates
            out[i + 4:i + n] = b'\0' * (n - 4)
        i += n
    return bytes(out)


def norm_exact(data):
    return data


def normaliser(path):
    if path.endswith('.gds'):
        return norm_gds
    if path.endswith('.json'):
        return norm_json
    if path.endswith(('.txt',)):
        return norm_text
    return norm_exact


def digest(b):
    return hashlib.sha256(b).hexdigest()


def compare_files(label, committed, fresh, mode=None):
    """(ok, line): committed vs regenerated, masked per file type (`mode` overrides, 'exact' = no masking)."""
    norm = norm_exact if mode == 'exact' else normaliser(committed)
    with open(committed, 'rb') as f:
        a = f.read()
    with open(fresh, 'rb') as f:
        b = f.read()
    da, db = digest(norm(a)), digest(norm(b))
    ok = da == db
    raw = 'raw-identical' if a == b else 'raw-differs(masked fields)'
    print(f"  {'OK  ' if ok else 'FAIL'} {label:<44} {da[:16]} {db[:16]}  {raw}")
    return ok

# ---------------------------------------------------------------- policy

def tracked(root):
    out = subprocess.run(['git', 'ls-files', '-z'], cwd=root, check=True, capture_output=True).stdout
    return [p for p in out.decode().split('\0') if p]


def expand_braces(s):
    m = re.search(r'\{([^{}]*)\}', s)
    if not m:
        return [s]
    res = []
    for alt in m.group(1).split(','):
        res += expand_braces(s[:m.start()] + alt + s[m.end():])
    return res


def check_policy(root):
    with open(POLICY) as f:
        pol = json.load(f)
    ok = True
    files = [p for p in tracked(root) if any(p.startswith(d) for d in pol['governed_dirs'])]
    total = 0
    print(f"policy: {len(files)} governed files")
    for p in files:
        size = os.path.getsize(os.path.join(root, p))
        total += size
        cap = pol['grandfathered'].get(p)
        if cap is not None:
            if size > cap:
                print(f"  FAIL {p}: {size} bytes exceeds its grandfathered ceiling {cap}")
                ok = False
        elif size > pol['max_file_bytes']:
            print(f"  FAIL {p}: {size} bytes exceeds the {pol['max_file_bytes']} byte limit (docs/RESULTS_POLICY.md, size limits)")
            ok = False
    for p in pol['grandfathered']:
        if p not in files:
            print(f"  note {p}: grandfathered entry no longer tracked; remove it from tools/results_policy.json")
    if total > pol['max_total_bytes']:
        print(f"  FAIL governed files total {total} bytes, budget {pol['max_total_bytes']}")
        ok = False
    print(f"  total {total} bytes of {pol['max_total_bytes']}; limit {pol['max_file_bytes']} per file; {len(pol['grandfathered'])} grandfathered")
    with open(os.path.join(root, pol['provenance_file'])) as f:
        prov = f.read()
    names = set()
    for m in re.finditer(r'`([^`]+)`', prov):
        for e in expand_braces(m.group(1)):
            names.add(e.strip())
    tracked_set = set(tracked(root))
    for p in pol['generated_outputs']:
        if p not in tracked_set:
            print(f"  FAIL generated output {p} is not tracked")
            ok = False
        elif not any(n == p or n == os.path.basename(p) or fnmatch.fnmatch(p, n) for n in names):
            print(f"  FAIL {p} has no row in {pol['provenance_file']}")
            ok = False
    print(f"  provenance: {len(pol['generated_outputs'])} generated outputs checked against {pol['provenance_file']}")
    return ok

# ---------------------------------------------------------------- regeneration

def run(cmd, cwd, env=None, out=None, ok_codes=(0,)):
    e = dict(os.environ)
    e.update(env or {})
    print('  $', ' '.join(cmd), f'(cwd {cwd})', flush=True)
    t0 = time.time()
    r = subprocess.run(cmd, cwd=cwd, env=e, stdout=out if out else subprocess.DEVNULL, stderr=subprocess.STDOUT if out else None)
    print(f'    ({time.time() - t0:.1f} s)', flush=True)
    if r.returncode not in ok_codes:
        print(f'  FAIL: exit status {r.returncode}')
        return False
    return True


def committed_gds_epoch(path):
    """BGNLIB modification date of a committed GDS file as a Unix time, so a regenerated file differs in no field."""
    import calendar
    with open(path, 'rb') as f:
        d = f.read(64)
    i = 0
    while i + 4 <= len(d):
        n, rec = struct.unpack('>HB', d[i:i + 3])
        if rec == 0x01:
            y, mo, dd, h, mi, s = struct.unpack('>6h', d[i + 4:i + 16])
            return calendar.timegm((y, mo, dd, h, mi, s))
        i += n
    raise ValueError('no BGNLIB record in ' + path)


def regen(root, build, tier):
    ok = True
    build = os.path.abspath(build)
    with open(POLICY) as f:
        pol = json.load(f)
    tmp = tempfile.mkdtemp(prefix='qpuforge-repro-')
    print(f'regenerating into {tmp} (tier {tier}); results are compared to the committed files, which are not touched')
    try:
        # 1. design_study: 300 trajectories, ~30 s. It writes results/design_study.csv relative to its working directory.
        work = os.path.join(tmp, 'ds')
        os.makedirs(os.path.join(work, 'results'))
        with open(os.path.join(work, 'design_study.txt'), 'wb') as out:
            ok &= run([os.path.join(build, 'design_study')], work, out=out)
        ok &= compare_files('results/design_study.csv', os.path.join(root, 'results/design_study.csv'),
                            os.path.join(work, 'results/design_study.csv'), 'exact')
        ok &= compare_files('results/design_study.txt', os.path.join(root, 'results/design_study.txt'),
                            os.path.join(work, 'design_study.txt'), 'exact')
        # 2. GDSII from the committed design + cell files (documented command in results/PROVENANCE.md). Only AF-17 is cheap enough
        #    for this tier; BF49.gds and GF97.gds are rewritten by the review programs and checked in the full tier.
        for g in pol['gds']:
            out_gds = os.path.join(tmp, os.path.basename(g['out']))
            epoch = str(committed_gds_epoch(os.path.join(root, g['out'])))
            cmd = [os.path.join(build, 'qpuforge'), 'gds', g['chip'], out_gds, '--cells', g['cells'], '--cell', g['cell']]
            # `qpuforge gds` exits 2 when the layout fails DRC (all three first-cut layouts do: known issues in the CHANGELOG)
            # and 0 when it is clean; the file is written either way.
            ok &= run(cmd, root, env={'SOURCE_DATE_EPOCH': epoch}, ok_codes=(0, 2))
            ok &= compare_files(g['out'], os.path.join(root, g['out']), out_gds, 'exact')
        # 3. specs from review JSON (pure function of the committed review data and the template)
        for kind, d, spec in (('beta', 'designs/beta1', 'BETA-1.md'), ('gamma', 'designs/gamma1', 'GAMMA-1.md')):
            sd = os.path.join(tmp, 'spec-' + kind)
            os.makedirs(sd)
            shutil.copy(os.path.join(root, d, f'{kind}1_review.json'), sd)
            ok &= run([sys.executable, 'tools/fill_spec.py', kind, sd], root)
            ok &= compare_files(f'{d}/{spec}', os.path.join(root, d, spec), os.path.join(sd, spec), 'exact')
        # 4. full tier: the review programs at their committed shot counts
        if tier == 'full':
            for prog, d, json_out, txt_out, gds_out, extra in pol['reviews']:
                od = os.path.join(tmp, prog)
                os.makedirs(od)
                with open(os.path.join(od, 'stdout.txt'), 'wb') as out:
                    ok &= run([os.path.join(build, prog), '100000', '--out', od], root, out=out)
                ok &= compare_files(f'{d}/{json_out}', os.path.join(root, d, json_out), os.path.join(od, json_out))
                ok &= compare_files(f'{d}/{txt_out}', os.path.join(root, d, txt_out), os.path.join(od, 'stdout.txt'))
                for x in extra:
                    ok &= compare_files(f'{d}/{x}', os.path.join(root, d, x), os.path.join(od, x))
                if gds_out:
                    ok &= compare_files(f'{d}/{gds_out}', os.path.join(root, d, gds_out), os.path.join(od, gds_out))
    finally:
        if os.environ.get('KEEP_REPRO_DIR'):
            print('kept', tmp)
        else:
            shutil.rmtree(tmp, ignore_errors=True)
    return ok


def compare_trees(a, b):
    """Committed outputs of two checkouts, masked: evidence that two regenerations agree."""
    with open(POLICY) as f:
        pol = json.load(f)
    files = list(pol['compare_files'])
    ok = True
    print(f'comparing {len(files)} committed outputs: A={a}  B={b}')
    for p in files:
        fa, fb = os.path.join(a, p), os.path.join(b, p)
        if not (os.path.exists(fa) and os.path.exists(fb)):
            print(f'  SKIP {p}: missing on one side')
            continue
        ok &= compare_files(p, fa, fb)
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('policy')
    r = sub.add_parser('regen')
    r.add_argument('--tier', choices=['fast', 'full'], default='fast')
    r.add_argument('--build', required=True)
    c = sub.add_parser('compare')
    c.add_argument('a')
    c.add_argument('b')
    args = ap.parse_args()
    root = os.getcwd()
    if args.cmd == 'policy':
        ok = check_policy(root)
    elif args.cmd == 'regen':
        ok = regen(root, args.build, args.tier)
    else:
        ok = compare_trees(args.a, args.b)
    print('RESULT:', 'PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
