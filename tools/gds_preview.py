#!/usr/bin/env python3
"""Independent GDSII reader (pure struct parsing) + PNG preview.

Usage: gds_preview.py chip.gds out.png [--zoom x0,y0,x1,y1]
Parses every record, checks structure (HEADER/BGNLIB/UNITS/.../ENDLIB), counts elements per layer and renders
layers 1 (etch), 2 (junctions), 3 (bumps), 4/5 (ground-strap spans/feet), 10 (outline) and text labels with matplotlib.
"""
import struct, sys
from collections import Counter

def real8(b):
    v = int.from_bytes(b, 'big')
    if v == 0: return 0.0
    sign = -1 if v >> 63 else 1
    exp = ((v >> 56) & 0x7f) - 64
    mant = (v & ((1 << 56) - 1)) / float(1 << 56)
    return sign * mant * 16.0 ** exp

def parse(path):
    data = open(path, 'rb').read()
    i, recs = 0, []
    while i < len(data):
        n, rt, dt = struct.unpack('>HBB', data[i:i + 4])
        if n < 4: raise ValueError(f'bad record length at {i}')
        recs.append((rt, dt, data[i + 4:i + n])); i += n
    lib = {'structs': {}, 'units': None}
    cur = None; el = None
    for rt, dt, d in recs:
        if rt == 0x03: lib['units'] = (real8(d[:8]), real8(d[8:16]))
        elif rt == 0x06: cur = {'polys': [], 'texts': []}; lib['structs'][d.rstrip(b'\0').decode()] = cur
        elif rt in (0x08, 0x0C): el = {'kind': 'B' if rt == 0x08 else 'T', 'layer': None}
        elif rt == 0x0D: el['layer'] = struct.unpack('>h', d)[0]
        elif rt == 0x10:
            xy = struct.unpack('>%di' % (len(d) // 4), d)
            el['xy'] = [(xy[k] / 1000.0, xy[k + 1] / 1000.0) for k in range(0, len(xy), 2)]
        elif rt == 0x19: el['s'] = d.rstrip(b'\0').decode()
        elif rt == 0x11:
            (cur['polys'] if el['kind'] == 'B' else cur['texts']).append(el); el = None
    if recs[0][0] != 0x00 or recs[-1][0] != 0x04: raise ValueError('missing HEADER/ENDLIB')
    return lib

def main():
    path, out = sys.argv[1], sys.argv[2]
    zoom = None
    if '--zoom' in sys.argv: zoom = [float(v) for v in sys.argv[sys.argv.index('--zoom') + 1].split(',')]
    lib = parse(path)
    print('units (user/db, db in m):', lib['units'])
    for name, s in lib['structs'].items():
        c = Counter(p['layer'] for p in s['polys'])
        closed = all(p['xy'][0] == p['xy'][-1] for p in s['polys'])
        print(f'structure {name}: {len(s["polys"])} boundaries {dict(sorted(c.items()))}, {len(s["texts"])} labels, rings closed: {closed}')
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    s = next(iter(lib['structs'].values()))
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    ax.set_facecolor('#c9a227')   # ground metal (Nb/Al) everywhere not etched
    colors = {1: '#1b2a49', 2: '#e63946', 3: '#2a9d8f', 10: 'none'}
    for layer in (1, 2, 3):
        polys = [p['xy'][:-1] for p in s['polys'] if p['layer'] == layer]
        if polys: ax.add_collection(PolyCollection(polys, facecolors=colors[layer], edgecolors='none', linewidths=0))
    spans = [p['xy'][:-1] for p in s['polys'] if p['layer'] == 4]
    if spans: ax.add_collection(PolyCollection(spans, facecolors='#9b5de580', edgecolors='#9b5de5', linewidths=0.4))
    feet = [p['xy'][:-1] for p in s['polys'] if p['layer'] == 5]
    if feet: ax.add_collection(PolyCollection(feet, facecolors='#f15bb5', edgecolors='none', linewidths=0))
    for p in s['polys']:
        if p['layer'] == 10:
            xs, ys = zip(*p['xy']); ax.plot(xs, ys, color='k', lw=1)
    xs = [x for p in s['polys'] if p['layer'] == 10 for x, _ in p['xy']]; ys = [y for p in s['polys'] if p['layer'] == 10 for _, y in p['xy']]
    if zoom: ax.set_xlim(zoom[0], zoom[2]); ax.set_ylim(zoom[1], zoom[3])
    else: ax.set_xlim(min(xs), max(xs)); ax.set_ylim(min(ys), max(ys))
    for t in s['texts']:
        x, y = t['xy'][0]
        if not zoom or (zoom[0] <= x <= zoom[2] and zoom[1] <= y <= zoom[3]):
            ax.text(x, y, t['s'], fontsize=6 if not zoom else 9, color='white')
    ax.set_aspect('equal'); ax.set_xlabel('um'); ax.set_ylabel('um')
    ax.set_title(f'{next(iter(lib["structs"]))}  (gold: metal, navy: etch, red: junctions, teal: bumps, violet: ground straps)', fontsize=9)
    fig.tight_layout(); fig.savefig(out); print('wrote', out)

if __name__ == '__main__':
    main()
