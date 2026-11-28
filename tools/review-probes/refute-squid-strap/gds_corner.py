# List layer 2 (junction), 4 (strap span), 5 (strap feet) boxes of a GDSII file inside a window around a qubit's
# (-x,-y) corner, relative to the qubit centre. Usage: gds_corner.py file.gds cx cy
import struct, sys

def real8(b):
    s = -1 if b[0] & 0x80 else 1
    e = (b[0] & 0x7f) - 64
    m = int.from_bytes(b[1:8], 'big') / float(1 << 56)
    return s * m * 16.0 ** e

path, cx, cy = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
data = open(path, 'rb').read()
i, cur, polys, dbu_um = 0, None, [], 1e-3
while i + 4 <= len(data):
    n, rt, _ = struct.unpack('>HBB', data[i:i + 4])
    if n < 4:
        break
    body = data[i + 4:i + n]
    i += n
    if rt == 0x03:
        dbu_um = real8(body[8:16]) / 1e-6
    elif rt == 0x08:
        cur = {'layer': None, 'xy': []}
    elif rt == 0x0D and cur is not None:
        cur['layer'] = struct.unpack('>h', body[:2])[0]
    elif rt == 0x10 and cur is not None:
        v = struct.unpack('>%di' % (len(body) // 4), body)
        cur['xy'] = [(v[k] * dbu_um, v[k + 1] * dbu_um) for k in range(0, len(v), 2)]
    elif rt == 0x11 and cur is not None:
        polys.append(cur)
        cur = None
counts = {}
for p in polys:
    counts[p['layer']] = counts.get(p['layer'], 0) + 1
print('boundaries per layer:', dict(sorted(counts.items())))
for p in polys:
    if p['layer'] not in (2, 4, 5):
        continue
    xs = [x - cx for x, _ in p['xy']]
    ys = [y - cy for _, y in p['xy']]
    if min(xs) > -200 and max(xs) < 0 and min(ys) > -200 and max(ys) < 0:
        print(f"layer {p['layer']}: ({min(xs):.1f},{min(ys):.1f})-({max(xs):.1f},{max(ys):.1f}) rel. qubit")
