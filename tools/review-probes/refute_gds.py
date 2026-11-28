# Minimal GDSII reader: count boundaries per layer and list layer-4/5 boxes near a point of the committed AF-17 mask.
import struct, sys
path = sys.argv[1]; qx, qy = float(sys.argv[2]), float(sys.argv[3])
data = open(path, 'rb').read()
i = 0; layer = None; uu = 1e-3; polys = []
def real8(b):
    s = -1 if b[0] & 0x80 else 1; e = (b[0] & 0x7f) - 64
    m = int.from_bytes(b[1:8], 'big') / (1 << 56)
    return s * m * 16 ** e
cur = None
while i < len(data):
    n, rt, dt = struct.unpack('>HBB', data[i:i+4]); body = data[i+4:i+n]; i += n
    if n == 0: break
    if rt == 0x03: uu = real8(body[0:8]); dbu_m = real8(body[8:16])
    elif rt == 0x08: cur = {'layer': None, 'xy': None}
    elif rt == 0x0D and cur is not None: cur['layer'] = struct.unpack('>h', body[:2])[0]
    elif rt == 0x10 and cur is not None:
        v = struct.unpack('>%di' % (len(body)//4), body); cur['xy'] = list(zip(v[0::2], v[1::2]))
    elif rt == 0x11 and cur is not None: polys.append(cur); cur = None
scale = dbu_m / 1e-6   # database unit in um
from collections import Counter
print('polygons per layer:', dict(sorted(Counter(p['layer'] for p in polys).items())))
for p in polys:
    if p['layer'] in (4, 5):
        xs = [x*scale for x, y in p['xy']]; ys = [y*scale for x, y in p['xy']]
        if min(xs) > qx - 200 and max(xs) < qx + 50 and min(ys) > qy - 250 and max(ys) < qy + 50:
            print(f"layer {p['layer']}: ({min(xs)-qx:.1f},{min(ys)-qy:.1f})-({max(xs)-qx:.1f},{max(ys)-qy:.1f}) rel. qubit")
