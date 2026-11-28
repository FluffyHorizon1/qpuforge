import shutil, sys, os
SRC = '/home/user/wt-fab/include'
def make(name, edits):
    d = f'mut_{name}/include'
    if os.path.exists(f'mut_{name}'): shutil.rmtree(f'mut_{name}')
    shutil.copytree(SRC, d)
    p = f'{d}/qpuforge/layout.hpp'
    s = open(p).read()
    for old, new in edits:
        n = s.count(old)
        assert n == 1, (name, n, old[:80])
        s = s.replace(old, new)
    open(p, 'w').write(s)
    print('made', name)

# M1: XY line drawn back into the qubit's (+x,-y) cut-out (the drawing on origin/main)
make('xy_back', [
 ("""    cpw_path({{cx + xy_end, cy - xy_level}, {cx + xy_end + 20, cy - xy_level}}, 10, 6, xm, xh, false);
    xm.push_back(rect_c(cx + xy_end + 32, cy - xy_level, 36, 36)); xh.push_back(rect_c(cx + xy_end + 32, cy - xy_level, 52, 52));""",
  """    (void)xy_end; (void)xy_level;
    stair(cx + 52, cy - 52, cx + 156, cy - 156, xm, xh, false);
    cpw_path({{cx + 156, cy - 156}, {cx + 250, cy - 156}}, 10, 6, xm, xh, false);
    xm.push_back(rect_c(cx + 262, cy - 156, 36, 36)); xh.push_back(rect_c(cx + 262, cy - 156, 52, 52));"""),
 ("L.bumps.push_back(rect_c(cx + kXYEnd + 32, cy - kXYLevel, 30, 30));", "L.bumps.push_back(rect_c(cx + 262, cy - 156, 30, 30));")])

# M2: each qubit SQUID drawn from the plane to the open end of the qubit's XY line instead of stub -> tongue
make('jj_on_xy', [
 ("        xmon(X, cx, cy, qm, qh, &L.jj);",
  "        { std::vector<Rect> jq; xmon(X, cx, cy, qm, qh, &jq); }\n        L.jj.push_back(Rect{cx + kXYEnd - 11, cy - kXYLevel - 2, cx + kXYEnd - 5, cy - kXYLevel + 2});")])

# M2b: the first qubit gets two SQUIDs, the second none (junction count unchanged)
make('jj_dup', [
 ("        xmon(X, cx, cy, qm, qh, &L.jj);",
  "        { std::vector<Rect> jq; xmon(X, cx, cy, qm, qh, &jq); static int qn = 0;\n          if (qn != 1) L.jj.insert(L.jj.end(), jq.begin(), jq.end()); if (qn == 0) L.jj.insert(L.jj.end(), jq.begin(), jq.end()); ++qn; }")])

# M3: no ground straps
make('no_straps', [("        L.straps = detail::place_ground_straps(map1, I, eligible, plane);", "        L.straps.clear(); (void)I;")])

# M4: qubit flux line drawn open again
make('flux_open', [("        qubit_control_lines(cx, cy, fm, fh, xm, xh, true);", "        qubit_control_lines(cx, cy, fm, fh, xm, xh, false);")])

# M5: SQUID ground lead (tongue) not cut out of the qubit's cut-outs
make('no_tongue', [("auto pcs = detail::subtract(holes[i], tongue);", "std::vector<Rect> pcs{holes[i]}; (void)tongue;")])

# M6: strap placer may span qubits
make('strap_over_q', [
 ("bool bad = false; I.each(sp.grow(-1e-6), [&](int i) { if (!I.kind[std::size_t(i)]) bad = true; });",
  "bool bad = false; I.each(sp.grow(-1e-6), [&](int) {});"),
 ("for (int i : crossed) { const int kd = I.kind[std::size_t(i)]; if (!kd) bad = true; else (kd == 1 ? s.lines_crossed : s.couplers_crossed)++; }",
  "for (int i : crossed) { const int kd = I.kind[std::size_t(i)]; if (!kd) {} else (kd == 1 ? s.lines_crossed : s.couplers_crossed)++; }")])

# M7: junction rule weakened back to the 0.4 rule (two different nets, plane not required)
make('jj_old_rule', [("if (a >= 0 && b >= 0 && a != b && (a == plane || b == plane)) C.junctions_landed++;",
                      "if (a >= 0 && b >= 0 && a != b) C.junctions_landed++;")])

# M9 (combined): XY line back in the cut-out (as M1) AND each qubit SQUID drawn from the XY line's bump pad (far edge,
# x = cx+280) across its trench to the ground at x = cx+288, instead of stub -> tongue
make('combo', [
 ("""    cpw_path({{cx + xy_end, cy - xy_level}, {cx + xy_end + 20, cy - xy_level}}, 10, 6, xm, xh, false);
    xm.push_back(rect_c(cx + xy_end + 32, cy - xy_level, 36, 36)); xh.push_back(rect_c(cx + xy_end + 32, cy - xy_level, 52, 52));""",
  """    (void)xy_end; (void)xy_level;
    stair(cx + 52, cy - 52, cx + 156, cy - 156, xm, xh, false);
    cpw_path({{cx + 156, cy - 156}, {cx + 250, cy - 156}}, 10, 6, xm, xh, false);
    xm.push_back(rect_c(cx + 262, cy - 156, 36, 36)); xh.push_back(rect_c(cx + 262, cy - 156, 52, 52));"""),
 ("L.bumps.push_back(rect_c(cx + kXYEnd + 32, cy - kXYLevel, 30, 30));", "L.bumps.push_back(rect_c(cx + 262, cy - 156, 30, 30));"),
 ("        xmon(X, cx, cy, qm, qh, &L.jj);",
  "        { std::vector<Rect> jq; xmon(X, cx, cy, qm, qh, &jq); }\n        L.jj.push_back(Rect{cx + 280, cy - 158, cx + 288, cy - 154});")])

# M10: strap feet allowed on island metal (foot_comp's island test removed)
make('foot_on_island', [("    I.each(f.grow(-1e-6), [&](int) { hit = true; });\n    if (hit) return -1;", "    (void)I; (void)hit;")])
