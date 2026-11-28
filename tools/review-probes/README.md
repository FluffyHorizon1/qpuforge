# Review probes and findings for PR #28 (mask connectivity)

An independent adversarial review of this PR at `56ba3e4` (code identical to `ceb2b84`) ran three lenses: electrical and
physical correctness, claims against code and outputs, and tests and documentation. The reviewers wrote the probe
programs in this directory against the PR's headers and ran them on the assembled masks. This file lists every finding
with its evidence. The review was stopped before every finding had been checked by a second reviewer: the first
(shared SQUID/flux pocket on one airbridge) was upheld by two independent verifiers; the others rest on the probe
output quoted here.

**These files are working material for the redesign, not part of the product.** Promote what is useful into library
checks and tests, and delete this directory before the PR is merged.

Build a probe against the tree, for example:

```
g++ -std=c++17 -O2 -fopenmp -ffp-contract=off -Iinclude tools/review-probes/probe1.cpp -o /tmp/probe1
```

Some probes take the chip JSON and cell file as arguments or expect `-DPROBE_DIR=...`; read the top of each file.
Paths inside the `.log` files and in the evidence text below refer to the reviewers' scratch directory, which is this
directory's content.

## Findings

### [HIGH] The SQUID ground lead and the flux-line short share one ~100 um modelled airbridge, which couples the qubit galvanically to the 50 ohm flux line; line_coupling misses this, so the 22.4x / 9.7x flux T1 limits do not hold

*Lens: claims. File: `include/qpuforge/layout.hpp`, line 1106.*

The PR puts the SQUID's ground lead (tongue, xmon l.324) and the flux-line short (qubit_control_lines l.1071) in the same sealed (-x,-y) pocket. The only connection from that pocket to the rest of the ground is a single ground strap. The strap therefore carries the qubit's junction (mode) current for all of C_sigma except the pocket's own share (C_ret = C_sigma - C_qp), and it is also the flux line's termination to its CPW ground. A shared inductance L_s between the junction current and a 50 ohm port puts a series resistance (wL_s)^2/Z0 in the junction loop: Gamma = w^4 L_s^2 C_ret^2 / (Z0 C_sigma). line_coupling extracts with one ideal ground conductor (l.1106-1107; the pocket is part of G) and only weights capacitive panel charges by sin(beta z) from an ideal short. So it leaves out the dominant channel. The reported flux-line limits (22.4x AF-17, 9.7x Beta-1) and 'SQUIDs land island-to-plane' describe an ideal ground the drawn mask does not have. The design note says the line is 'shorted at the SQUID' and the pocket 'tied to the plane' (l.39-40, 57); its not-modelled list names only strap capacitance and the SQUID-loop mutual. Two related points. First, the flux current runs short -> pocket -> strap foot below the short. It never passes the SQUID, which sits at the tip of a dead-end tongue about 40 um from the short, so the (unmodelled) flux-bias coupling is also geometrically weak. Second, the straps form a Kruskal spanning forest, so each enclosed region hangs on exactly one bridge. 'Shorted into the plane' is therefore a DC statement only.

**Evidence.** Scratch probe probe_straps.cpp (output probe_straps_af17.log): all 17 qubits print 'pocket comp N (tongue same, ...): 1 strap(s) on it', e.g. 'Q16 ... strap rel. qubit (-81,-164)-(-71,-46), free span 98 um, crosses 1 line(s) 0 coupler/bus; other end comp 127 (flux-line CPW ground)'; summary 'pockets: 17; straps in series to the plane: 1:9 4:2 6:1 7:1 8:1 9:1 10:1 12:1; worst summed span 516 um'. On BF-49, 49 of 49 pockets have exactly one strap (free span 88-90 um), and up to 22 straps sit in series to the border plane ('worst summed span 850 um'). A split extraction of the line_coupling cell with the pocket as its own conductor (probe_lines.cpp, AF-17 measure cell) gives 'C_qubit-pocket 4.899 fF, C_sigma 95.90 fF', so C_ret = 91.0 fF. With Gamma = w^4 L_s^2 C_ret^2/(Z0 C_sigma) at Q16 (f01 5.99 GHz, design T1 65.5 us): L_s = 1 pH gives 289 us (4.4x, already below the 10x budget); 10 pH gives 2.9 us; 70 pH gives 0.059 us. Beta-1 Q10 at 1 pH: 285 us = 2.1x its 136.4 us. A 10 um x 100-118 um bridge has mu0 l/2pi [ln(2l/w)+1/2] = 70-88 pH. The same series L_s also detunes f01 by about -0.5 % (L_J is about 7.4 nH), which is not modelled either.

**Suggested fix.** Do not let a strap be in the SQUID's current path or the flux line's return. Give the (-x,-y) corner a wide galvanic connection to the flux line's CPW ground, for example by opening the T-bar cut-out overlap at that corner (T/2 + gap < L - t_w - gap there, or notch one cut-out), and land the SQUID lead on that continuous ground. Add a check that fails when a junction's ground side sits on a strap-connected component that also terminates a line. Model the strap inductance, or at least report the galvanic term, in line_budget. List strap inductance and RF grounding (single-bridge tree) in the design note's not-modelled section, and stop calling the flux line 'shorted at the SQUID'.

### [HIGH] The moved XY bump sits beside two coupler pads with the ground between them etched away; line_coupling draws no pads, so the XY budget (23.8x) misses about 1.5 fF of XY-to-coupler coupling

*Lens: claims. File: `include/qpuforge/layout.hpp`, line 1091.*

kXYEnd/kXYLevel put the XY bump 'in the free corner between the two coupler pads' (design note l.55). The bump's 52 um cut-out ends at cx+305 / cy-305. The +x and -y coupler-pad cut-outs start at p/2 - pad/2 - gap = 306.97 um. The 1.97 um ground slivers between them are removed by fill_ground_slivers (l.1362), so the bump faces both coupler pads across 30 um of bare etch. line_coupling draws the couplers with far_end = false (l.1091), and coupler_and_bus draws the pad only 'if (far_end)'. So the budget contains no XY-to-coupler capacitance. On the real mask the XY line couples about 45x more strongly to each coupler than to its qubit. Through the coupler (V_c/V_q = (C_qc/C_c) w^2/(w^2 - w_c^2)) this changes the qubit's effective XY coupling. It also gives the couplers an XY decay channel that no budget covers.

**Evidence.** Scratch probe probe_xy_window.cpp: a window of the assembled AF-17 mask around Q2 using the chip's own cut-outs after the fill, fine mesh 1.5/50/1.6, 7931 panels. Output: 'XY bump pad metal (2061.0,1503.0)-(2097.0,1539.0); ground at ... outside its +x and -y faces: [+x +4 um: etch ...] ... [+x +29 um: etch, -y +29 um: etch]'; 'C(XY2, C0-2) = 1.6368 fF', 'C(XY2, C6-2) = 1.5135 fF', 'C(XY2, Q2) = 0.0360 fF'. Inputs: C_qc 5.61 fF, C_c 155.25 fF (cell file), couplers parked at 8.398 GHz, so V_c/V_q = -0.0374 at 5.99 GHz. For Q5 (5.99 GHz, both pads present, design T1 65.2 us): C_eff = 0.036 - 0.0374 x 3.15 = -0.082 fF, giving T1_XY = 0.20 ms = 3.1x the design T1, against the reported >= 23.8x and the 10x budget. If the two paths add instead of cancel, it is 0.88x. The coupler's T1 limit from the XY line is C_c/(w_c^2 C_xc^2 Z0) = 0.42 us.

**Suggested fix.** Move the XY end and bump away from the coupler pads, or keep a ground strip of at least 4 um between their cut-outs (exclude it from the sliver fill). Extract line_coupling with the pads of the couplers that neighbour the XY end (far_end geometry, or a chip-window extraction). Include the coupler-mediated term in the qubit budget, and add a coupler T1 budget for the XY lines.

### [HIGH] Qubit SQUID ground and flux-line short share one enclosed pocket that reaches the plane only through a ~100 um airbridge: galvanic qubit-to-flux-line T1 limit of about 0.1-1 us

*Lens: electrical. File: `include/qpuforge/layout.hpp`, line 1071.*

The new tongue in xmon (l.320-328) ties the SQUID's ground electrode to the (-x,-y) inner-corner pocket. qubit_control_lines (l.1071) shorts the flux line into that same pocket. The pocket is sealed by the arm and T-bar cut-outs, so it reaches the plane only through one ground strap (a modelled airbridge).

About 95% of the island's charge ends outside the pocket. So the qubit's plasma current runs island -> SQUID -> pocket -> strap -> plane, and the 50-ohm flux line sits directly across the strap inductance L_s. That puts Re Z = (wL_s)^2/Z0 in series with L_J, giving Gamma1 ~ (wL_s)^2/(Z0 L_J).

Nothing in the PR sees this:
- line_coupling puts the pocket inside the single ideal 'ground' conductor (outline minus holes, l.1106-1108). Its comment calls the short 'at ground potential in the real line' (l.1097-1098).
- The sin(beta z) weighting assumes an ideal short.
- ConnectivityReport merges the strap as a perfect net (l.1413), and the junction check accepts it (l.1442).

So the claims 'SQUIDs land island-to-plane', 'flux lines shorted into the plane' and 'flux T1 limit 22.4x / 9.7x the design T1' are not supported electrically.

**Evidence.** Probe sources are in /tmp/claude-0/-home-user/a8a4cb10-31eb-5817-b046-9e2692d52262/scratchpad/fab-review/.

probe1 on the assembled masks (from the committed cells):
- AF-17: 'SQUID ground on the plane itself (no strap): 0 / 17; flux short in the SQUID-ground fragment: 17 / 17'.
- BF-49: 0/49 and 49/49. GF-97: 0/97 and 97/97.
- 'Q0 pocket strap 106: foot (-51.0,-76.5) -> foot (-159.0,-76.5) rel. qubit; bridge 98 um'. Q1 and Q2: 100 um.
- Per-qubit chains reach 12 straps / 516 um of bridge (Q15).

probe3, using the PR's own line_coupling cell with the pocket split out as its own conductor:
- The tongue is metal component 4; the cell plane is component 0; pocket area 757 um^2.
- 'C_sigma 93.94 fF; qubit-pocket 4.74 fF (5.0% of C_sigma)'.

Circuit island -L_J- pocket -(L_s || 50 ohm)- plane, solved for its complex resonance, with L_J from analyze() (E_J 17.72 / 23.73 GHz, i.e. 9.22 / 6.89 nH):
- L_s 2 pH: T1 111 / 62 us
- L_s 5 pH: T1 17.7 / 9.9 us
- L_s 20 pH: T1 1.11 / 0.62 us
- L_s 50 pH: T1 0.18 / 0.10 us

Grover partial self-inductance of a 98 um x 10 um strip is 68 pH. Design T1 is 65-76 us (AF-17) and 135-157 us (BF-49/GF-97). Meeting the PR's 10x budget would need L_s of about 0.6-0.8 pH.

**Suggested fix.** Land the SQUID ground lead on continuous plane metal: open the corner pocket to the outer plane with a wide metal path (do not let the T-bar cut-outs overlap), or run the lead out of the pocket. Short the flux line on continuous plane, so that no airbridge carries both the flux return and the SQUID's plasma current. Gate that each qubit SQUID's ground end lies on the plane component before straps are counted. Add the strap inductance to the control-line model.

### [HIGH] Flux-line short does not couple flux into the SQUID loop (M about 0.004-0.011 pH): the shorted flux lines cannot practically tune the qubits

*Lens: electrical. File: `include/qpuforge/layout.hpp`, line 1071.*

The short is the last stair square at (-61,-61) relative to the qubit: square [-66,-56]^2, the exact spot where the line ended open on main. The SQUID region is x[-34,-24], y[-23,-19] (l.328), 40-50 um away.

The return current runs from the short to the pocket-strap foot and over the bridge, away from the loop. The tongue is a dead end, so no current passes beside the loop. The loop region also shrank from 30x12 um to 10x4 um.

The design note says M is 'not modelled'. But the PR counts the flux lines as fixed: 'qubit flux lines 17/17' in the gate and test, and 'each qubit flux line is shorted to the plane' in the CHANGELOG. The device model treats the qubits as flux-tunable (df_dflux, flux-noise T_phi), and ROADMAP relies on flux-tunable qubits to fix frequency collisions.

**Evidence.** Biot-Savart filament model of the closed current path, integrated over the SQUID region (mloop.cpp and mseg.cpp in the scratch dir):
- Strap running left (Q0-Q5 type): M = -0.0036 pH, so I(Phi0/2) = 290 mA.
- Strap running down (Q14 type, strap 119): M = 0.0113 pH, so I(Phi0/2) = 91 mA.
- Every individual segment is at most 0.051 pH, so at least 20 mA would be needed even for the most favourable segment.
- Reference: current running up the tongue past the loop, (-39,-54) -> (-39,-18), gives M = -0.567 pH (about 1.8 mA for Phi0/2).

**Suggested fix.** Route the flux-line current past the loop: short the centre conductor onto the SQUID ground lead, so the current flows along the loop's side and returns on continuous plane. Compute M (filament or sheet-current model) and gate I(Phi0/2), e.g. at most 1-2 mA. Until then, do not report the qubit flux lines as fixed.

### [HIGH] New XY pad sits 30 um from 1-2 coupler pads in one merged cut-out, but line_coupling leaves the pads out: the '0.03 fF / 23.8x' XY result is wrong, the drive is coupler-dominated, neighbour crosstalk is 60-160% and Q6's XY T1 limit is 5.4x

*Lens: electrical. File: `include/qpuforge/layout.hpp`, line 1091.*

qubit_control_lines puts the XY end at (+247,-279) and its pad at (+279,-279), 'in the free corner between the two coupler pads' (design note l.55). In the mask, the 2 um ground strip between the XY pad cut-out and each coupler-pad cut-out is etched by the sliver fill, leaving 30 um of bare etch to the coupler pads.

line_coupling builds its cell with truncated couplers (coupler_and_bus(..., false): no pads), so in the cell the nearest coupler metal is about 200 um away. Below the coupler frequency, the coupler-mediated path C_qc*C_xc/(C_c(w_c^2/w^2 - 1)) has the opposite sign to the direct coupling and is about as large. It cancels most of the drive and drives the neighbour qubit as strongly as the target. It is also an unmodelled coupler T1 channel.

On its own, the 0.030 fF direct coupling would be adequate: a 25 MHz Rabi rate needs about -74 dBm on chip. So the move did not starve the drive; the defect is the unmodelled coupler path.

**Evidence.** xyhole: for Q1, along y = -280 rel, the metal map reads '295:15 300:-1 ... 325:-1 330:16' (XY pad, then etch, then coupler pad), and at(+306) = -1. Q6 is the same along both x and y.

probe5, using the PR's extract_capacitance on windows of the assembled AF-17 mask (h_min 3 and 2 um agree within 3%), coupler parked at 8.398 GHz:
- Q1: C(XY1,Q1) = 0.038 fF, C(XY1,C5-1) = 1.40-1.43 fF, C(Q1,C5-1) = 5.29-5.43 fF, C_coupler = 164-166 fF.
- At 5.11 GHz: Q1 net +0.010 to +0.012 fF; neighbour Q5 gets -0.015 to -0.016 fF (133-161% of Q1's).
- At 5.5 GHz the target's net coupling is only +0.0026 to +0.0042 fF.
- Q6 (C6-3 at 1.60 fF, C6-9 at 1.48 fF): net -0.0617 fF, so its own-XY T1 limit is 355 us = 5.4x its design T1 of 66 us. The PR reports at least 23.8x, 'budget met'.
- Neighbours of Q6: Q3 -0.043 fF, Q9 -0.038 fF.
- Coupler T1 limit through a single XY line: 0.46-0.60 us.

**Suggested fix.** Extract the control lines with full couplers (pads included), or from a window of the assembled mask. Report the effective drive including coupler-mediated paths, the crosstalk to neighbours and the coupler T1. Move the XY end away from the coupler pads, or keep a ground shield between the cut-outs wide enough that the sliver fill cannot remove it.

### [HIGH] Docs say the control-line T1 budget is gated and passed by all three chips; the code does not gate it and BF-49/GF-97 miss it

*Lens: tests-docs. File: `CHANGELOG.md`, line 115.*

Commit 5369c1e removed the `lines` criterion from all three reviews (the flux line misses the 10x budget on the Beta-1 cells). The docs commit 56ba3e4 came after it and still describes the budget as a passed gate:
- CHANGELOG.md:115-116: "The mask width/space DRC, metal connectivity and the control-line T1 budget are gated in the Alpha, Beta and Gamma reviews ... Alpha-1 v1.1: 8 of 8 criteria". At HEAD the Alpha review has 7 criteria.
- CHANGELOG.md:38 lists gate criteria "`connectivity` and `lines`"; no review has `lines`.
- CHANGELOG.md:156: "The masks pass the in-house checks (... metal connectivity, control-line budget)".
- CHANGELOG.md:114 and 117 still contain the placeholders `BF49_GF97_CONN` and `BETA_GAMMA_GATES`.
- README.md:255-257: "All three chips pass the in-house checks: ... and the control-line T1 budget".
- docs/WHITEPAPER.md:15: "... and the control-line T1 budget. BF-49 rev C and GF-97 pass the same checks." WHITEPAPER.md:169 says "Control lines (gated)".
- ROADMAP.md:18-19 say Beta and Gamma "criteria met in simulation, including ... control lines". ROADMAP.md:66/130/169 add "Control lines" as a simulation-gate row, but the reviews that print "Alpha simulation gate (ROADMAP.md)" no longer check it.
- designs/alpha1-v11/ALPHA-1-v1.1.md:186 and :201 present control lines as a gate row. This spec is hand-written (fill_spec.py does not generate it), so the planned regeneration of outputs will not fix it.
- The committed alpha1_v11_review.txt (6cd35a0) still shows a lines [PASS] row that HEAD no longer emits.
- docs/design/mask-connectivity.md:60-61 says line_budget "asks for at least 10x" but never says that the budget is not gated or that Beta-1 misses it.
fill_spec.py and the BETA/GAMMA templates do handle this correctly ("reported, not gated", "**not met**").

**Evidence.** Code: `git show 5369c1e -- examples/` shows `crit("lines", ...)` deleted in alpha1_v11_review.cpp, beta1_review.cpp and gamma1_review.cpp. beta1_review.cpp now reads: "The control-line T1 budget (layout.control_lines) is reported, not gated: the flux line misses it on the Beta-1 cells."
Independent fine-mesh run (MeshOptions{1.5,50,1.6}) of line_budget with probe /tmp/claude-0/-home-user/a8a4cb10-31eb-5817-b046-9e2692d52262/scratchpad/fab-review/tests-docs/lc.cpp, log line_budget_fine.log:
- BF-49: "T1 limits >= 1.56 ms (XY), 1.32 ms (flux): 11.4x (Q10) and 9.7x (Q10) the design T1 (budget 10x: NOT met)"
- GF-97: "... 9.7x (Q43) ... NOT met"
- Qubits below 10x: 17 of 49 (BF-49) and 36 of 97 (GF-97).
- AF-17: 23.8x / 22.4x, met.

**Suggested fix.** Say everywhere that the control-line budget is reported, not gated. AF-17 meets it (23.8x XY / 22.4x flux, as coded). The Beta-1 cells (BF-49, GF-97) miss it on the flux line (9.7x; 17/49 and 36/97 qubits below 10x). Make the ROADMAP rows "reported" or remove them. Drop `lines` from the Added list and correct "8 of 8". Fill or delete the BF49_GF97_CONN and BETA_GAMMA_GATES placeholders. Fix README:257, WHITEPAPER:15/169 and ALPHA-1:186/201 by hand. Add the outcome to the design note.

### [MEDIUM] Docs say the control-line T1 budget is gated and passed on all three chips; the code does not gate it, and Beta-1 / GF-97 fail it (9.7x)

*Lens: claims. File: `ROADMAP.md`, line 18.*

Commit 5369c1e removed the 'lines' criterion from all three reviews because the flux line reaches only 9.7x on the Beta-1 cells. The docs commit after it (56ba3e4) still describes the budget as gated and met everywhere. ROADMAP lists 'Control lines' as a gate row under 'all must hold' (l.66, 130, 169) and reports Beta and Gamma criteria as met 'including ... control lines' (l.18-19). README says all three chips pass 'the control-line T1 budget'. CHANGELOG claims a `lines` gate criterion and '8 of 8 criteria'. WHITEPAPER marks control lines '(gated)' and says BF-49/GF-97 'pass the same checks'. The AF-17 outputs regenerated in 6cd35a0 (ALPHA-1-v1.1.md l.7, 186, 201; alpha1_v11_review.json exit_criteria.lines = true; review.txt '[PASS] each qubit's XY and flux line T1 limits ...') show a criterion the head code no longer computes; it now produces 7 criteria, not 8. Only the code, its output text and the Beta/Gamma templates state 'reported, not gated'.

**Evidence.** examples/alpha1_v11_review.cpp, beta1_review.cpp, gamma1_review.cpp at 56ba3e4 have no crit("lines", ...): 5369c1e deleted it, and lines_ok in review_util.hpp:104 is unused. README.md:255-257: 'All three chips pass the in-house checks: ... and the control-line T1 budget.' CHANGELOG.md:38: 'Gate criteria `connectivity` and `lines` in the three reviews'. CHANGELOG.md:115-116: 'the control-line T1 budget are gated ... Alpha-1 v1.1: 8 of 8 criteria'. CHANGELOG.md:155-156: 'The masks pass the in-house checks (... control-line budget)'. WHITEPAPER.md:15: '... and the control-line T1 budget. BF-49 rev C and GF-97 pass the same checks.' WHITEPAPER.md:169: 'Control lines (gated)'. probe_lines BF-49 prints '11.4x (Q10) and 9.7x (Q10) the design T1 (budget 10x: NOT met)'.

**Suggested fix.** State in README, ROADMAP (move the 'Control lines' rows out of the 'all must hold' tables or mark them reported-only), CHANGELOG and WHITEPAPER that the control-line budget is reported, not gated, and that BF-49/GF-97 miss it. Fix '8 of 8' and the `lines` criterion entry. Regenerate the AF-17 spec and review so exit_criteria match the head code.

### [MEDIUM] Connectivity gate treats a minimum spanning tree of single airbridges as ideal ground (chains up to 34 straps; one 36 um bridge carries the returns of 77 circuits)

*Lens: electrical. File: `include/qpuforge/layout.hpp`, line 1058.*

place_ground_straps keeps exactly one spanning tree: 122 straps for 123 ground nets in AF-17, 373 for 374 in BF-49, 774 for 775 in GF-97. Every fragment hangs on a single bridge with no redundancy. ConnectivityReport merges each landed strap into an ideal net (l.1413), and clean() only checks net membership.

The result passes as 'clean' (layout_presets: CHECK(L.conn.clean())) even though, at microwave and flux-pulse frequencies:
- SQUID grounds reach the plane through long chains of bridges;
- a single bridge is shared by the return currents of many SQUIDs, flux lines and resonators;
- coupler flux lines and readout resonators are shorted into the same strapped fragment.
Tens of pH of shared bridge inductance couple all of these circuits galvanically. The 'clean' result says nothing about ground quality.

**Evidence.** probe1:
- Strap tree depth max 12 (AF-17), 22 (BF-49), 34 (GF-97).
- Qubit-SQUID-to-plane bridge length: max 516 / 850 / 1284 um.
- Most-shared strap is a 36 um bridge carrying the return of 28 / 49 / 77 circuits. AF-17 example: 'SQ6 SQ8 SQ11 SQ15 SC0-2 ... FC11-8 ... FQ6 FQ8 R0 R2 R6 R9'.
- Fragments holding more than one shorted line: 8 / 32 / 72 (e.g. 'fragment 4 holds FC4-7 R1').
- Coupler SQUIDs on a strapped fragment that also holds their own flux short: 17 of 24, 66 of 80, 147 of 168.

**Suggested fix.** Gate ground quality, not only net membership. SQUID grounds and line shorts should sit on the plane itself, or within a bounded strap inductance. Require several spaced straps per fragment boundary, a maximum chain length, and no strap shared by the returns of different circuits. Describe the straps in the docs as a connectivity model only.

### [MEDIUM] Junction-landing gate never checks that a SQUID lands on its own qubit/coupler, or that each has exactly one; moving every qubit SQUID onto the XY line passes layout_presets

*Lens: tests-docs. File: `include/qpuforge/layout.hpp`, line 1442.*

The rule `a >= 0 && b >= 0 && a != b && (a == plane || b == plane)` accepts any non-plane net at the island end. That includes an XY line, a bus, another qubit, or a second junction on the same island. The test comment (tests/test_qpuforge.cpp:1920: "every SQUID from its island to the plane") and the CHANGELOG (line 40: layout_presets requires "one junction per qubit and coupler") claim more than is checked. The only count check (test line 1922) compares the total number of jj rectangles with qubits+couplings. No negative test pins the new island-to-plane rule either. The current mask is correct: a probe finds every junction on its own island and exactly one per Q/C island, 41/41 and 129/129. Only the gate is weaker than documented.

**Evidence.** Mutation runs on copies of include/ in /tmp/claude-0/-home-user/a8a4cb10-31eb-5817-b046-9e2692d52262/scratchpad/fab-review/tests-docs/:
(a) `combo` mutant: each qubit SQUID drawn from the XY line's bump pad to ground instead of stub to tongue, and the XY line drawn back in the cut-out. The real `test_qpuforge layout_presets` gives "21 checks passed, 0 failed", exit 0, with "junctions landing island-to-plane 41/41 ... 129/129 ... 265/265" (test_combo_layout_presets.log).
(b) `jj_dup` mutant: Q0 gets two SQUIDs and Q1 none. jj_owner.cpp reports "Q/C islands with exactly one junction: 39 of 41", yet the layout_presets mask/connectivity assertions (probe presets_mask.cpp, a copy of test lines 1915-1923) PASS.
(c) `jj_old_rule` mutant: the 0.4 rule `a != b` restored. PASS (mutants_af17.log).

**Suggested fix.** When drawing (xmon, coupler loop), store each junction's owner island index. Count a junction as landed only if one end's net is the plane and the other end's net is the owner's net. Require exactly one junction per Q*/C* island. Add a synthetic layout_geometry case where a junction between two islands, between an island and a floating fragment, and between an XY line and the plane must not count as landed.

### [MEDIUM] The XY-line move and the whole control-line model are untested: drawing the XY line back into the qubit cut-out (3.2 fF, T1 limit 0.15 us) passes every layout test

*Lens: tests-docs. File: `tests/test_qpuforge.cpp`, line 1921.*

The PR's control-line fix (XY 3.0 fF -> 0.03 fF) has no regression protection. No test calls line_coupling, line_budget, qubit_control_lines or the new extract_capacitance panel_charge output, and none checks the XY position. The extraction cells used by layout_presets (star_cell/pair_cell) contain no XY or flux line. The reviews only print the budget "(reported, not gated)". The connectivity check cannot see the XY line's position, because an open XY line is its own net wherever it is.

**Evidence.** `grep -rn 'line_coupling|line_budget|kXYEnd|qubit_control_lines|PanelCharge' tests/` finds nothing.
Mutant `xy_back` restores the origin/main drawing (`stair(cx + 52, cy - 52, cx + 156, cy - 156, ...)`):
- line_coupling, coarse mesh: "C_xy 3.1963 fF ... T1_xy 0.153 us", against 0.0295 fF / 1.85e3 us for the PR geometry on the review mesh.
- `test_xy_back layout layout_geometry`: "34 checks passed, 0 failed".
- layout_presets mask/connectivity assertions: PASS.
- The combined mutant, which includes this change, passes the full layout_presets group 21/21.

**Suggested fix.** Add a check on the committed cells. A coarse mesh {3.0, 50, 2.0} takes about 2 s per cell, and the probe matches the review to 1 %. Assert that line_coupling gives C_xy and C_flux_eff below, say, 0.1 fF, and that the XY metal lies outside the qubit cut-out. Add a unit test of the panel_charge output (the panel charges of conductor j sum to the Maxwell entry) and of the T1 / sin(beta z) formula.

### [MEDIUM] line_budget reuses the flux coupling computed at 5.5 GHz at every qubit frequency, so flux-line T1 limits are overstated for f01 > 5.5 GHz

*Lens: tests-docs. File: `include/qpuforge/layout.hpp`, line 1149.*

line_coupling is called with f = 5.5 GHz (lines 1144-1145). Line 1117 builds beta from that frequency, and C_flux_fF = sum q_k sin(beta z_k) (line 1124). line_budget then applies 1/T1 = w^2 C^2 Z0 / C_sigma with w = 2 pi f01 but the 5.5 GHz C_flux (lines 1149-1150). On this line beta z <= ~0.08 rad, so sin(beta z) is beta z to 0.1 % and the effective capacitance scales with f. The shorted-line limit therefore goes as f^-4, while the code gives f^-2. For each qubit above 5.5 GHz the flux T1 limit is overstated by (f01/5.5)^2. The XY (open line) limit is correct. The design note's "evaluates the T1 limit of each line at every qubit's idle frequency" and the doc figures (22.4x, 1.47 ms, Beta-1 9.7x) inherit the error.

**Evidence.** Probe lc.cpp, fine mesh, log line_budget_fine.log, comparing each qubit's ratio as coded with C_flux scaled by f01/5.5:
- BF-49: "worst flux ratio 9.71x (Q10) as coded, 8.19x (Q10) with C_flux scaled to f01". Q10 is at f01 = 5.990 GHz. Qubits below 10x: 17 -> 24 of 49.
- GF-97: 9.71x -> 8.19x (Q43). Below 10x: 36 -> 48 of 97.
- AF-17: "22.42x (Q16) as coded, 18.90x (Q16)". Still above 10x.

**Suggested fix.** Evaluate the flux coupling at each qubit's f01. For example, keep the panel charges and z (or sum q z) in LineCoupling and recompute sin(beta(f01) z) inside line_budget, or call line_coupling once per distinct frequency. Then regenerate the documented ratios.

### [MEDIUM] Docs present the strapped pockets as tied to the plane and dropped the ideal-ground caveat; SQUID ground leads and flux returns actually run through chains of up to 12/22/34 modelled airbridges, shared by several qubits

*Lens: tests-docs. File: `docs/design/mask-connectivity.md`, line 39.*

The LVS-lite (layout.hpp:1408-1417) merges any chain of straps into the plane net. place_ground_straps builds a minimum-cost Kruskal forest, so the enclosed regions are chained to each other rather than each tied to the plane. Each qubit's (-x,-y) pocket carries both the SQUID ground lead and the flux-line short, and it hangs on one strap. Nothing checks or reports how many bridges are in series, or that return currents share bridges. The extraction cells still model these pockets as ideal ground: cell_ground_holes produces one ground conductor.

The PR nonetheless deleted the caveat "the extraction treats every enclosed ground fragment as grounded" from the CHANGELOG Known issues, the WHITEPAPER and ALPHA-1-v1.1.md (diff lines with '-' at CHANGELOG ~137 and the ALPHA-1/WHITEPAPER connectivity bullets). The design note says "the pocket is tied to the plane by ground straps" (lines 39-40). Its "Not modelled" list (lines 71-76) names only strap capacitance, the flux-SQUID mutual inductance and the foundry deck. Strap inductance and shared return paths are not mentioned anywhere, so the docs overstate what "connectivity clean" and "SQUIDs land island-to-plane" verify.

**Evidence.** Probe strap_stats.cpp on the HEAD masks (log strap_stats.log). In every qubit, the SQUID ground lead and the flux-line short are the same metal component (17/17, 49/49, 97/97), touched by exactly one strap. Bridges in series from that component to the plane:
- AF-17: 9 qubits through 1 bridge, the rest 4-12; up to 4 qubits' flux returns share one bridge.
- BF-49: up to 22 bridges; up to 7 share one bridge.
- GF-97: up to 34 bridges; up to 11 share one bridge.
No return path passes through another qubit's SQUID pocket.

**Suggested fix.** Report the depth and sharing of the strap chains in ConnectivityReport/connectivity_summary (maximum bridges from each SQUID/flux short to the plane, maximum returns per strap). Restore a caveat: the pockets are DC-connected through modelled airbridges; their inductance, the shared returns and the extraction's ideal-ground assumption are not modelled. Consider tying each SQUID/flux pocket directly to the plane, with a bound on chain length and some redundancy.

### [LOW] Flux-line standing-wave weights are computed at 5.5 GHz but rescaled to each qubit's f01 by w^2 only, which overstates the flux T1 limit by (f01/5.5)^2

*Lens: claims. File: `include/qpuforge/layout.hpp`, line 1150.*

line_budget calls line_coupling(..., 5.5, ...) (l.1144-1145). That bakes sin(beta z) at 5.5 GHz into C_flux_fF. It then evaluates T1 = C_sigma/(w^2 C_flux^2 Z0) with w = 2 pi f01 (l.1150). For the flux line C_eff is about beta * sum(C_i z_i), which is proportional to f, so T1 scales as f^-4, not f^-2. The worst qubits are measure qubits at 5.81-5.99 GHz, so the reported minima are optimistic. The design note (l.60-61) says the limit is evaluated 'at every qubit's idle frequency'.

**Evidence.** probe_lines.cpp output. AF-17: 'qubit f01 range 5.110-5.990 GHz; flux ratio with the weights at each qubit's f01: 18.90x (Q16) vs reported 22.42x (Q16)'. Beta-1: '8.19x (Q10) vs reported 9.71x (Q10)'.

**Suggested fix.** Compute the sin(beta z) weights at each qubit's f01. Return the per-panel (z, q) list from line_coupling and weight it per qubit, or call line_coupling per distinct frequency.

### [LOW] The new control-line model (line_coupling, line_budget, panel charges) has no test, and its z-map hard-codes the flux-line geometry

*Lens: claims. File: `tests/test_qpuforge.cpp`, line 1921.*

The added checks cover the strap placer and the chips' connectivity and DRC. Nothing exercises line_coupling, line_budget or extract_capacitance's panel_charge output: no check of the 1/T1 = w^2 C^2 Z0/C_sigma form, the sin(beta z) weighting, or that the AF-17 budget holds. The z-coordinate map (l.1101, 1123) re-encodes the stair of qubit_control_lines as literals (95, -61, -156, -161, -140). A change to the line geometry, or a sign error in the weighting, would silently change the published T1 limits. The comment at l.1066 already disagrees with the code ('y = cy - 156'; the code uses cy - kXYLevel = cy - 279).

**Evidence.** grep -n 'line_coupling\|line_budget\|PanelCharge\|panel_charge' tests/*.cpp returns no match. The only callers are examples/alpha1_v11_review.cpp:184 and examples/review_util.hpp:103. CHANGELOG: 'Four checks more than before (one and three)', none of them on control lines.

**Suggested fix.** Add a unit test. For example: the panel charges summed over the flux conductor equal -C(qubit, flux) from the matrix; C_eff/C_open is about beta * <z> for a known geometry; an analytic T1 check; and a regression bound on the AF-17 cell values. Derive z from the drawn stair (pass it out of qubit_control_lines) instead of literals.

### [LOW] The junction-landing check accepts any non-plane conductor on the far end, not the junction's own island

*Lens: claims. File: `include/qpuforge/layout.hpp`, line 1442.*

The claim is 'SQUIDs island-to-plane'. The predicate is 'one end on the plane net, the other on any other net'. It never checks that the non-plane end is the qubit or coupler island the junction belongs to. A SQUID region drawn from the ground lead onto an XY line, a bus or a neighbouring island would count as landed, so conn.clean(), the layout_presets checks and the gated 'connectivity' criterion would pass a mis-landed SQUID. The current masks are correct: the qubit region spans tongue x = s0-10 to stub x = s0, and coupler regions span pad to ground.

**Evidence.** layout.hpp:1438-1442: 'const int a = net(horiz ? map1.at(j.x0 - d, j.cy()) : ...); const int b = net(...); if (a >= 0 && b >= 0 && a != b && (a == plane || b == plane)) C.junctions_landed++;'. L.jj holds only rectangles, with no owner.

**Suggested fix.** Store the owning island with each junction region (qubit id / coupler tag). Require the non-plane end's net to equal that island's net, and that net to hold exactly that island.

### [LOW] Straps have no clearance or foot-enclosure rules and no DRC on GDS layers 4/5: bridges flush (0.00-0.22 um) with qubit T-bar ends

*Lens: electrical. File: `include/qpuforge/layout.hpp`, line 1042.*

The only rule is that a strap must not strictly intersect qubit rectangles (l.1042), so a bridge may run flush along a qubit's T-bar end, where the qubit's field is strongest. Feet are tested at only 9 sample points, 0.25 um inside the edge (l.988), so they can sit flush with, or slightly over, the ground edge. The DRC runs on the etch layer only (l.1448); layers 4/5 are written (l.1499) but never checked. The capacitance and loss of bridges next to the qubits are not in the extraction.

**Evidence.** probe2 (AF-17):
- 'strap 119 ... 0.00 um from Q14 rect x[-60.0,60.0] y[-125.7,-113.7]; side-by-side run 12.0 um'.
- Straps 31, 32, 39, 41, 42, 45, 47 and 49 are 0.22 um from the -x T-bar end of Q1/Q4/Q6/Q10/Q13/Q11/Q16/Q12.
- 29 of 122 straps are within 10 um of qubit metal (BF-49: 147 of 373; GF-97: 291 of 774, five at 0.00 um).
- Feet: 107 of 244 flush with an etch edge; 44 overhang etch by up to 2.5 um^2 (0.25 x 10 um). BF-49: 380 and 124 of 746. GF-97: 786 and 270 of 1548.

**Suggested fix.** Add strap rules: lateral clearance from qubit islands and junctions; foot enclosure of at least 1-2 um, checked on the whole foot rectangle; strap-to-strap spacing. Run width/space on layers 4/5, or keep straps out of the qubit cut-outs.

### [LOW] line_budget applies the 5.5 GHz standing-wave weighting at every qubit's f01, overestimating the flux T1 by about 19% at 5.99 GHz

*Lens: electrical. File: `include/qpuforge/layout.hpp`, line 1144.*

C_flux_eff is the panel charge weighted by sin(beta z) with beta taken at 5.5 GHz (l.1144-1145). line_budget then evaluates T1 at each qubit's own f01 (l.1149-1151) with that same coupling. Since beta scales with f, the weighted coupling grows with frequency and T1 falls faster than the formula used.

**Evidence.** lcf (coarse mesh): 'C_flux_eff weighted at 5.50 GHz: 0.02993 fF; at 5.99 GHz: 0.03260 fF (ratio 1.089)'. 'Q16 ... as line_budget computes it 1480 us = 22.6x; with the 5.99 GHz weighting 1248 us = 19.1x'. Beta-1's 9.7x drops further for qubits above 5.5 GHz.

**Suggested fix.** Store the panel charges once and re-weight them at each qubit's frequency. The ideal-short assumption from the first finding also needs fixing.

### [LOW] ROADMAP and CHANGELOG say the control-line budget is gated and met on Beta and Gamma; the code reports it ungated because Beta-1's flux line misses 10x

*Lens: electrical. File: `ROADMAP.md`, line 18.*

The docs added in this PR say Beta and Gamma met their criteria 'including ... control lines' and that the control-line T1 budget is gated. The last code commit un-gated it because the Beta-1 flux line is at 9.7x. With the first and third findings, the budget is not met on any chip.

**Evidence.** - ROADMAP.md:18 and :19: 'criteria met in simulation, including mask DRC, metal connectivity and control lines'.
- CHANGELOG.md:115: 'The mask width/space DRC, metal connectivity and the control-line T1 budget are gated'.
- Commit 5369c1e: 'The flux line meets 9.7x the design T1 on the Beta-1 cells (budget 10x), so the criterion is not gated'.
- examples/beta1_review.cpp: 'The control-line T1 budget (layout.control_lines) is reported, not gated'.

**Suggested fix.** State that the control-line budget is reported, not gated, and is missed. Do not call control lines 'met' on any chip until the coupling paths and strap inductance are modelled.

### [LOW] The safeguard that keeps strap feet off island metal is untested; if it regresses, airbridge feet land on resonator and flux-line centre conductors while connectivity stays clean

*Lens: tests-docs. File: `include/qpuforge/layout.hpp`, line 987.*

foot_comp rejects a foot that touches any island rectangle. That is the only thing stopping a strap from landing on a shorted resonator's or flux line's centre conductor, which is in the same metal component as its pocket and would short the line mid-length. The connectivity gate cannot see this, because the foot lands on the 'right' net. The new layout_geometry strap case has no shorted line inside a pocket, and layout_presets asserts nothing about strap feet or spans on the real chips. The current code is correct: 0 feet on island metal and 0 spans over a qubit or junction region on AF-17 and BF-49.

**Evidence.** Mutant `foot_on_island` (`if (hit) return -1;` removed), probe feet.cpp: "AF-17: 122 straps, 47 with a foot on island metal (R: 30, F*: 17); conn.clean yes, DRC 0+0". The real `layout_geometry` group passes 26/26 (test_foot_on_island.log), and the layout_presets mask/connectivity assertions PASS on AF-17, BF-49 and GF-97.

**Suggested fix.** In layout_presets, assert that no strap foot intersects island metal and that no strap span intersects a qubit island or a junction region. Add a geometry case with a shorted line ending in a pocket, where the foot must avoid the line.

### [LOW] Stale or wrong comments about the new geometry and checks

*Lens: tests-docs. File: `include/qpuforge/layout.hpp`, line 1066.*

Three comments no longer match the code:
- qubit_control_lines' comment says the XY line ends "at x = cx + xy_end, y = cy - 156". It is drawn at y = cy - kXYLevel = cy - 279, as the design note correctly states.
- The ConnectivityReport header (line 1190) still describes the 0.4 junction rule ("lands on metal at both ends ..., on two different components").
- tests/test_qpuforge.cpp:1916 says "the only DRC markers left are the ground necks the connectivity guard kept (one per undone fill)". assemble_chip no longer sets fills_reverted (always 0), so the check at 1917 is simply space == 0. connectivity_json still exports `sliver_fills_undone`.

**Evidence.** layout.hpp:1066 vs `cpw_path({{cx + xy_end, cy - xy_level}, ...})` with kXYLevel = 279.0. layout.hpp:1190 vs the rule at 1442. `grep -n fills_reverted include/qpuforge/layout.hpp` finds only the declaration (1195) and the JSON export (1470).

**Suggested fix.** Correct the y coordinate in the comment and update the ConnectivityReport header to the island-to-plane rule. Rewrite the test comment, or compare against 0 directly, and drop or rename the always-zero field.
