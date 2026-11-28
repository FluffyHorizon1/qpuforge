// flux through the SQUID region (x[-34,-24] y[-23,-19]) from each open filament segment of the flux-line current path
#include <cmath>
#include <cstdio>
#include <vector>
int main() {
    struct Seg { const char* name; double ax, ay, bx, by; };
    const std::vector<Seg> segs{{"stair centre line (-156,-156)->(-61,-61)", -156, -156, -61, -61},
                                {"pocket: short (-61,-61) -> strap foot (-51,-76.5)", -61, -61, -51, -76.5},
                                {"bridge (-51,-76.5) -> (-159,-76.5)", -51, -76.5, -159, -76.5},
                                {"return (-159,-76.5) -> (-156,-156)", -159, -76.5, -156, -156},
                                {"worst case: pocket current hugging the tongue base (-60,-54)->(-38,-54)", -60, -54, -38, -54},
                                {"reference: current up the tongue centre (-39,-54)->(-39,-18)", -39, -54, -39, -18}};
    for (const auto& s : segs) {
        double phi = 0; const int g = 40, ns = 4000;
        for (int i = 0; i < g; ++i) for (int j = 0; j < g; ++j) {
            const double x = -34 + (i + 0.5) * 10.0 / g, y = -23 + (j + 0.5) * 4.0 / g;
            double b = 0;
            for (int k = 0; k < ns; ++k) {
                const double t = (k + 0.5) / ns, px = s.ax + (s.bx - s.ax) * t, py = s.ay + (s.by - s.ay) * t;
                const double dlx = (s.bx - s.ax) / ns, dly = (s.by - s.ay) / ns, rx = x - px, ry = y - py, r2 = rx * rx + ry * ry;
                b += (dlx * ry - dly * rx) / std::pow(r2, 1.5);
            }
            phi += 1e-7 * b * 1e6 * (10.0 / g) * (4.0 / g) * 1e-12;
        }
        std::printf("%-75s M = %+.3f pH\n", s.name, phi * 1e12);
    }
    std::printf("Phi0/2 / (0.2 pH) = %.1f mA\n", 2.067833848e-15 / 2 / 0.2e-12 * 1e3);
}
