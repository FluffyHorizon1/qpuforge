#include <cmath>
#include <cstdio>
#include <vector>
#include <utility>
static double flux(const std::vector<std::pair<double,double>>& P) {   // closed in-plane filament loop -> flux through x[-34,-24] y[-23,-19] per A
    double phi = 0; const int g = 40, ns = 4000;
    for (int i = 0; i < g; ++i) for (int j = 0; j < g; ++j) {
        const double x = -34 + (i + 0.5) * 10.0 / g, y = -23 + (j + 0.5) * 4.0 / g; double b = 0;
        for (std::size_t s = 0; s + 1 < P.size(); ++s) for (int k = 0; k < ns; ++k) {
            const double ax = P[s].first, ay = P[s].second, bx = P[s + 1].first, by = P[s + 1].second, t = (k + 0.5) / ns;
            const double px = ax + (bx - ax) * t, py = ay + (by - ay) * t, dlx = (bx - ax) / ns, dly = (by - ay) / ns, rx = x - px, ry = y - py;
            b += (dlx * ry - dly * rx) / std::pow(rx * rx + ry * ry, 1.5);
        }
        phi += 1e-7 * b * 1e6 * (10.0 / g) * (4.0 / g) * 1e-12;
    }
    return phi;
}
int main() {
    const double Phi0 = 2.067833848e-15;
    const double a = flux({{-156, -156}, {-61, -61}, {-51, -76.5}, {-159, -76.5}, {-156, -156}});
    const double b = flux({{-156, -156}, {-61, -61}, {-65, -51}, {-65, -161}, {-156, -156}});
    std::printf("Q0-type pocket strap (left):  M = %.4f pH, I(Phi0/2) = %.0f mA\nQ14-type pocket strap (down): M = %.4f pH, I(Phi0/2) = %.0f mA\n", a * 1e12, Phi0 / 2 / std::abs(a) * 1e3, b * 1e12, Phi0 / 2 / std::abs(b) * 1e3);
}
