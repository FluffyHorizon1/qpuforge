// Build-flag guard: the review programs must not fuse a*b+c into one FMA instruction.
// With contraction (-march=native / -mfma with GCC's default -ffp-contract=fast, or an FMA-by-default target),
// the Alpha-1 Gen-0 noise model changes in the last bits and the frame Monte Carlo takes different branches:
// Gen-0 eps 0.0120023693144 becomes 0.0122286733965 (8 threads, 100000 shots). CMakeLists.txt therefore builds every
// target with -ffp-contract=off; this test fails if that flag is lost and the compiler starts contracting.
//
// a = 1 + 2^-30:  a*a = 1 + 2^-29 + 2^-60.  Rounded product then subtract (no FMA): exactly 0.
// Fused multiply-add keeps the 2^-60 term.  The inputs come through volatile so nothing is constant-folded.
#include <cstdio>

static double mul_sub(double a, double b, double c) { return a * b - c; }

int main() {
    volatile double one = 1.0, tiny = 0x1p-30;
    const double a = one + tiny;
    const double c = one + 2 * tiny;
    const double r = mul_sub(a, a, c);
    if (r != 0.0) {
        std::printf("FAIL fp_contract: a*a-c = %a (expected 0): the compiler fused multiply and add; build with -ffp-contract=off\n", r);
        return 1;
    }
    std::printf("fp_contract: ok (no FMA contraction)\n");
    return 0;
}
