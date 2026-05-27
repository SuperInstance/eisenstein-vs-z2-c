#include "eisenstein_vs_z2.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* ── Simple LCG PRNG (deterministic for reproducibility) ──────────── */
static uint64_t rng_state = 42;

static double rng_uniform(double lo, double hi) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    double t = (rng_state >> 11) / (double)(1ULL << 53);
    return lo + t * (hi - lo);
}

/* ── Tests ─────────────────────────────────────────────────────────── */

static void test_eisenstein_roundtrip(void) {
    /* Known Eisenstein integer (3, 2) → Cartesian → back */
    double x, y;
    eisenstein_to_cartesian(3, 2, &x, &y);
    /* x = 3 - 2/2 = 2.0, y = 2 * √3/2 = √3 ≈ 1.7321 */
    assert(fabs(x - 2.0) < 1e-10);
    assert(fabs(y - 1.7320508076) < 1e-10);

    /* Round-trip */
    EisensteinInt snapped = snap_to_eisenstein(x, y);
    assert(snapped.a == 3);
    assert(snapped.b == 2);
    printf("  ✓ Eisenstein roundtrip (3,2)\n");
}

static void test_z2_snap(void) {
    Z2Int z = snap_to_z2(1.4, 2.6);
    assert(z.x == 1);
    assert(z.y == 3);
    printf("  ✓ Z² snap (1.4, 2.6) → (1, 3)\n");
}

static void test_eisenstein_snap(void) {
    /* Point close to origin should snap to (0,0) */
    EisensteinInt e = snap_to_eisenstein(0.1, 0.05);
    assert(e.a == 0);
    assert(e.b == 0);
    printf("  ✓ Eisenstein snap near origin\n");
}

static void test_eisenstein_norm(void) {
    /* Norm of (1,0) should be 1 */
    assert(eisenstein_norm(1, 0) == 1);
    /* Norm of (1,1) should be 1 - 1 + 1 = 1 */
    assert(eisenstein_norm(1, 1) == 1);
    /* Norm of (2,1) should be 4 - 2 + 1 = 3 */
    assert(eisenstein_norm(2, 1) == 3);
    printf("  ✓ Eisenstein norm\n");
}

static void test_covering_radius(void) {
    /* Eisenstein covering radius should be < Z² covering radius */
    double e = max_theoretical_error_eisenstein();
    double z = max_theoretical_error_z2();
    assert(e < z);
    assert(fabs(e - 0.5774) < 0.001);
    assert(fabs(z - 0.7071) < 0.001);
    printf("  ✓ Covering radii: Eisenstein (%.4f) < Z² (%.4f)\n", e, z);
}

static void test_voronoi_areas(void) {
    double ve = voronoi_area_eisenstein();
    double vz = voronoi_area_z2();
    assert(fabs(ve - 0.8660) < 0.001);
    assert(fabs(vz - 1.0) < 0.001);
    printf("  ✓ Voronoi areas: Eisenstein (%.4f), Z² (%.4f)\n", ve, vz);
}

static void test_snap_error_eisenstein_better(void) {
    /* Eisenstein should have lower mean error over random points */
    int n = 10000;
    double *px = malloc(n * sizeof(double));
    double *py = malloc(n * sizeof(double));

    rng_state = 12345;
    for (int i = 0; i < n; i++) {
        px[i] = rng_uniform(-10, 10);
        py[i] = rng_uniform(-10, 10);
    }

    TrialResult e_result = run_eisenstein_trial(px, py, n);
    rng_state = 12345;
    for (int i = 0; i < n; i++) {
        px[i] = rng_uniform(-10, 10);
        py[i] = rng_uniform(-10, 10);
    }
    TrialResult z_result = run_z2_trial(px, py, n);

    printf("  ✓ Eisenstein mean error: %.6f\n", e_result.mean_error);
    printf("  ✓ Z² mean error:         %.6f\n", z_result.mean_error);
    printf("  ✓ Eisenstein / Z² ratio: %.4f\n", e_result.mean_error / z_result.mean_error);
    assert(e_result.mean_error < z_result.mean_error);

    free(px);
    free(py);
}

static void test_packing_eisenstein_better(void) {
    /* Eisenstein should have better packing (more unique snaps) */
    int n = 5000;
    double *px = malloc(n * sizeof(double));
    double *py = malloc(n * sizeof(double));

    rng_state = 99;
    for (int i = 0; i < n; i++) {
        px[i] = rng_uniform(-10, 10);
        py[i] = rng_uniform(-10, 10);
    }

    TrialResult e_result = run_eisenstein_trial(px, py, n);
    rng_state = 99;
    for (int i = 0; i < n; i++) {
        px[i] = rng_uniform(-10, 10);
        py[i] = rng_uniform(-10, 10);
    }
    TrialResult z_result = run_z2_trial(px, py, n);

    printf("  ✓ Eisenstein packing ratio: %.4f\n", e_result.packing_ratio);
    printf("  ✓ Z² packing ratio:         %.4f\n", z_result.packing_ratio);
    assert(e_result.packing_ratio >= z_result.packing_ratio);

    free(px);
    free(py);
}

static void test_recovery_rates(void) {
    int n = 1000;
    double *px = malloc(n * sizeof(double));
    double *py = malloc(n * sizeof(double));

    rng_state = 777;
    for (int i = 0; i < n; i++) {
        px[i] = rng_uniform(-5, 5);
        py[i] = rng_uniform(-5, 5);
    }

    TrialResult e = run_eisenstein_trial(px, py, n);
    printf("  ✓ Eisenstein recovery @0.5: %.4f\n", e.recovery_05);
    assert(e.recovery_05 > 0.5); /* Most points should be within 0.5 */

    free(px);
    free(py);
}

int main(void) {
    printf("=== eisenstein-vs-z2 tests ===\n\n");

    test_eisenstein_roundtrip();
    test_z2_snap();
    test_eisenstein_snap();
    test_eisenstein_norm();
    test_covering_radius();
    test_voronoi_areas();
    test_snap_error_eisenstein_better();
    test_packing_eisenstein_better();
    test_recovery_rates();

    printf("\n✅ All tests passed!\n");
    return 0;
}
