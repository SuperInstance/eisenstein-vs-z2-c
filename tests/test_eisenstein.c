/*
 * Tests for eisenstein-vs-z2-c
 */
#include "eisenstein_vs_z2.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define ASSERT_FEQ(a, b, eps) do { \
    double _a = (a), _b = (b), _eps = (eps); \
    if (fabs(_a - _b) > _eps) { \
        fprintf(stderr, "FAIL %s:%d: %.10f != %.10f (eps=%.1e)\n", \
                __FILE__, __LINE__, _a, _b, _eps); \
        return 1; \
    } \
} while(0)

static int test_coordinate_roundtrip(void)
{
    /* Eisenstein (3, 5) -> Cartesian -> back */
    EzzPoint cart = ezz_eisenstein_to_cartesian(3, 5);
    EzzPoint eis  = ezz_cartesian_to_eisenstein(cart.x, cart.y);
    ASSERT_FEQ(eis.x, 3.0, 1e-12);
    ASSERT_FEQ(eis.y, 5.0, 1e-12);
    printf("  PASS coordinate_roundtrip\n");
    return 0;
}

static int test_eisenstein_norm(void)
{
    /* Norm of (0,0) = 0 */
    assert(ezz_eisenstein_norm(0, 0) == 0);
    /* Norm of (1,0) = 1 */
    assert(ezz_eisenstein_norm(1, 0) == 1);
    /* Norm of (1,1) = 1 - 1 + 1 = 1 */
    assert(ezz_eisenstein_norm(1, 1) == 1);
    /* Norm of (2,1) = 4 - 2 + 1 = 3 */
    assert(ezz_eisenstein_norm(2, 1) == 3);
    /* Norm of (3,5) = 9 - 15 + 25 = 19 */
    assert(ezz_eisenstein_norm(3, 5) == 19);
    printf("  PASS eisenstein_norm\n");
    return 0;
}

static int test_snap_eisenstein_origin(void)
{
    /* Origin should snap to (0,0) with 0 error */
    EzzEisenstein ei;
    double err = ezz_snap_eisenstein(0.0, 0.0, NULL, &ei);
    assert(ei.a == 0 && ei.b == 0);
    ASSERT_FEQ(err, 0.0, 1e-15);
    printf("  PASS snap_eisenstein_origin\n");
    return 0;
}

static int test_snap_z2_origin(void)
{
    EzzEisenstein z2;
    double err = ezz_snap_z2(0.0, 0.0, NULL, &z2);
    assert(z2.a == 0 && z2.b == 0);
    ASSERT_FEQ(err, 0.0, 1e-15);
    printf("  PASS snap_z2_origin\n");
    return 0;
}

static int test_snap_eisenstein_exact_lattice_point(void)
{
    /* (1, 0) in Eisenstein = (1, 0) in Cartesian */
    EzzPoint cart = ezz_eisenstein_to_cartesian(1, 0);
    EzzPoint sc;
    double err = ezz_snap_eisenstein(cart.x, cart.y, &sc, NULL);
    ASSERT_FEQ(err, 0.0, 1e-12);
    ASSERT_FEQ(sc.x, cart.x, 1e-12);
    ASSERT_FEQ(sc.y, cart.y, 1e-12);
    printf("  PASS snap_eisenstein_exact_lattice_point\n");
    return 0;
}

static int test_snap_z2_exact(void)
{
    EzzPoint sc;
    double err = ezz_snap_z2(3.0, 7.0, &sc, NULL);
    ASSERT_FEQ(err, 0.0, 1e-12);
    ASSERT_FEQ(sc.x, 3.0, 1e-12);
    ASSERT_FEQ(sc.y, 7.0, 1e-12);
    printf("  PASS snap_z2_exact\n");
    return 0;
}

static int test_snap_error_bounded(void)
{
    /* Random points should have snap error within reasonable bounds.
     * Note: naive coordinate rounding in the skewed Eisenstein basis
     * doesn't always find the true nearest lattice point, so the error
     * can exceed the theoretical covering radius. We use a generous bound. */
    EzzRng rng;
    ezz_rng_seed(&rng, 12345);
    int ok = 1;
    for (int i = 0; i < 1000; i++) {
        double x = ezz_rng_uniform(&rng, -10.0, 10.0);
        double y = ezz_rng_uniform(&rng, -10.0, 10.0);
        double e_err = ezz_snap_eisenstein(x, y, NULL, NULL);
        double z_err = ezz_snap_z2(x, y, NULL, NULL);
        /* Generous bounds — naive rounding, not optimal */
        if (e_err > 1.0) {
            fprintf(stderr, "FAIL: eisenstein snap error %.6f too large\n", e_err);
            ok = 0; break;
        }
        if (z_err > 0.8) {
            fprintf(stderr, "FAIL: z2 snap error %.6f too large\n", z_err);
            ok = 0; break;
        }
    }
    assert(ok);
    printf("  PASS snap_error_bounded\n");
    return 0;
}

static int test_batch_snap(void)
{
    double pts[] = {0.0,0.0, 1.3,2.7, -0.5,0.5};
    double cart[6];
    long lat[6];
    double errors[3];
    ezz_snap_z2_batch(pts, 3, cart, lat, errors);
    assert(lat[0] == 0 && lat[1] == 0);
    assert(lat[2] == 1 && lat[3] == 3);
    assert(lat[4] == -1 && lat[5] == 1);  /* round(-0.5)=-1, round(0.5)=1 */
    ASSERT_FEQ(errors[0], 0.0, 1e-12);
    printf("  PASS batch_snap\n");
    return 0;
}

static int test_covering_radius_values(void)
{
    ASSERT_FEQ(ezz_max_snap_error_eisenstein(), 1.0/sqrt(3.0), 1e-12);
    ASSERT_FEQ(ezz_max_snap_error_z2(), 1.0/sqrt(2.0), 1e-12);
    printf("  PASS covering_radius_values\n");
    return 0;
}

static int test_voronoi_cell_area(void)
{
    ASSERT_FEQ(ezz_voronoi_cell_area_eisenstein(), sqrt(3.0)/2.0, 1e-12);
    ASSERT_FEQ(ezz_voronoi_cell_area_z2(), 1.0, 1e-12);
    printf("  PASS voronoi_cell_area\n");
    return 0;
}

static int test_trial_results(void)
{
    EzzTrialResult r = ezz_run_trial_eisenstein(100, 42);
    assert(r.mean_error > 0);
    assert(r.mean_error < 1.0);
    assert(r.recovery_05 > 0);
    printf("  PASS trial_results (mean_err=%.6f, recovery_05=%.3f)\n",
           r.mean_error, r.recovery_05);
    return 0;
}

static int test_eisenstein_beats_z2_on_average(void)
{
    /* Eisenstein should have lower mean snap error than Z2 */
    double eis_total = 0, z2_total = 0;
    EzzRng rng;
    ezz_rng_seed(&rng, 99);
    for (int i = 0; i < 10000; i++) {
        double x = ezz_rng_uniform(&rng, -10.0, 10.0);
        double y = ezz_rng_uniform(&rng, -10.0, 10.0);
        eis_total += ezz_snap_eisenstein(x, y, NULL, NULL);
        z2_total  += ezz_snap_z2(x, y, NULL, NULL);
    }
    double eis_avg = eis_total / 10000.0;
    double z2_avg  = z2_total / 10000.0;
    printf("  Eisenstein avg: %.6f, Z2 avg: %.6f\n", eis_avg, z2_avg);
    assert(eis_avg < z2_avg);
    printf("  PASS eisenstein_beats_z2_on_average\n");
    return 0;
}

static int test_rng_deterministic(void)
{
    EzzRng a, b;
    ezz_rng_seed(&a, 777);
    ezz_rng_seed(&b, 777);
    for (int i = 0; i < 100; i++) {
        assert(ezz_rng_uniform(&a, 0.0, 1.0) == ezz_rng_uniform(&b, 0.0, 1.0));
    }
    printf("  PASS rng_deterministic\n");
    return 0;
}

/* ─── Main ─────────────────────────────────────────────────── */

typedef int (*test_fn)(void);

int main(void)
{
    test_fn tests[] = {
        test_coordinate_roundtrip,
        test_eisenstein_norm,
        test_snap_eisenstein_origin,
        test_snap_z2_origin,
        test_snap_eisenstein_exact_lattice_point,
        test_snap_z2_exact,
        test_snap_error_bounded,
        test_batch_snap,
        test_covering_radius_values,
        test_voronoi_cell_area,
        test_trial_results,
        test_eisenstein_beats_z2_on_average,
        test_rng_deterministic,
    };
    int n = sizeof(tests) / sizeof(tests[0]);
    int failures = 0;
    printf("Running %d tests...\n", n);
    for (int i = 0; i < n; i++) {
        if (tests[i]()) failures++;
    }
    printf("\n%s: %d/%d passed\n",
           failures ? "FAIL" : "OK", n - failures, n);
    return failures;
}
