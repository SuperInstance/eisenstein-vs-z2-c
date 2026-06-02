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

static int test_count = 0;
#define TEST(name) do { printf("  ✓ %s\n", name); test_count++; } while(0)

/* ── Original tests (9) ────────────────────────────────────────────── */

static void test_eisenstein_roundtrip(void) {
    double x, y;
    eisenstein_to_cartesian(3, 2, &x, &y);
    assert(fabs(x - 2.0) < 1e-10);
    assert(fabs(y - 1.7320508076) < 1e-10);
    EisensteinInt snapped = snap_to_eisenstein(x, y);
    assert(snapped.a == 3);
    assert(snapped.b == 2);
    TEST("Eisenstein roundtrip (3,2)");
}

static void test_z2_snap(void) {
    Z2Int z = snap_to_z2(1.4, 2.6);
    assert(z.x == 1);
    assert(z.y == 3);
    TEST("Z² snap (1.4, 2.6) → (1, 3)");
}

static void test_eisenstein_snap(void) {
    EisensteinInt e = snap_to_eisenstein(0.1, 0.05);
    assert(e.a == 0);
    assert(e.b == 0);
    TEST("Eisenstein snap near origin");
}

static void test_eisenstein_norm(void) {
    assert(eisenstein_norm(1, 0) == 1);
    assert(eisenstein_norm(1, 1) == 1);
    assert(eisenstein_norm(2, 1) == 3);
    assert(eisenstein_norm(0, 0) == 0);
    TEST("Eisenstein norm");
}

static void test_covering_radius(void) {
    double e = max_theoretical_error_eisenstein();
    double z = max_theoretical_error_z2();
    assert(e < z);
    assert(fabs(e - 0.5774) < 0.001);
    assert(fabs(z - 0.7071) < 0.001);
    TEST("Covering radii: Eisenstein < Z²");
}

static void test_voronoi_areas(void) {
    double ve = voronoi_area_eisenstein();
    double vz = voronoi_area_z2();
    assert(fabs(ve - 0.8660) < 0.001);
    assert(fabs(vz - 1.0) < 0.001);
    TEST("Voronoi areas");
}

static void test_snap_error_eisenstein_better(void) {
    int n = 10000;
    double *px = malloc((size_t)n * sizeof(double));
    double *py = malloc((size_t)n * sizeof(double));

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

    assert(e_result.mean_error < z_result.mean_error);
    free(px);
    free(py);
    TEST("Eisenstein mean error < Z² mean error");
}

static void test_packing_eisenstein_better(void) {
    int n = 5000;
    double *px = malloc((size_t)n * sizeof(double));
    double *py = malloc((size_t)n * sizeof(double));

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
    assert(e_result.packing_ratio >= z_result.packing_ratio);
    free(px);
    free(py);
    TEST("Eisenstein packing ratio >= Z²");
}

static void test_recovery_rates(void) {
    int n = 1000;
    double *px = malloc((size_t)n * sizeof(double));
    double *py = malloc((size_t)n * sizeof(double));

    rng_state = 777;
    for (int i = 0; i < n; i++) {
        px[i] = rng_uniform(-5, 5);
        py[i] = rng_uniform(-5, 5);
    }
    TrialResult e = run_eisenstein_trial(px, py, n);
    assert(e.recovery_05 > 0.5);
    free(px);
    free(py);
    TEST("Recovery rates @0.5 > 0.5");
}

/* ── New tests (target: 18+) ──────────────────────────────────────── */

static void test_eisenstein_units(void) {
    /* The 6 Eisenstein units have norm 1 */
    assert(eisenstein_is_unit(1, 0));
    assert(eisenstein_is_unit(-1, 0));
    assert(eisenstein_is_unit(0, 1));
    assert(eisenstein_is_unit(0, -1));
    assert(eisenstein_is_unit(1, 1));
    assert(eisenstein_is_unit(-1, -1));
    /* Non-units */
    assert(!eisenstein_is_unit(2, 0));
    assert(!eisenstein_is_unit(1, 2));
    assert(!eisenstein_is_unit(0, 0));
    TEST("Eisenstein units (6 units, norm=1)");
}

static void test_eisenstein_eq(void) {
    EisensteinInt a = {3, 5};
    EisensteinInt b = {3, 5};
    EisensteinInt c = {3, 6};
    assert(eisenstein_eq(a, b));
    assert(!eisenstein_eq(a, c));
    TEST("Eisenstein equality");
}

static void test_z2_eq(void) {
    Z2Int a = {1, 2};
    Z2Int b = {1, 2};
    Z2Int c = {2, 1};
    assert(z2_eq(a, b));
    assert(!z2_eq(a, c));
    TEST("Z² equality");
}

static void test_eisenstein_add(void) {
    EisensteinInt a = {1, 2};
    EisensteinInt b = {3, -1};
    EisensteinInt c = eisenstein_add(a, b);
    assert(c.a == 4 && c.b == 1);
    TEST("Eisenstein addition");
}

static void test_eisenstein_distance(void) {
    EisensteinInt a = {0, 0};
    EisensteinInt b = {1, 0};
    double d = eisenstein_distance(a, b);
    assert(fabs(d - 1.0) < 1e-10);
    /* Distance to self is zero */
    assert(fabs(eisenstein_distance(a, a)) < 1e-10);
    TEST("Eisenstein distance");
}

static void test_z2_distance(void) {
    Z2Int a = {0, 0};
    Z2Int b = {3, 4};
    double d = z2_distance(a, b);
    assert(fabs(d - 5.0) < 1e-10);
    TEST("Z² distance (3-4-5)");
}

static void test_null_pointer_handling(void) {
    /* Trial with NULL pointers should return zeroed result */
    TrialResult r = run_eisenstein_trial(NULL, NULL, 100);
    assert(r.n == 100);
    assert(r.mean_error == 0.0);

    TrialResult r2 = run_z2_trial(NULL, NULL, 100);
    assert(r2.n == 100);
    TEST("NULL pointer handling");
}

static void test_zero_points(void) {
    TrialResult r = run_eisenstein_trial(NULL, NULL, 0);
    assert(r.n == 0);
    TEST("Zero points trial");
}

static void test_negative_points(void) {
    double px[] = {0.0};
    double py[] = {0.0};
    TrialResult r = run_eisenstein_trial(px, py, -5);
    assert(r.n == -5);
    TEST("Negative count handling");
}

static void test_exact_lattice_snap(void) {
    /* Points exactly on lattice should have zero error */
    EisensteinInt ei = {5, 3};
    double x, y;
    eisenstein_to_cartesian(ei.a, ei.b, &x, &y);
    double err = eisenstein_snap_error(x, y);
    assert(fabs(err) < 1e-10);

    Z2Int zi = {7, 11};
    double zerr = z2_snap_error((double)zi.x, (double)zi.y);
    assert(fabs(zerr) < 1e-10);
    TEST("Exact lattice snap (zero error)");
}

static void test_large_coordinates(void) {
    EisensteinInt e = snap_to_eisenstein(1e6, 1e6);
    assert(e.a > 0 || e.b > 0);
    Z2Int z = snap_to_z2(1e8, -1e8);
    assert(z.x == 100000000);
    assert(z.y == -100000000);
    TEST("Large coordinates");
}

static void test_norm_symmetry(void) {
    /* Norm(a,b) == Norm(b,a) for Eisenstein? No: Norm(a,b)=a²-ab+b², Norm(b,a)=b²-ab+a² = same */
    assert(eisenstein_norm(3, 7) == eisenstein_norm(7, 3));
    assert(eisenstein_norm(-2, 5) == eisenstein_norm(5, -2));
    TEST("Eisenstein norm symmetry");
}

static void test_timed_trial(void) {
    int n = 10000;
    double *px = malloc((size_t)n * sizeof(double));
    double *py = malloc((size_t)n * sizeof(double));

    rng_state = 42;
    for (int i = 0; i < n; i++) {
        px[i] = rng_uniform(-100, 100);
        py[i] = rng_uniform(-100, 100);
    }
    TrialResult r = run_timed_trial(px, py, n, 1);
    assert(r.n == n);
    assert(r.elapsed_s >= 0.0);
    assert(r.mean_error > 0.0);
    printf("    Eisenstein benchmark: %d pts in %.4f s\n", n, r.elapsed_s);

    TrialResult r2 = run_timed_trial(px, py, n, 0);
    assert(r2.n == n);
    assert(r2.elapsed_s >= 0.0);
    printf("    Z² benchmark: %d pts in %.4f s\n", n, r2.elapsed_s);

    free(px);
    free(py);
    TEST("Timed benchmark trial");
}

static void test_cartesian_null_pointers(void) {
    /* Should not crash with NULL output pointers */
    eisenstein_to_cartesian(3, 2, NULL, NULL);
    cartesian_to_eisenstein(2.0, 1.732, NULL, NULL);
    TEST("NULL output pointer safety");
}

static void test_negative_eisenstein(void) {
    EisensteinInt e = snap_to_eisenstein(-2.0, -1.732);
    assert(e.a <= 0 || e.b <= 0);
    double err = eisenstein_snap_error(-2.0, -1.732);
    assert(err < 1.0);
    TEST("Negative coordinate snap");
}

/* ── Main ──────────────────────────────────────────────────────────── */

int main(void) {
    printf("=== eisenstein-vs-z2 tests ===\n\n");

    /* Original tests */
    test_eisenstein_roundtrip();
    test_z2_snap();
    test_eisenstein_snap();
    test_eisenstein_norm();
    test_covering_radius();
    test_voronoi_areas();
    test_snap_error_eisenstein_better();
    test_packing_eisenstein_better();
    test_recovery_rates();

    /* New tests */
    test_eisenstein_units();
    test_eisenstein_eq();
    test_z2_eq();
    test_eisenstein_add();
    test_eisenstein_distance();
    test_z2_distance();
    test_null_pointer_handling();
    test_zero_points();
    test_negative_points();
    test_exact_lattice_snap();
    test_large_coordinates();
    test_norm_symmetry();
    test_timed_trial();
    test_cartesian_null_pointers();
    test_negative_eisenstein();

    printf("\n✅ All %d tests passed!\n", test_count);
    return 0;
}
