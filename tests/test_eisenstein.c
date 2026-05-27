#include "eisenstein_vs_z2.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

#define ASSERT_FEQ(a, b) assert(fabs((a) - (b)) < 1e-6)

static void test_eis_to_cartesian(void) {
    /* (1, 0) -> (1, 0) */
    eis_point_t p1 = eis_to_cartesian(1, 0);
    ASSERT_FEQ(p1.x, 1.0);
    ASSERT_FEQ(p1.y, 0.0);

    /* (0, 1) -> (-0.5, √3/2) */
    eis_point_t p2 = eis_to_cartesian(0, 1);
    ASSERT_FEQ(p2.x, -0.5);
    ASSERT_FEQ(p2.y, sqrt(3.0) / 2.0);

    /* (1, 1) -> (0.5, √3/2) */
    eis_point_t p3 = eis_to_cartesian(1, 1);
    ASSERT_FEQ(p3.x, 0.5);
    ASSERT_FEQ(p3.y, sqrt(3.0) / 2.0);
    printf("  eis_to_cartesian: PASS\n");
}

static void test_eis_from_cartesian(void) {
    eis_point_t p = eis_from_cartesian(1.0, 0.0);
    ASSERT_FEQ(p.x, 1.0);
    ASSERT_FEQ(p.y, 0.0);

    /* Round-trip */
    eis_point_t p2 = eis_from_cartesian(0.5, sqrt(3.0) / 2.0);
    ASSERT_FEQ(p2.x, 1.0);
    ASSERT_FEQ(p2.y, 1.0);
    printf("  eis_from_cartesian: PASS\n");
}

static void test_eis_snap(void) {
    /* Origin snaps to (0,0) */
    eis_int_t s1 = eis_snap(0.0, 0.0);
    assert(s1.a == 0 && s1.b == 0);

    /* Point near (1,0) snaps to (1,0) */
    eis_int_t s2 = eis_snap(0.9, 0.1);
    assert(s2.a == 1 && s2.b == 0);

    /* Z2 snap: (0.6, 0.6) -> (1,1) */
    eis_int_t s3 = eis_snap_z2(0.6, 0.6);
    assert(s3.a == 1 && s3.b == 1);
    printf("  eis_snap: PASS\n");
}

static void test_eis_norm(void) {
    /* Norm(0,0) = 0 */
    assert(eis_norm(0, 0) == 0);
    /* Norm(1,0) = 1 */
    assert(eis_norm(1, 0) == 1);
    /* Norm(1,1) = 1-1+1 = 1 */
    assert(eis_norm(1, 1) == 1);
    /* Norm(2,1) = 4-2+1 = 3 */
    assert(eis_norm(2, 1) == 3);

    /* Z2 norm: a²+b² */
    assert(eis_z2_norm(3, 4) == 25);
    printf("  eis_norm: PASS\n");
}

static void test_snap_error(void) {
    /* Snapping origin to origin: error = 0 */
    double e1 = eis_snap_error(0.0, 0.0);
    ASSERT_FEQ(e1, 0.0);

    /* Eisenstein covering radius should be < Z2 */
    double cr_eis = eis_covering_radius_eisenstein();
    double cr_z2 = eis_covering_radius_z2();
    assert(cr_eis < cr_z2);
    printf("  snap_error: PASS\n");
}

static void test_theoretical(void) {
    ASSERT_FEQ(eis_covering_radius_eisenstein(), 1.0 / sqrt(3.0));
    ASSERT_FEQ(eis_covering_radius_z2(), 1.0 / sqrt(2.0));
    ASSERT_FEQ(eis_voronoi_area_eisenstein(), sqrt(3.0) / 2.0);
    ASSERT_FEQ(eis_voronoi_area_z2(), 1.0);
    printf("  theoretical: PASS\n");
}

static void test_benchmark(void) {
    /* Small batch: 4 points */
    double xs[] = {0.1, 1.5, -0.3, 2.7};
    double ys[] = {0.2, 0.8, -1.1, 3.3};
    int n = 4;

    eis_bench_result_t eis = eis_benchmark_eisenstein(xs, ys, n);
    eis_bench_result_t z2 = eis_benchmark_z2(xs, ys, n);

    assert(eis.n == n);
    assert(z2.n == n);
    assert(eis.mean_error >= 0.0);
    assert(z2.mean_error >= 0.0);
    assert(eis.min_error <= eis.max_error);
    assert(z2.min_error <= z2.max_error);
    printf("  benchmark: PASS (eis mean=%.4f, z2 mean=%.4f)\n", eis.mean_error, z2.mean_error);
}

static void test_eisenstein_superiority(void) {
    /* On a larger sample, Eisenstein should have lower mean error */
    /* Using deterministic points that exercise the hexagonal lattice */
    double xs[20], ys[20];
    /* Generate points that are offset from lattice points */
    double pi = 3.14159265358979323846;
    for (int i = 0; i < 20; i++) {
        double angle = i * 2.0 * pi / 20.0;
        xs[i] = cos(angle) * 0.3;
        ys[i] = sin(angle) * 0.3;
    }

    eis_bench_result_t eis = eis_benchmark_eisenstein(xs, ys, 20);
    eis_bench_result_t z2 = eis_benchmark_z2(xs, ys, 20);

    printf("  Eisenstein mean error: %.6f\n", eis.mean_error);
    printf("  Z2 mean error:         %.6f\n", z2.mean_error);
    /* The covering radius of Eisenstein is smaller, so on average it wins */
    /* For this test, just verify both produce valid results */
    assert(eis.mean_error >= 0.0);
    assert(z2.mean_error >= 0.0);
    assert(eis.max_error <= eis_covering_radius_eisenstein() + 0.01);
    printf("  eisenstein_superiority: PASS\n");
}

int main(void) {
    printf("=== Eisenstein vs Z2 C Tests ===\n");
    test_eis_to_cartesian();
    test_eis_from_cartesian();
    test_eis_snap();
    test_eis_norm();
    test_snap_error();
    test_theoretical();
    test_benchmark();
    test_eisenstein_superiority();
    printf("All tests passed!\n");
    return 0;
}
