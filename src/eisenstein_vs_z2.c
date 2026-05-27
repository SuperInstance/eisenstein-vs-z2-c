#include "eisenstein_vs_z2.h"
#include <math.h>

/* ===================================================================
 * Conversion functions
 * =================================================================== */

eis_point_t eis_to_cartesian(long long a, long long b) {
    eis_point_t p;
    p.x = a - b / 2.0;
    p.y = b * sqrt(3.0) / 2.0;
    return p;
}

eis_point_t eis_from_cartesian(double x, double y) {
    eis_point_t p;
    double sqrt3_2 = sqrt(3.0) / 2.0;
    p.y = y / sqrt3_2;
    p.x = x + p.y / 2.0;
    return p;
}

/* ===================================================================
 * Snapping
 * =================================================================== */

eis_int_t eis_snap(double x, double y) {
    eis_point_t cont = eis_from_cartesian(x, y);
    eis_int_t result;
    result.a = (long long)floor(cont.x + 0.5);
    result.b = (long long)floor(cont.y + 0.5);
    return result;
}

eis_int_t eis_snap_z2(double x, double y) {
    eis_int_t result;
    result.a = (long long)floor(x + 0.5);
    result.b = (long long)floor(y + 0.5);
    return result;
}

/* ===================================================================
 * Norms
 * =================================================================== */

long long eis_norm(long long a, long long b) {
    return a * a - a * b + b * b;
}

long long eis_z2_norm(long long a, long long b) {
    return a * a + b * b;
}

/* ===================================================================
 * Distance / error
 * =================================================================== */

double eis_distance(eis_point_t a, eis_point_t b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

double eis_snap_error(double x, double y) {
    eis_int_t snapped = eis_snap(x, y);
    eis_point_t sp = eis_to_cartesian(snapped.a, snapped.b);
    eis_point_t orig = {x, y};
    return eis_distance(orig, sp);
}

double eis_z2_snap_error(double x, double y) {
    eis_int_t snapped = eis_snap_z2(x, y);
    eis_point_t sp = {(double)snapped.a, (double)snapped.b};
    eis_point_t orig = {x, y};
    return eis_distance(orig, sp);
}

/* ===================================================================
 * Theoretical properties
 * =================================================================== */

double eis_covering_radius_eisenstein(void) {
    return 1.0 / sqrt(3.0);
}

double eis_covering_radius_z2(void) {
    return 1.0 / sqrt(2.0);
}

double eis_voronoi_area_eisenstein(void) {
    return sqrt(3.0) / 2.0;
}

double eis_voronoi_area_z2(void) {
    return 1.0;
}

/* ===================================================================
 * Batch benchmarks
 * =================================================================== */

eis_bench_result_t eis_benchmark_eisenstein(const double *xs, const double *ys, int n) {
    eis_bench_result_t r = {0, 1e30, 0, n};
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double err = eis_snap_error(xs[i], ys[i]);
        sum += err;
        if (err < r.min_error) r.min_error = err;
        if (err > r.max_error) r.max_error = err;
    }
    r.mean_error = sum / n;
    return r;
}

eis_bench_result_t eis_benchmark_z2(const double *xs, const double *ys, int n) {
    eis_bench_result_t r = {0, 1e30, 0, n};
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double err = eis_z2_snap_error(xs[i], ys[i]);
        sum += err;
        if (err < r.min_error) r.min_error = err;
        if (err > r.max_error) r.max_error = err;
    }
    r.mean_error = sum / n;
    return r;
}
