#ifndef EISENSTEIN_VS_Z2_H
#define EISENSTEIN_VS_Z2_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

/* ── Eisenstein integers ──────────────────────────────────────────── */

/* Eisenstein integer: a + b*ω where ω = e^(2πi/3) = -1/2 + i√3/2
 * Cartesian mapping: (a - b/2, b*√3/2)
 * Norm: a² - ab + b²
 */

typedef struct {
    int64_t a;
    int64_t b;
} EisensteinInt;

/* Convert Eisenstein integer to Cartesian coordinates */
static inline void eisenstein_to_cartesian(int64_t a, int64_t b, double *x, double *y) {
    *x = a - b / 2.0;
    *y = b * 0.8660254037844386; /* √3/2 */
}

/* Convert Cartesian to Eisenstein coords (continuous, pre-rounding) */
static inline void cartesian_to_eisenstein(double x, double y, double *a, double *b) {
    *b = y / 0.8660254037844386; /* y / (√3/2) */
    *a = x + (*b) / 2.0;
}

/* Eisenstein norm: a² - ab + b² */
static inline int64_t eisenstein_norm(int64_t a, int64_t b) {
    return a * a - a * b + b * b;
}

/* Snap a 2D point to the nearest Eisenstein integer lattice point */
static inline EisensteinInt snap_to_eisenstein(double x, double y) {
    double a_cont, b_cont;
    cartesian_to_eisenstein(x, y, &a_cont, &b_cont);
    EisensteinInt result;
    result.a = (int64_t)round(a_cont);
    result.b = (int64_t)round(b_cont);
    return result;
}

/* Snap a 2D point to nearest Z² (square lattice) point */
typedef struct {
    int64_t x;
    int64_t y;
} Z2Int;

static inline Z2Int snap_to_z2(double x, double y) {
    Z2Int result;
    result.x = (int64_t)round(x);
    result.y = (int64_t)round(y);
    return result;
}

/* ── Metrics ──────────────────────────────────────────────────────── */

/* Euclidean distance between a point and its Eisenstein snap */
double eisenstein_snap_error(double x, double y);

/* Euclidean distance between a point and its Z² snap */
double z2_snap_error(double x, double y);

/* ── Benchmark ────────────────────────────────────────────────────── */

typedef struct {
    int n;
    double mean_error;
    double max_error;
    double median_error;
    double std_error;
    int packing_unique;
    double packing_ratio;
    double recovery_001;
    double recovery_01;
    double recovery_05;
    double elapsed_s;
} TrialResult;

/* Run a benchmark trial comparing Eisenstein vs Z² snapping */
TrialResult run_eisenstein_trial(const double *points_x, const double *points_y, int n);
TrialResult run_z2_trial(const double *points_x, const double *points_y, int n);

/* Maximum theoretical snap error (covering radius) */
static inline double max_theoretical_error_eisenstein(void) {
    return 1.0 / sqrt(3.0); /* ≈ 0.5774 */
}

static inline double max_theoretical_error_z2(void) {
    return 1.0 / sqrt(2.0); /* ≈ 0.7071 */
}

/* Voronoi cell area */
static inline double voronoi_area_eisenstein(void) {
    return 0.8660254037844386; /* √3/2 */
}

static inline double voronoi_area_z2(void) {
    return 1.0;
}

#endif /* EISENSTEIN_VS_Z2_H */
