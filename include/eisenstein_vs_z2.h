#ifndef EISENSTEIN_VS_Z2_H
#define EISENSTEIN_VS_Z2_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

/* ── Eisenstein integers ──────────────────────────────────────────── */

/**
 * Eisenstein integer: a + b*ω where ω = e^(2πi/3) = -1/2 + i√3/2
 *
 * Cartesian mapping: (a - b/2, b*√3/2)
 * Norm: a² - ab + b²
 *
 * The Eisenstein lattice provides the densest 2D packing, giving
 * better quantization (snap) properties than the square Z² lattice.
 */

/** Eisenstein integer lattice point (a, b) */
typedef struct {
    int64_t a;
    int64_t b;
} EisensteinInt;

/** Square lattice (Z²) integer point */
typedef struct {
    int64_t x;
    int64_t y;
} Z2Int;

/**
 * Convert Eisenstein integer to Cartesian coordinates.
 * @param a  Eisenstein a-component
 * @param b  Eisenstein b-component
 * @param x  [out] Cartesian x = a - b/2
 * @param y  [out] Cartesian y = b * √3/2
 */
static inline void eisenstein_to_cartesian(int64_t a, int64_t b, double *x, double *y) {
    if (x) *x = a - b / 2.0;
    if (y) *y = b * 0.8660254037844386; /* √3/2 */
}

/**
 * Convert Cartesian coordinates to Eisenstein coords (continuous, pre-rounding).
 * @param x  Cartesian x
 * @param y  Cartesian y
 * @param a  [out] Eisenstein a = x + b/2
 * @param b  [out] Eisenstein b = y / (√3/2)
 */
static inline void cartesian_to_eisenstein(double x, double y, double *a, double *b) {
    if (b) *b = y / 0.8660254037844386; /* y / (√3/2) */
    if (a) *a = x + (*b) / 2.0;
}

/**
 * Compute the Eisenstein norm: a² - ab + b².
 * This equals |a + bω|² and is always non-negative.
 * @param a  Eisenstein a-component
 * @param b  Eisenstein b-component
 * @return   Norm value (always >= 0)
 */
static inline int64_t eisenstein_norm(int64_t a, int64_t b) {
    return a * a - a * b + b * b;
}

/**
 * Snap a 2D point to the nearest Eisenstein integer lattice point.
 * @param x  Cartesian x coordinate
 * @param y  Cartesian y coordinate
 * @return   Nearest Eisenstein lattice point
 */
static inline EisensteinInt snap_to_eisenstein(double x, double y) {
    double a_cont, b_cont;
    cartesian_to_eisenstein(x, y, &a_cont, &b_cont);
    EisensteinInt result;
    result.a = (int64_t)round(a_cont);
    result.b = (int64_t)round(b_cont);
    return result;
}

/**
 * Snap a 2D point to nearest Z² (square lattice) point.
 * @param x  Cartesian x coordinate
 * @param y  Cartesian y coordinate
 * @return   Nearest square lattice point
 */
static inline Z2Int snap_to_z2(double x, double y) {
    Z2Int result;
    result.x = (int64_t)round(x);
    result.y = (int64_t)round(y);
    return result;
}

/* ── Metrics ──────────────────────────────────────────────────────── */

/**
 * Compute Euclidean distance between a point and its Eisenstein snap.
 * @param x  Cartesian x coordinate
 * @param y  Cartesian y coordinate
 * @return   Distance to nearest Eisenstein lattice point
 */
double eisenstein_snap_error(double x, double y);

/**
 * Compute Euclidean distance between a point and its Z² snap.
 * @param x  Cartesian x coordinate
 * @param y  Cartesian y coordinate
 * @return   Distance to nearest Z² lattice point
 */
double z2_snap_error(double x, double y);

/* ── Benchmark ────────────────────────────────────────────────────── */

/** Results from a lattice snap benchmark trial */
typedef struct {
    int n;                 /**< Number of points tested */
    double mean_error;     /**< Mean snap error */
    double max_error;      /**< Maximum snap error */
    double median_error;   /**< Median snap error */
    double std_error;      /**< Standard deviation of errors */
    int packing_unique;    /**< Number of unique snapped points */
    double packing_ratio;  /**< Fraction of unique snaps (unique/n) */
    double recovery_001;   /**< Fraction within error ≤ 0.01 */
    double recovery_01;    /**< Fraction within error ≤ 0.1 */
    double recovery_05;    /**< Fraction within error ≤ 0.5 */
    double elapsed_s;      /**< Wall-clock time in seconds (if timed) */
} TrialResult;

/**
 * Run a benchmark trial snapping points to the Eisenstein lattice.
 * @param points_x  Array of x coordinates
 * @param points_y  Array of y coordinates
 * @param n         Number of points
 * @return          Trial results with statistics
 */
TrialResult run_eisenstein_trial(const double *points_x, const double *points_y, int n);

/**
 * Run a benchmark trial snapping points to the Z² lattice.
 * @param points_x  Array of x coordinates
 * @param points_y  Array of y coordinates
 * @param n         Number of points
 * @return          Trial results with statistics
 */
TrialResult run_z2_trial(const double *points_x, const double *points_y, int n);

/**
 * Run a timed benchmark, measuring wall-clock elapsed time.
 * @param points_x  Array of x coordinates
 * @param points_y  Array of y coordinates
 * @param n         Number of points
 * @param use_eisenstein  1 for Eisenstein, 0 for Z²
 * @return          Trial results including elapsed_s
 */
TrialResult run_timed_trial(const double *points_x, const double *points_y, int n, int use_eisenstein);

/** Maximum theoretical snap error for Eisenstein lattice (covering radius ≈ 0.5774) */
static inline double max_theoretical_error_eisenstein(void) {
    return 1.0 / sqrt(3.0);
}

/** Maximum theoretical snap error for Z² lattice (covering radius ≈ 0.7071) */
static inline double max_theoretical_error_z2(void) {
    return 1.0 / sqrt(2.0);
}

/** Voronoi cell area for Eisenstein lattice (√3/2 ≈ 0.8660) */
static inline double voronoi_area_eisenstein(void) {
    return 0.8660254037844386;
}

/** Voronoi cell area for Z² lattice (1.0) */
static inline double voronoi_area_z2(void) {
    return 1.0;
}

/**
 * Compare two EisensteinInt for equality.
 * @return 1 if equal, 0 otherwise
 */
static inline int eisenstein_eq(EisensteinInt a, EisensteinInt b) {
    return a.a == b.a && a.b == b.b;
}

/**
 * Compare two Z2Int for equality.
 * @return 1 if equal, 0 otherwise
 */
static inline int z2_eq(Z2Int a, Z2Int b) {
    return a.x == b.x && a.y == b.y;
}

/**
 * Add two Eisenstein integers.
 */
static inline EisensteinInt eisenstein_add(EisensteinInt a, EisensteinInt b) {
    EisensteinInt r = { a.a + b.a, a.b + b.b };
    return r;
}

/**
 * Compute the Euclidean distance between two Eisenstein lattice points.
 */
double eisenstein_distance(EisensteinInt a, EisensteinInt b);

/**
 * Compute the Euclidean distance between two Z² lattice points.
 */
double z2_distance(Z2Int a, Z2Int b);

/**
 * Check if an Eisenstein integer is a unit (norm == 1).
 * The six Eisenstein units are: ±1, ±ω, ±ω²
 * @return 1 if unit, 0 otherwise
 */
static inline int eisenstein_is_unit(int64_t a, int64_t b) {
    return eisenstein_norm(a, b) == 1;
}

#endif /* EISENSTEIN_VS_Z2_H */
