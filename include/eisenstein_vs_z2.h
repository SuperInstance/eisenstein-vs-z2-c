#ifndef EISENSTEIN_VS_Z2_H
#define EISENSTEIN_VS_Z2_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ─── Types ─────────────────────────────────────────────────── */

typedef struct {
    double x, y;
} EzzPoint;

typedef struct {
    long a, b;  /* Eisenstein coordinates */
} EzzEisenstein;

typedef struct {
    double mean_error;
    double std_error;
    double max_error;
    double median_error;
    double p25;
    double p75;
    double p95;
    double p99;
    size_t packing_unique;
    double packing_ratio;      /* unique / total */
    double recovery_001;       /* fraction within 0.01 */
    double recovery_01;        /* fraction within 0.1  */
    double recovery_05;        /* fraction within 0.5  */
    double elapsed_s;
} EzzTrialResult;

/* ─── Coordinate conversion ────────────────────────────────── */

/* Eisenstein (a,b) -> Cartesian (x,y) */
EzzPoint ezz_eisenstein_to_cartesian(long a, long b);

/* Cartesian (x,y) -> continuous Eisenstein (a,b) */
EzzPoint ezz_cartesian_to_eisenstein(double x, double y);

/* ─── Snapping ─────────────────────────────────────────────── */

/* Snap a single point to the nearest Eisenstein integer lattice point.
   Stores the snapped Cartesian point in *out_cart and Eisenstein coords
   in *out_eis. Returns the Euclidean snap error. */
double ezz_snap_eisenstein(double x, double y,
                           EzzPoint *out_cart, EzzEisenstein *out_eis);

/* Snap a single point to the nearest Z2 (square) lattice point. */
double ezz_snap_z2(double x, double y,
                   EzzPoint *out_cart, EzzEisenstein *out_z2);

/* Snap an array of n points in-place.
   points_xy is [x0,y0, x1,y1, ...] (length 2*n).
   out_cartesian and out_lattice_coords are similarly packed.
   errors_out receives n snap errors (may be NULL). */
void ezz_snap_eisenstein_batch(const double *points_xy, size_t n,
                               double *out_cartesian,
                               long *out_lattice_coords,
                               double *errors_out);

void ezz_snap_z2_batch(const double *points_xy, size_t n,
                       double *out_cartesian,
                       long *out_lattice_coords,
                       double *errors_out);

/* ─── Norms and metrics ────────────────────────────────────── */

/* Eisenstein norm: a^2 - a*b + b^2 */
long ezz_eisenstein_norm(long a, long b);

/* Maximum theoretical snap error (covering radius) for each lattice */
double ezz_max_snap_error_eisenstein(void);  /* 1/sqrt(3) */
double ezz_max_snap_error_z2(void);          /* 1/sqrt(2) */

/* Voronoi cell area */
double ezz_voronoi_cell_area_eisenstein(void);  /* sqrt(3)/2 */
double ezz_voronoi_cell_area_z2(void);          /* 1.0 */

/* ─── Benchmark ────────────────────────────────────────────── */

/* Run a single benchmark trial. Points are uniform in [-10,10]^2.
   Uses a deterministic PRNG seeded by seed. */
EzzTrialResult ezz_run_trial_eisenstein(size_t n_points, uint64_t seed);
EzzTrialResult ezz_run_trial_z2(size_t n_points, uint64_t seed);

/* ─── Utility ──────────────────────────────────────────────── */

/* Simple xoshiro256** PRNG */
typedef struct { uint64_t s[4]; } EzzRng;
void   ezz_rng_seed(EzzRng *rng, uint64_t seed);
double ezz_rng_uniform(EzzRng *rng, double lo, double hi);

#ifdef __cplusplus
}
#endif

#endif /* EISENSTEIN_VS_Z2_H */
