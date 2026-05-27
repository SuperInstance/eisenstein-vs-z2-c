/*
 * eisenstein_vs_z2 — lattice snapping and comparison
 *
 * C99 port of the Python benchmark.
 * Eisenstein integers (hexagonal lattice) vs Z2 (square lattice).
 */

#include "eisenstein_vs_z2.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ─── Constants ────────────────────────────────────────────── */

#define SQRT3_2 0.8660254037844386  /* sqrt(3)/2 */
#define PI      3.14159265358979323846

/* ─── Coordinate conversion ────────────────────────────────── */

EzzPoint ezz_eisenstein_to_cartesian(long a, long b)
{
    EzzPoint p;
    p.x = (double)a - (double)b * 0.5;
    p.y = (double)b * SQRT3_2;
    return p;
}

EzzPoint ezz_cartesian_to_eisenstein(double x, double y)
{
    EzzPoint p;
    double b = y / SQRT3_2;
    double a = x + b * 0.5;
    p.x = a;
    p.y = b;
    return p;
}

/* ─── Snapping ─────────────────────────────────────────────── */

double ezz_snap_eisenstein(double x, double y,
                           EzzPoint *out_cart, EzzEisenstein *out_eis)
{
    EzzPoint cont = ezz_cartesian_to_eisenstein(x, y);
    long a = (long)round(cont.x);
    long b = (long)round(cont.y);
    if (out_eis) { out_eis->a = a; out_eis->b = b; }
    EzzPoint snapped = ezz_eisenstein_to_cartesian(a, b);
    if (out_cart) { *out_cart = snapped; }
    double dx = x - snapped.x;
    double dy = y - snapped.y;
    return sqrt(dx * dx + dy * dy);
}

double ezz_snap_z2(double x, double y,
                   EzzPoint *out_cart, EzzEisenstein *out_z2)
{
    long a = (long)round(x);
    long b = (long)round(y);
    if (out_z2) { out_z2->a = a; out_z2->b = b; }
    if (out_cart) { out_cart->x = (double)a; out_cart->y = (double)b; }
    double dx = x - (double)a;
    double dy = y - (double)b;
    return sqrt(dx * dx + dy * dy);
}

/* ─── Batch snapping ──────────────────────────────────────── */

void ezz_snap_eisenstein_batch(const double *xy, size_t n,
                               double *out_cart, long *out_lat,
                               double *errors)
{
    for (size_t i = 0; i < n; i++) {
        EzzPoint sc;
        EzzEisenstein ei;
        double err = ezz_snap_eisenstein(xy[2*i], xy[2*i+1], &sc, &ei);
        if (out_cart) { out_cart[2*i] = sc.x; out_cart[2*i+1] = sc.y; }
        if (out_lat)  { out_lat[2*i] = ei.a; out_lat[2*i+1] = ei.b; }
        if (errors)   { errors[i] = err; }
    }
}

void ezz_snap_z2_batch(const double *xy, size_t n,
                       double *out_cart, long *out_lat,
                       double *errors)
{
    for (size_t i = 0; i < n; i++) {
        EzzPoint sc;
        EzzEisenstein z2;
        double err = ezz_snap_z2(xy[2*i], xy[2*i+1], &sc, &z2);
        if (out_cart) { out_cart[2*i] = sc.x; out_cart[2*i+1] = sc.y; }
        if (out_lat)  { out_lat[2*i] = z2.a; out_lat[2*i+1] = z2.b; }
        if (errors)   { errors[i] = err; }
    }
}

/* ─── Norms and metrics ────────────────────────────────────── */

long ezz_eisenstein_norm(long a, long b)
{
    return a * a - a * b + b * b;
}

double ezz_max_snap_error_eisenstein(void) { return 1.0 / sqrt(3.0); }
double ezz_max_snap_error_z2(void)         { return 1.0 / sqrt(2.0); }
double ezz_voronoi_cell_area_eisenstein(void) { return SQRT3_2; }
double ezz_voronoi_cell_area_z2(void)         { return 1.0; }

/* ─── PRNG (xoshiro256**) ──────────────────────────────────── */

static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

void ezz_rng_seed(EzzRng *rng, uint64_t seed)
{
    /* SplitMix64 to initialize state */
    uint64_t z;
    for (int i = 0; i < 4; i++) {
        seed += 0x9e3779b97f4a7c15ULL;
        z = seed;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        rng->s[i] = z ^ (z >> 31);
    }
}

static uint64_t rng_next(EzzRng *rng)
{
    uint64_t *s = rng->s;
    uint64_t result = rotl64(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl64(s[3], 45);
    return result;
}

double ezz_rng_uniform(EzzRng *rng, double lo, double hi)
{
    double u = (rng_next(rng) >> 11) * 0x1.0p-53;  /* [0, 1) */
    return lo + u * (hi - lo);
}

/* ─── Comparison for qsort ─────────────────────────────────── */

static int cmp_double(const void *a, const void *b)
{
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

/* ─── Benchmark trial ──────────────────────────────────────── */

static EzzTrialResult run_trial(const double *points_xy, size_t n,
                                int is_eisenstein)
{
    EzzTrialResult r;
    memset(&r, 0, sizeof(r));

    double *errors = (double *)malloc(n * sizeof(double));
    if (!errors) return r;

    /* Snap and collect errors */
    if (is_eisenstein) {
        for (size_t i = 0; i < n; i++) {
            errors[i] = ezz_snap_eisenstein(points_xy[2*i], points_xy[2*i+1],
                                            NULL, NULL);
        }
    } else {
        for (size_t i = 0; i < n; i++) {
            errors[i] = ezz_snap_z2(points_xy[2*i], points_xy[2*i+1],
                                    NULL, NULL);
        }
    }

    /* Sort for percentile calculations */
    qsort(errors, n, sizeof(double), cmp_double);

    /* Compute statistics */
    double sum = 0, sum2 = 0;
    for (size_t i = 0; i < n; i++) { sum += errors[i]; sum2 += errors[i] * errors[i]; }
    r.mean_error = sum / (double)n;
    r.std_error  = sqrt(sum2 / (double)n - r.mean_error * r.mean_error);
    r.max_error  = errors[n - 1];
    r.median_error = (n % 2 == 0)
        ? (errors[n/2 - 1] + errors[n/2]) * 0.5
        : errors[n/2];
    r.p25 = errors[(size_t)(n * 0.25)];
    r.p75 = errors[(size_t)(n * 0.75)];
    r.p95 = errors[(size_t)(n * 0.95)];
    r.p99 = errors[(size_t)(n * 0.99)];

    /* Recovery rates */
    size_t cnt001 = 0, cnt01 = 0, cnt05 = 0;
    for (size_t i = 0; i < n; i++) {
        if (errors[i] <= 0.01) cnt001++;
        if (errors[i] <= 0.1)  cnt01++;
        if (errors[i] <= 0.5)  cnt05++;
    }
    r.recovery_001 = (double)cnt001 / (double)n;
    r.recovery_01  = (double)cnt01  / (double)n;
    r.recovery_05  = (double)cnt05  / (double)n;

    free(errors);
    return r;
}

EzzTrialResult ezz_run_trial_eisenstein(size_t n_points, uint64_t seed)
{
    EzzRng rng;
    ezz_rng_seed(&rng, seed);
    double *pts = (double *)malloc(n_points * 2 * sizeof(double));
    for (size_t i = 0; i < n_points; i++) {
        pts[2*i]   = ezz_rng_uniform(&rng, -10.0, 10.0);
        pts[2*i+1] = ezz_rng_uniform(&rng, -10.0, 10.0);
    }
    EzzTrialResult r = run_trial(pts, n_points, 1);
    free(pts);
    return r;
}

EzzTrialResult ezz_run_trial_z2(size_t n_points, uint64_t seed)
{
    EzzRng rng;
    ezz_rng_seed(&rng, seed);
    double *pts = (double *)malloc(n_points * 2 * sizeof(double));
    for (size_t i = 0; i < n_points; i++) {
        pts[2*i]   = ezz_rng_uniform(&rng, -10.0, 10.0);
        pts[2*i+1] = ezz_rng_uniform(&rng, -10.0, 10.0);
    }
    EzzTrialResult r = run_trial(pts, n_points, 0);
    free(pts);
    return r;
}
