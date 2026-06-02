#define _POSIX_C_SOURCE 199309L
#include "eisenstein_vs_z2.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ── Snap error calculations ──────────────────────────────────────── */

double eisenstein_snap_error(double x, double y) {
    EisensteinInt snapped = snap_to_eisenstein(x, y);
    double sx, sy;
    eisenstein_to_cartesian(snapped.a, snapped.b, &sx, &sy);
    double dx = x - sx;
    double dy = y - sy;
    return sqrt(dx * dx + dy * dy);
}

double z2_snap_error(double x, double y) {
    Z2Int snapped = snap_to_z2(x, y);
    double dx = x - (double)snapped.x;
    double dy = y - (double)snapped.y;
    return sqrt(dx * dx + dy * dy);
}

/* ── Distance functions ──────────────────────────────────────────── */

double eisenstein_distance(EisensteinInt a, EisensteinInt b) {
    double ax, ay, bx, by;
    eisenstein_to_cartesian(a.a, a.b, &ax, &ay);
    eisenstein_to_cartesian(b.a, b.b, &bx, &by);
    double dx = ax - bx;
    double dy = ay - by;
    return sqrt(dx * dx + dy * dy);
}

double z2_distance(Z2Int a, Z2Int b) {
    double dx = (double)(a.x - b.x);
    double dy = (double)(a.y - b.y);
    return sqrt(dx * dx + dy * dy);
}

/* ── Internal: comparison for qsort ───────────────────────────────── */

static int double_cmp(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/* ── Internal: compute trial stats from errors ────────────────────── */

static TrialResult compute_stats(const double *errors, int n,
                                  const int64_t *snap_a, const int64_t *snap_b,
                                  int use_eisenstein) {
    TrialResult result;
    memset(&result, 0, sizeof(result));
    result.n = n;

    if (n <= 0) return result;

    /* Mean and max */
    double sum = 0.0, max_err = 0.0;
    for (int i = 0; i < n; i++) {
        sum += errors[i];
        if (errors[i] > max_err) max_err = errors[i];
    }
    result.mean_error = sum / n;
    result.max_error = max_err;

    /* Median */
    double *sorted = (double *)malloc((size_t)n * sizeof(double));
    memcpy(sorted, errors, (size_t)n * sizeof(double));
    qsort(sorted, (size_t)n, sizeof(double), double_cmp);
    result.median_error = (n % 2 == 0)
        ? (sorted[n/2 - 1] + sorted[n/2]) / 2.0
        : sorted[n/2];

    /* Standard deviation */
    double var_sum = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = errors[i] - result.mean_error;
        var_sum += diff * diff;
    }
    result.std_error = sqrt(var_sum / n);

    /* Recovery rates */
    int r001 = 0, r01 = 0, r05 = 0;
    for (int i = 0; i < n; i++) {
        if (errors[i] <= 0.01) r001++;
        if (errors[i] <= 0.1)  r01++;
        if (errors[i] <= 0.5)  r05++;
    }
    result.recovery_001 = (double)r001 / n;
    result.recovery_01 = (double)r01 / n;
    result.recovery_05 = (double)r05 / n;

    /* Packing uniqueness: count unique snap points */
    typedef struct { int64_t x; int64_t y; } Pair;
    Pair *pairs = (Pair *)malloc((size_t)n * sizeof(Pair));
    for (int i = 0; i < n; i++) {
        if (use_eisenstein) {
            pairs[i].x = snap_a[i];
            pairs[i].y = snap_b[i];
        } else {
            pairs[i].x = snap_a[i];
            pairs[i].y = snap_b[i];
        }
    }
    int unique = 0;
    for (int i = 0; i < n; i++) {
        int is_dup = 0;
        for (int j = 0; j < i; j++) {
            if (pairs[i].x == pairs[j].x && pairs[i].y == pairs[j].y) {
                is_dup = 1;
                break;
            }
        }
        if (!is_dup) unique++;
    }
    result.packing_unique = unique;
    result.packing_ratio = (double)unique / n;

    free(pairs);
    free(sorted);

    return result;
}

/* ── Trial runners ────────────────────────────────────────────────── */

static TrialResult run_trial_internal(
    const double *points_x, const double *points_y, int n,
    int use_eisenstein
) {
    TrialResult result;
    memset(&result, 0, sizeof(result));
    result.n = n;

    if (n <= 0 || !points_x || !points_y) return result;

    double *errors = (double *)malloc((size_t)n * sizeof(double));
    int64_t *snap_a = (int64_t *)malloc((size_t)n * sizeof(int64_t));
    int64_t *snap_b = (int64_t *)malloc((size_t)n * sizeof(int64_t));

    for (int i = 0; i < n; i++) {
        if (use_eisenstein) {
            EisensteinInt ei = snap_to_eisenstein(points_x[i], points_y[i]);
            double sx, sy;
            eisenstein_to_cartesian(ei.a, ei.b, &sx, &sy);
            double dx = points_x[i] - sx;
            double dy = points_y[i] - sy;
            errors[i] = sqrt(dx * dx + dy * dy);
            snap_a[i] = ei.a;
            snap_b[i] = ei.b;
        } else {
            Z2Int zi = snap_to_z2(points_x[i], points_y[i]);
            double dx = points_x[i] - (double)zi.x;
            double dy = points_y[i] - (double)zi.y;
            errors[i] = sqrt(dx * dx + dy * dy);
            snap_a[i] = zi.x;
            snap_b[i] = zi.y;
        }
    }

    result = compute_stats(errors, n, snap_a, snap_b, use_eisenstein);

    free(snap_b);
    free(snap_a);
    free(errors);

    return result;
}

TrialResult run_eisenstein_trial(const double *px, const double *py, int n) {
    return run_trial_internal(px, py, n, 1);
}

TrialResult run_z2_trial(const double *px, const double *py, int n) {
    return run_trial_internal(px, py, n, 0);
}

TrialResult run_timed_trial(const double *points_x, const double *points_y, int n, int use_eisenstein) {
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    TrialResult result = run_trial_internal(points_x, points_y, n, use_eisenstein);

    clock_gettime(CLOCK_MONOTONIC, &end);
    result.elapsed_s = (double)(end.tv_sec - start.tv_sec) +
                        (double)(end.tv_nsec - start.tv_nsec) / 1e9;
    return result;
}
