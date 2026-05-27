#include "eisenstein_vs_z2.h"
#include <stdlib.h>
#include <string.h>

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

/* ── Internal: comparison for qsort ───────────────────────────────── */

static int double_cmp(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/* ── Trial runners ────────────────────────────────────────────────── */

static TrialResult run_trial_internal(
    const double *points_x, const double *points_y, int n,
    int use_eisenstein
) {
    TrialResult result = {0};
    result.n = n;

    double *errors = (double *)malloc(n * sizeof(double));
    double *snap_x = (double *)malloc(n * sizeof(double));
    double *snap_y = (double *)malloc(n * sizeof(double));

    /* Snap all points */
    for (int i = 0; i < n; i++) {
        if (use_eisenstein) {
            EisensteinInt ei = snap_to_eisenstein(points_x[i], points_y[i]);
            eisenstein_to_cartesian(ei.a, ei.b, &snap_x[i], &snap_y[i]);
        } else {
            Z2Int zi = snap_to_z2(points_x[i], points_y[i]);
            snap_x[i] = (double)zi.x;
            snap_y[i] = (double)zi.y;
        }
        double dx = points_x[i] - snap_x[i];
        double dy = points_y[i] - snap_y[i];
        errors[i] = sqrt(dx * dx + dy * dy);
    }

    /* Compute stats */
    double sum = 0.0, max_err = 0.0;
    for (int i = 0; i < n; i++) {
        sum += errors[i];
        if (errors[i] > max_err) max_err = errors[i];
    }
    result.mean_error = sum / n;
    result.max_error = max_err;

    /* Median: sort and pick middle */
    double *sorted = (double *)malloc(n * sizeof(double));
    memcpy(sorted, errors, n * sizeof(double));
    qsort(sorted, n, sizeof(double), double_cmp);
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

    /* Packing uniqueness */
    /* Simple hash set via sorted pairs */
    typedef struct { int64_t x; int64_t y; } Pair;
    Pair *pairs = (Pair *)malloc(n * sizeof(Pair));
    for (int i = 0; i < n; i++) {
        if (use_eisenstein) {
            EisensteinInt ei = snap_to_eisenstein(points_x[i], points_y[i]);
            pairs[i].x = ei.a;
            pairs[i].y = ei.b;
        } else {
            Z2Int zi = snap_to_z2(points_x[i], points_y[i]);
            pairs[i].x = zi.x;
            pairs[i].y = zi.y;
        }
    }
    /* Count unique by naive O(n²) — fine for benchmark sizes */
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
    free(snap_y);
    free(snap_x);
    free(errors);

    return result;
}

TrialResult run_eisenstein_trial(const double *px, const double *py, int n) {
    return run_trial_internal(px, py, n, 1);
}

TrialResult run_z2_trial(const double *px, const double *py, int n) {
    return run_trial_internal(px, py, n, 0);
}
