#ifndef EISENSTEIN_VS_Z2_H
#define EISENSTEIN_VS_Z2_H

/* ===================================================================
 * Eisenstein integer lattice and Z2 (square lattice) operations
 * =================================================================== */

/* Eisenstein integer: a + bω where ω = e^(2πi/3) = -1/2 + i√3/2 */

typedef struct {
    long long a;  /* real Eisenstein component */
    long long b;  /* imaginary Eisenstein component */
} eis_int_t;

typedef struct {
    double x;
    double y;
} eis_point_t;

/* --- Eisenstein operations --- */

/* Convert Eisenstein (a,b) to Cartesian (x,y) */
eis_point_t eis_to_cartesian(long long a, long long b);

/* Convert Cartesian (x,y) to continuous Eisenstein (a,b) */
eis_point_t eis_from_cartesian(double x, double y);

/* Snap a 2D point to the nearest Eisenstein integer lattice point */
eis_int_t eis_snap(double x, double y);

/* Snap a 2D point to the nearest Z2 (square) lattice point */
eis_int_t eis_snap_z2(double x, double y);

/* Eisenstein norm: a² - ab + b² */
long long eis_norm(long long a, long long b);

/* Z2 norm: a² + b² */
long long eis_z2_norm(long long a, long long b);

/* Euclidean distance from original point to snapped Eisenstein point */
double eis_snap_error(double x, double y);

/* Euclidean distance from original point to snapped Z2 point */
double eis_z2_snap_error(double x, double y);

/* Distance between two Cartesian points */
double eis_distance(eis_point_t a, eis_point_t b);

/* Theoretical covering radius */
double eis_covering_radius_eisenstein(void);  /* 1/√3 ≈ 0.5774 */
double eis_covering_radius_z2(void);          /* 1/√2 ≈ 0.7071 */

/* Voronoi cell area */
double eis_voronoi_area_eisenstein(void);     /* √3/2 ≈ 0.8660 */
double eis_voronoi_area_z2(void);             /* 1.0 */

/* Batch snap: snap n points, compute mean/min/max errors */
typedef struct {
    double mean_error;
    double min_error;
    double max_error;
    int n;
} eis_bench_result_t;

eis_bench_result_t eis_benchmark_eisenstein(const double *xs, const double *ys, int n);
eis_bench_result_t eis_benchmark_z2(const double *xs, const double *ys, int n);

#endif /* EISENSTEIN_VS_Z2_H */
