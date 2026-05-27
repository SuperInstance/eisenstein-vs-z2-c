# eisenstein-vs-z2-c

C99 port of the [Eisenstein vs ℤ²](https://github.com/SuperInstance/eisenstein-vs-z2) benchmark — comparing hexagonal (Eisenstein integer) and square (ℤ²) lattice snapping.

Designed for embedded constraint-solving and lattice geometry applications.

## What It Does

Eisenstein integers form the hexagonal lattice (proven densest packing in 2D via Thue's theorem). This library provides:

- **Coordinate conversion** between Cartesian and Eisenstein coordinates
- **Lattice snapping** — snap arbitrary 2D points to nearest Eisenstein or ℤ² lattice points
- **Batch operations** for high-throughput snapping
- **Benchmark utilities** for comparing the two lattices

The hexagonal lattice has a smaller Voronoi cell and lower covering radius than the square lattice, yielding ~13–17% lower snap error on average.

## C API

```c
#include <eisenstein_vs_z2.h>

/* Snap a point to the Eisenstein lattice */
EzzPoint snapped;
EzzEisenstein coords;
double error = ezz_snap_eisenstein(1.5, 2.3, &snapped, &coords);

/* Snap to Z2 for comparison */
double error_z2 = ezz_snap_z2(1.5, 2.3, &snapped, &(EzzEisenstein){0});

/* Batch snap */
double points[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
double out_cart[6];
long out_lat[6];
double errors[3];
ezz_snap_eisenstein_batch(points, 3, out_cart, out_lat, errors);

/* Run benchmark */
EzzTrialResult r = ezz_run_trial_eisenstein(10000, 42);
printf("Mean error: %.6f\n", r.mean_error);
```

## Build

```sh
make        # builds libeisenstein_vs_z2.a
make test   # builds and runs tests
make clean
```

No dependencies beyond C99, `libc`, and `libm`.

## Port Notes

Ported from the Python benchmark which uses NumPy for vectorized operations. The C version:

- Replaces NumPy with plain C loops (cache-friendly, no allocation in hot path)
- Uses xoshiro256\*\* PRNG instead of NumPy's `default_rng`
- Implements percentile calculation via sorting instead of NumPy's `percentile`
- All math is pure `double` precision

## License

MIT
