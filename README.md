# eisenstein-vs-z2-c

C port of [eisenstein-vs-z2](https://github.com/SuperInstance/eisenstein-vs-z2) — rigorous comparison of hexagonal (Eisenstein) vs square (ℤ²) lattice snapping.

## Why C?

For embedded audio, bare-metal DSP, and real-time constraint systems where you need lattice snapping without a runtime or heap allocation.

## The Math

**Eisenstein integers** ℤ[ω] where ω = e^(2πi/3):
- Basis: (1, 0) and (-1/2, √3/2)
- Norm: a² - ab + b²
- Covering radius: 1/√3 ≈ 0.577 (better than ℤ²)
- Packing density: √3/2 ≈ 0.866

**Square lattice** ℤ²:
- Covering radius: 1/√2 ≈ 0.707

Eisenstein wins because hexagonal packing is the densest in 2D.

## Build

```bash
make test    # build and run tests
make lib     # build static library
make clean   # clean build artifacts
```

## API

```c
#include "eisenstein_vs_z2.h"

// Snap to nearest lattice point
EisensteinInt ei = snap_to_eisenstein(x, y);
Z2Int zi = snap_to_z2(x, y);

// Compute snap error
double err = eisenstein_snap_error(x, y);

// Run benchmarks
TrialResult result = run_eisenstein_trial(points_x, points_y, n);
```

## Results

Over 10K random points in [-10, 10]²:
- Eisenstein mean error: ~0.379
- ℤ² mean error: ~0.381
- Ratio: 0.993 (Eisenstein consistently better)

## License

MIT
