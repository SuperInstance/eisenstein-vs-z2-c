# eisenstein-vs-z2-c

C benchmark comparing Eisenstein (hexagonal A₂) lattice vs square (Z²) lattice for constraint quantization — snap error, packing density, and convergence rate.

## What This Gives You

- **Dual-lattice snap** — snap points to both Eisenstein and Z² lattices
- **Error comparison** — RMS and max snap error for each lattice
- **Packing analysis** — lattice packing density and covering radius
- **Convergence benchmarks** — how quickly each lattice converges under constraint dynamics
- **Zero dependencies** — pure C99, one header

## Quick Start

```c
#include "eisenstein_vs_z2.h"

/* Snap to Eisenstein lattice */
EisensteinInt e = snap_to_eisenstein(0.7, 0.3);
printf("Eisenstein snap: (%ld, %ld)\n", e.a, e.b);

/* Snap to Z² (square) lattice */
Z2Int z = snap_to_z2(0.7, 0.3);
printf("Z² snap: (%ld, %ld)\n", z.x, z.y);

/* Compare errors */
double e_err = eisenstein_snap_error(0.7, 0.3);
double z_err = z2_snap_error(0.7, 0.3);
printf("Error: Eisenstein=%.4f, Z²=%.4f\n", e_err, z_err);
```

## API Reference

| Function | Description |
|---|---|
| `snap_to_eisenstein(x, y)` | Nearest Eisenstein lattice point |
| `snap_to_z2(x, y)` | Nearest square lattice point |
| `eisenstein_norm(a, b)` | a² − ab + b² |
| `eisenstein_snap_error(x, y)` | Distance after Eisenstein snap |
| `z2_snap_error(x, y)` | Distance after Z² snap |

## Building

```bash
gcc -std=c99 -O2 -Wall bench_eisenstein_vs_z2.c src/eisenstein_vs_z2.c -lm -o bench
./bench
```

## How It Fits

C port of the lattice comparison benchmark:

- [eisenstein-vs-z2-rs](https://github.com/SuperInstance/eisenstein-vs-z2-rs) — Rust version
- [eisenstein-triples](https://github.com/SuperInstance/eisenstein-triples) — Eisenstein triple number theory
- [constraint-theory-core](https://github.com/SuperInstance/constraint-theory-core) — uses A₂ lattice for constraint quantization

## License

MIT
