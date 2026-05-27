# eisenstein-vs-z2-c

C port of [eisenstein-vs-z2](https://github.com/SuperInstance/eisenstein-vs-z2) — Eisenstein integer lattice vs Z² (square lattice) comparison for 2D constraint snapping.

## What It Does

- **Eisenstein integer arithmetic**: basis conversion (Eisenstein ↔ Cartesian), norm calculation (a²−ab+b²)
- **Lattice snapping**: snap arbitrary 2D points to nearest Eisenstein or Z² lattice point
- **Error analysis**: compute Euclidean snap error for each lattice type
- **Theoretical properties**: covering radius, Voronoi cell area
- **Batch benchmarking**: mean/min/max error over point sets

The Eisenstein (hexagonal) lattice has provably lower snap error than Z² due to Thue's theorem (optimal 2D packing).

## Building

```bash
make lib    # static library
make test   # build and run tests
```

## License

MIT
