# LDPM + MDOM compression-prism template

This external example runs the repository's existing 100 x 100 x 200 mm
LDPM compression prism with Chrono's MPI multidomain (MDOM) module. It uses
the mesh and material inputs in
`template_project_ldpm/Unconfined Compression/LDPMgeo000Box000`, so the large
mesh files are not duplicated.

The MPI launch always uses one coordinator rank in addition to the requested
worker domains. For example, two worker domains require three MPI ranks.

## Requirements

- A C++17 compiler and CMake 3.18 or newer
- An MPI implementation with C++ support (Open MPI or MPICH)
- The normal Chrono dependencies required by this repository

## 1. Build Chrono with LDPM and MDOM

From the repository root:

```bash
cmake -S . -B build-mdom \
  -DCMAKE_BUILD_TYPE=Release \
  -DCH_ENABLE_MODULE_LDPM=ON \
  -DCH_ENABLE_MODULE_MULTIDOMAIN=ON
cmake --build build-mdom -j
```

Check the CMake summary before building. It should show both the LDPM and
multidomain modules enabled and should report a detected MPI C++ installation.

## 2. Build the external example

```bash
cmake -S benchmark/mdom_compression -B build-mdom-compression \
  -DCMAKE_BUILD_TYPE=Release \
  -DChrono_DIR="$PWD/build-mdom/cmake"
cmake --build build-mdom-compression -j
```

Keeping this as an external CMake project demonstrates the same way another
application can consume installed or locally built Chrono libraries with
`find_package(Chrono REQUIRED COMPONENTS LDPM MULTIDOMAIN CONFIG)`.

## 3. Run the prism

Run a short verification with two worker domains plus one coordinator:

```bash
benchmark/mdom_compression/run_prism.sh 2 quick
```

Run the complete compression history with the same loading used for the
validated prism campaign:

```bash
benchmark/mdom_compression/run_prism.sh 2 full
```

The full case uses a 5e-7 s step through 0.4 s (800,000 steps), so start with
`quick`. Output is written under `build-mdom-compression`, including `hist.dat`
and one per-rank output directory. The template disables VTK output by default
to keep the first run small; remove `--vtk_stride=0` in `run_prism.sh` if you
want snapshots.

If your build directories are elsewhere, set these environment variables:

```bash
MDOM_CHRONO_BUILD=/path/to/chrono-build \
MDOM_CASE_BUILD=/path/to/example-build \
benchmark/mdom_compression/run_prism.sh 2 quick
```

You can also launch the executable directly. The number passed to MPI is
`worker domains + 1`:

```bash
mpiexec -n 3 build-mdom-compression/ldpm_uc \
  --dt=5e-7 --t_end=0.0001 --slice_axis=x --partition_update \
  --no_diagonal_lumping --no_tag_reaction_selector \
  --fix_top_rotations --fix_boundary --vtk_stride=0
```

## Expected quick-run checks

A healthy run exits with status zero and prints a `RUN_SUMMARY` line. The
example also audits the global node and element ownership, the top boundary,
and specimen mass after partitioning. MPI errors, missing shared entities, or
non-finite response values should be treated as a failed run.

## Adapting the template

- Change `LDPM_EXAMPLE_DIR` in `CMakeLists.txt` to point at another LDPM input
  directory.
- Adjust the specimen bounds and top/bottom planes in `LDPM_UC.cpp` when the
  geometry dimensions change.
- Use `--slice_axis=x`, `y`, or `z` for sliced decomposition.
- Use `--grid=2x2x2` with eight worker domains (nine MPI ranks total) for a
  Cartesian decomposition.
- Add `--partition_update_each_step` only when testing dynamic repartitioning;
  it adds substantial communication cost.

The default direct boundary prescription is intentional. Migrating the serial
example's plate/tie links across MPI domains is not supported by this template.
