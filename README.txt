=============================================================================================
TTI-FDModeling Project README
=============================================================================================

1. Overview

This project is a C-language finite-difference modeling codebase for 2D TTI elastic wave
simulation. It is built around the Seismic Unix (SU / CWP) ecosystem and uses `getpar`
style command-line parameters together with SU traces and raw binary model files.

The repository currently contains:
- full source code under `src/`
- installed headers under `include/`
- built static libraries under `lib/`
- built executables under `bin/`
- demo scripts, sample data, generated demo outputs, and plotting scripts under `demo/`

Copyright (c) 2026 WaveTomo. All rights reserved.


2. Current Project Status

2.1 Maintained main source programs

The current source tree under `src/main/suTTI-FDModeling2dTTIelastic/` contains 3 maintained main
program variants:

- `suTTI-FDModeling2dTTIelastic_coordinate_rotation.c`
  2D TTI elastic modeling with the rotated-operator / coordinate-rotation formulation
  Source `sdoc` still labels this variant as `SUTTI-FDModeling2DTTIELASTIC.V1.0`

- `suTTI-FDModeling2dTTIelastic_tensor_rotation.c`
  2D TTI elastic modeling with the rotated-tensor formulation
  Source `sdoc` still labels this variant as `SUTTI-FDModeling2DTTIELASTIC.V3.0`

- `suTTI-FDModeling2dTTIelastic_bond_transformation.c`
  2D TTI elastic modeling with the Bond-transformation formulation
  Source `sdoc` still labels this variant as `SUTTI-FDModeling2DTTIELASTIC.V4.0`

Important: the previously documented self-joint `v2.0` source file is no longer present in
the current source tree. However, related legacy build artifacts still exist in `bin/` and
`src/main/suTTI-FDModeling2dTTIelastic/build/`.

2.2 Main functionality

The codebase currently supports:
- reading `vp`, `epsilon`, `delta`, and `theta` model files
- optional `vs` and `rho` input files
- reading source and receiver geometry from SU-format files
- forward modeling on a rotated staggered grid for 2D TTI elastic media
- effective absorbing-layer handling and optional free-surface handling
- saving snapshots and writing seismic traces
- helper utilities for source generation, receiver generation, and binary model cropping

2.3 Common main-program parameters

According to the current main-program source files, the key parameters include:
- grid and time: `nz nx nt dx dz dt`
- model files: `vpfile epsilonfile deltafile thetafile`
- optional model files: `vsfile rhofile`
- source file: `sourcefile`
- receiver file: `recposfile`
- source controls: `stype fpeak scalar`
- recording controls: `seistype mt`
- boundary controls: `ealThick eal_alpha freesurface`
- output controls: `savew wsave snap_sample snap_prefix trace_prefix`
- other controls: `fx fz mode shotID energy`


3. Directory Layout

3.1 Root directories

- `src/`
  full source tree
- `include/`
  installed public headers
- `lib/`
  built static libraries
- `bin/`
  built executables
- `demo/`
  sample models, demo scripts, generated outputs, and plotting scripts

3.2 Source layout

- `src/Makefile`
  top-level build entry
- `src/include/`
  original header source location
- `src/lib/`
  library source code
- `src/main/`
  main programs and utility tools

3.3 Library modules

- `src/lib/TTI-FDModeling/`
  core numerical modeling code
- `src/lib/absorbc/`
  absorbing-boundary related code
- `src/lib/rockmod/`
  model-parameter and rock-property conversion
- `src/lib/supro/`
  SU and binary I/O helpers
- `src/lib/wavelet/`
  wavelet generation


4. Library Responsibilities

`TTI-FDModeling`
- finite-difference derivatives, gradient, divergence, and Laplacian
- source and receiver structures
- source and receiver file reading
- staggered-grid parameter transforms
- stability checks
- snapshot output and trace output
- isotropic and anisotropic elastic parameter conversion

`absorbc`
- model padding and extension
- MEAL absorbing-boundary coefficient generation
- Higdon boundary handling
- `eal.c` and `pml.c` currently remain placeholder or minimal files

`rockmod`
- `vp -> vs`, `vp -> rho`, and `vp -> slowness` conversion
- interface gridding helpers

`supro`
- binary file read/write helpers
- SU header scaling helpers
- SU trace and shot counting helpers
- snapshot writing helpers

`wavelet`
- Ricker wavelet
- derivative Ricker wavelet
- Gaussian wavelet
- derivative Gaussian wavelet


5. Utility Programs

The utility programs under `src/main/wavelet/` are unchanged and remain:

- `sumksrc2d.c`
  generates an SU-format source file
- `sumkrcev2d.c`
  generates an SU-format receiver-position file
- `cutbin2d.c`
  crops a 2D subregion from a larger binary model


6. Build Information

6.1 Compiler and flags

The top-level build file is `src/Makefile`.

Default compiler:
- `icc`

Default compile flags:
- `-Wall -pedantic -std=c11 -qopenmp -O3 -g`

6.2 External dependencies

The current Makefiles reference:
- `CWPROOT` include and library directories
- `lsu`
- `lpar`
- `lcwp`
- `lm`

This project is intended for a Unix-like build environment with Seismic Unix configured
before compilation.

6.3 Install targets

Build outputs are installed into the project root:
- `include/`
- `lib/`
- `bin/`


7. Current Demo Layout

7.1 Demo scripts

The current demo directory contains 3 shell scripts:

- `demo/test_bp01.sh`
  runs `suTTI-FDModeling2dTTIelastic_coordinate_rotation`
- `demo/test_bp02.sh`
  runs `suTTI-FDModeling2dTTIelastic_tensor_rotation`
- `demo/test_bp03.sh`
  runs `suTTI-FDModeling2dTTIelastic_bond_transformation`

Compared with the previous project state, `test_bp04.sh` is no longer present.

7.2 Demo input model data

The sample input models remain in `demo/file/`:
- `vp_bp_1126_7872_10m.bin`
- `delta_bp_1126_7872_10m.bin`
- `epsilon_bp_1126_7872_10m.bin`
- `theta_bp_1126_7872_10m.bin`

7.3 Demo runtime parameters

All three demo scripts currently use:
- `nz=1126`
- `nx=1251`
- `nt=8001`
- `dx=10`
- `dz=10`
- `dt=0.001`
- `stype=1`
- `fpeak=10`
- `sdepth=5`
- `sx=625`
- `gdel=5`
- `offsetmin=0`
- `dtr=1`
- `offsetmax=2000`
- `seistype=0`
- `mt=2`
- `freesurface=0`
- `savew=0`
- `wsave=2`
- `snap_sample=1000`

The current demos are configured to save pressure snapshots rather than recorded traces.

7.4 Generated demo outputs currently present

The repository now includes generated outputs from previous demo runs:

- `demo/test_bp01/`
  cropped `vp`, `delta`, `epsilon`, `theta` models and pressure snapshots
- `demo/test_bp02/`
  pressure snapshots
- `demo/test_bp03/`
  pressure snapshots

These are generated artifacts, not source files.


8. Plotting and Post-processing

The repository now includes a plotting directory:

- `demo/plotFigure/plotBP.m`
  MATLAB plotting script for model panels and snapshot comparison
- `demo/plotFigure/lib/read2bin.m`
- `demo/plotFigure/lib/get_precision_bytes.m`
- `demo/plotFigure/lib/fnmemb.m`

Current generated figure files include:
- `demo/plotFigure/Fig10.pdf`
- `demo/plotFigure/Fig11.pdf`

From `plotBP.m`, the plotting workflow compares snapshots from:
- `test_bp01` coordinate rotation
- `test_bp02` tensor rotation
- `test_bp03` bond transformation


9. Current Executables and Artifacts

9.1 Built executables in `bin/`

Currently present executables include:
- `cutbin2d`
- `sumksrc2d`
- `sumkrcev2d`
- `suTTI-FDModeling2dTTIelastic_coordinate_rotation`
- `suTTI-FDModeling2dTTIelastic_tensor_rotation`
- `suTTI-FDModeling2dTTIelastic_bond_transformation`
- `suTTI-FDModeling2dTTIelastic.v1.0`
- `suTTI-FDModeling2dTTIelastic.v2.0`
- `suTTI-FDModeling2dTTIelastic.v3.0`
- `suTTI-FDModeling2dTTIelastic.v4.0`

The last four are legacy version-named executables that still remain in the repository as
built artifacts.

9.2 Built static libraries in `lib/`

Currently present static libraries:
- `libabsorbc.a`
- `libTTI-FDModeling.a`
- `librockmod.a`
- `libsupro.a`
- `libwavelet.a`

9.3 Legacy object artifacts

The directory `src/main/suTTI-FDModeling2dTTIelastic/build/` still contains legacy object files:
- `suTTI-FDModeling2dTTIelastic.v1.0.o`
- `suTTI-FDModeling2dTTIelastic.v2.0.o`
- `suTTI-FDModeling2dTTIelastic.v3.0.o`
- `suTTI-FDModeling2dTTIelastic.v4.0.o`

It also contains object files for the current descriptive source names:
- `suTTI-FDModeling2dTTIelastic_coordinate_rotation.o`
- `suTTI-FDModeling2dTTIelastic_tensor_rotation.o`
- `suTTI-FDModeling2dTTIelastic_bond_transformation.o`


10. Current File Snapshot

At the time of this README update, the repository contains:
- 126 total files
- 31 C source files
- 10 header files
- 11 Makefiles
- 3 shell scripts
- 4 MATLAB `.m` files
- 2 PDF files
- 5 static libraries
- 35 object files
- 14 `.bin` files including both sample inputs and generated outputs


11. Notes for Future Maintenance

The current repository state mixes source code with generated outputs and legacy build
artifacts. The main points to keep in mind are:

- the maintained source variants are now 3 descriptive programs, not 4 source files
- the legacy `v1.0` to `v4.0` executable names still remain in `bin/`
- `v2.0` is currently present as an artifact but not as an active source file
- demo output directories and figure PDFs are now part of the repository
- `demo/plotFigure/` should be documented as the current post-processing workflow

If the repository is intended for cleaner release or publication, it would be better to
separate source code from generated demo outputs and stale build products.
