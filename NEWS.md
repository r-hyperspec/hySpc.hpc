# hySpc.hpc 0.0.0.9000 (Development Version)

### Core Features & Native Kernels
- **Spatial Graph Smoothing**: Implemented `graphSmooth()` for 2D/3D spatial smoothing on `hyperSpec` datacubes using 4-neighbor and 8-neighbor grid graph Laplacians.
- **Native Rust Solvers**: Implemented zero-allocation `csc_matvec` Krylov Conjugate Gradient (`solve_cg`) and BiCGSTAB (`solve_bicgstab`) solvers in Rust via `extendr`.
- **Direct Sparse Cholesky**: Integrated `faer` and `dgcmatrix-faer-bridge` for direct sparse factorizations.
- **HPC Thread Management**: Added `set_hpc_threads()` and `get_hpc_thread_limit()` with auto-detection for Slurm, AWS Batch, and cgroup CPU quotas.

### Documentation & Site Enhancements
- Added **Get Started Guide** (`articles/get_started.html`).
- Added **Spatial Smoothing Vignette** (`articles/spatial_smoothing.html`).
- Added **Benchmarking Vignette** (`articles/benchmarking.html`).
- Added **Distributed Spectroscopy & Global HPC Clusters Guide** (`articles/hpc_cluster.html`).
- Re-structured `_pkgdown.yml` site navigation bar, reference indexing, and logical article ordering.
