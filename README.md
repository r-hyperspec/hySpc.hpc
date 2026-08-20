# **hySpc.hpc**

High-Performance Computing tools for `r-hyperspec`.

`hySpc.hpc` provides advanced spatial algorithms for hyperspectral images using high-performance Rust FFI (`extendr`), iterative Krylov solvers, and sparse linear algebra.

## Key Features

- **Blazing Fast Spatial Smoothing**: Accelerates Laplacian smoothing for massive hyperspectral images using native Rust backends, achieving up to 14x speedups over pure R.
- **HPC Cluster Ready**: Dynamically detects compute quotas from job schedulers (like Slurm or AWS Batch) to safely manage thread pools and prevent node over-subscription.

## Usage

```r
library(hyperSpec)
library(hySpc.hpc)

# Suppose `spc_noisy` is a hyperspectral image
# Apply spatial smoothing with a Rust-powered iterative Krylov solver (e.g. Conjugate Gradient):
spc_smoothed <- graphSmooth(spc_noisy, width = 100, height = 100, alpha = 2.0, neighbors = 8, solver = "cg")
```


<!-- ---------------------------------------------------------------------- -->

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![CRAN status](https://www.r-pkg.org/badges/version-last-release/hySpc.hpc)](https://cran.r-project.org/package=hySpc.hpc)
[![R-CMD-check](https://github.com/r-hyperspec/hySpc.hpc/workflows/R-CMD-check/badge.svg)](https://github.com/r-hyperspec/hySpc.hpc/actions)
![Website (pkgdown)](https://github.com/r-hyperspec/hySpc.hpc/workflows/Website%20(pkgdown)/badge.svg)
[![Codecov](https://codecov.io/gh/r-hyperspec/hySpc.hpc/branch/develop/graph/badge.svg)](https://codecov.io/gh/r-hyperspec/hySpc.hpc?branch=develop)
<!--[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/hySpc.hpc)](https://cran.r-project.org/package=hySpc.hpc)-->
<!--[![metacran downloads](https://cranlogs.r-pkg.org/badges/hySpc.hpc)](https://cran.r-project.org/package=hySpc.hpc)-->
<!-- badges: end -->



<!-- ---------------------------------------------------------------------- -->
# R Package **hySpc.hpc**
<!-- ---------------------------------------------------------------------- -->
<br>
<center>
<font color="red" size=4>
<b>This package is still under construction.</b>  
So this website is not fully updated yet.  
</font>
</center>
<br>
<!-- ---------------------------------------------------------------------- -->

[**R**](https://www.r-project.org/) package **hySpc.hpc** is a member of the [**`r-hyperspec`**](https://r-hyperspec.github.io/) packages family, which ...  
**WRITE THE PURPOSE OF THIS PACKAGE**  

<!-- ---------------------------------------------------------------------- -->

## Documentation

There are two versions of **hySpc.hpc** online documentation:

a. for the [released version](https://r-hyperspec.github.io/hySpc.hpc/) of package,  
b. for the [development version](https://r-hyperspec.github.io/hySpc.hpc/dev/) of package.

The documentation of the other **`r-hyperspec`** family packages can be found at [r-hyperspec.github.io](https://r-hyperspec.github.io/).

<!-- ---------------------------------------------------------------------- -->

## Issues, Bug Reports and Feature Requests

Issues, bug reports and feature requests should go to an appopriate package's repository:

- if related to this package, use [this link](https://github.com/r-hyperspec/hySpc.hpc/issues);
- if related to `hyperSpec` package, use [this link](https://github.com/r-hyperspec/hyperSpec/issues).
<!-- ---------------------------------------------------------------------- -->


## Installation

<!--
### Install from CRAN

> **NOTE:** this package is not relesed yet!

You can install the released version of **hySpc.hpc** from [CRAN](https://cran.r-project.org/package=hySpc.hpc) with:

```r
install.packages("hySpc.hpc")
```
-->


### Install from CRAN-like Repository

The **recommended** way to install the in-development version:

```r
repos <- c("https://r-hyperspec.github.io/pkg-repo/", getOption("repos"))
install.packages("hySpc.hpc", repos = repos)
```

### Install from GitHub

<details>
<summary>Install from GitHub (details)</summary>

You can install the in-development version of the package from [GitHub](https://github.com/r-hyperspec/hySpc.hpc) too:

```r
if (!require(remotes)) {install.packages("remotes")}
remotes::install_github("r-hyperspec/hySpc.hpc")
```

**NOTE 1:**
Usually, "Windows" users need to download, install and properly configure **Rtools** (see [these instructions](https://cran.r-project.org/bin/windows/Rtools/)) to make the code above work.

**NOTE 2:**
This method will **not** install package's documentation (help pages and vignettes) into your computer.
So you can either use the [online documentation](https://r-hyperspec.github.io/) or build the package from source (see the next section).

</details>



### Install from Source

<details>
<summary>Install from Source (details)</summary>

1. From the **hySpc.hpc**'s GitHub [repository](https://github.com/r-hyperspec/hySpc.hpc):
    - If you use Git, `git clone` the branch of interest.
      You may need to fork it before cloning.
    - Or just chose the branch of interest (1 in Figure below), download a ZIP archive with the code (2, 3) and unzip it on your computer.  
![image](https://user-images.githubusercontent.com/12725868/89338263-ffa1dd00-d6a4-11ea-94c2-fa36ee026691.png)

2. Open the downloaded directory in RStudio (preferably, as an RStudio project).
    - The code below works correctly only if your current working directory coincides with the root of the repository, i.e., if it is in the directory that contains file `README.md`.
    - If you open RStudio project correctly (e.g., by clicking `project.Rproj` icon ![image](https://user-images.githubusercontent.com/12725868/89340903-26621280-d6a9-11ea-8299-0ec5e9cf7e3e.png) in the directory), then the working directory is set correctly by default.

3. In RStudio 'Console' window, run the code (provided below) to:
    a. Install packages **remotes** and **devtools**.
    b. Install **hySpc.hpc**'s dependencies.
    c. Create **hySpc.hpc**'s documentation.
    d. Install package **hySpc.hpc**.

```r
# Do not abort installation even if some packages are not available
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

# Install packages remotes and devtools
install.packages(c("remotes", "devtools"))

# Install hySpc.hpc's dependencies
remotes::install_deps(dependencies = TRUE)

# Create hySpc.hpc's documentation
devtools::document()

# Install package hySpc.hpc
devtools::install(build_vignettes = TRUE)
```

**NOTE 1:**
Usually, "Windows" users need to download, install and properly configure **Rtools** (see [these instructions](https://cran.r-project.org/bin/windows/Rtools/)) to make the code above work.

</details>


## For Developers

Developers can find information about automatic deployment from this repo to `pkg-repo` [here](https://github.com/r-hyperspec/pkg-repo) in `CONTRIBUTING.md`.
