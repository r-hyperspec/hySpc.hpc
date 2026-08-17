#' Spatial Graph Smoothing
#'
#' Applies a graph-Laplacian spatial smoother to a [hyperSpec::hyperSpec]
#' object whose spectra are laid out on a regular `width` x `height`
#' image grid. For each wavelength band the smoothed image is the
#' solution of
#'
#' \deqn{(I + \alpha L)\, x = b,}
#'
#' where `L = D - W` is the combinatorial Laplacian of the pixel
#' adjacency graph (4- or 8-connectivity) and `b` is the original band.
#'
#' Two backends are available:
#'
#' * `backend = "rust"` (default): the high-performance Rust kernel
#'   ([graph_smooth_rust()]). The pixel neighborhood graph is rebuilt in
#'   Rust with `petgraph` from `width`/`height`/`neighbors` (no adjacency
#'   matrix crosses the FFI boundary), the Laplacian `L = D - W` is
#'   assembled in sparse (CSC) form, and each band is solved with an
#'   iterative Krylov method (`solver`: CG or BiCGSTAB). Requires the
#'   package's Rust extension to be compiled. Silently falls back to "r"
#'   if unavailable or fails.
#' * `backend = "r"`: the pure-R baseline ([graph_smooth_r()]) built on
#'   [Matrix] sparse routines. Here the Laplacian is assembled as a
#'   `dgCMatrix` and solved with [Matrix::solve()]. Used as the
#'   correctness reference and the pure-R side of the benchmark.
#' @param x a [hyperSpec::hyperSpec] object whose `spc` slot has
#'   `nrow(x) == width * height`.
#' @param width,height integer image dimensions.
#' @param alpha non-negative numeric smoothing strength. Larger values
#'   produce smoother output; `alpha = 0` returns `x` unchanged.
#' @param neighbors integer, either `4` (von Neumann) or `8` (Moore)
#'   neighborhood connectivity.
#' @param backend one of `"r"` or `"rust"`; selects the solver
#'   implementation.
#' @param solver iterative Krylov method used by the Rust backend, either
#'   `"cg"` (Conjugate Gradient, the default; `I + alpha L` is symmetric
#'   positive-definite) or `"bicgstab"` (BiCGSTAB, for the general case).
#'   Ignored by the R backend, which uses a direct [Matrix::solve()].
#'
#' @return A [hyperSpec::hyperSpec] object with the same metadata as
#'   `x` whose `spc` slot has been replaced by the smoothed spectra.
#'
#' @seealso [graph_smooth_r()], [graph_smooth_rust()].
#'
#' @export
setGeneric(
  "graphSmooth",
  function(x, width, height, alpha = 1.0, neighbors = 4L,
           backend = c("rust", "r"), solver = c("cg", "bicgstab")) {
    standardGeneric("graphSmooth")
  }
)

#' @rdname graphSmooth
#' @export
setMethod(
  "graphSmooth",
  signature(x = "hyperSpec"),
  function(x, width, height, alpha = 1.0, neighbors = 4L,
           backend = c("rust", "r"), solver = c("cg", "bicgstab")) {
    backend   <- match.arg(backend)
    solver    <- match.arg(solver)
    width     <- as.integer(width)
    height    <- as.integer(height)
    neighbors <- as.integer(neighbors)
    alpha     <- as.numeric(alpha)

    if (length(width) != 1L || length(height) != 1L ||
        is.na(width) || is.na(height) || width < 1L || height < 1L) {
      stop("`width` and `height` must be positive integers.")
    }
    if (nrow(x) != width * height) {
      stop(sprintf(
        "nrow(x) (%d) must equal width * height (%d).",
        nrow(x), width * height
      ))
    }
    if (!(neighbors %in% c(4L, 8L))) {
      stop("`neighbors` must be 4 or 8.")
    }
    if (length(alpha) != 1L || is.na(alpha) || alpha < 0) {
      stop("`alpha` must be a non-negative scalar.")
    }

    spc <- x@data$spc

    smoothed <- switch(
      backend,
      r    = graph_smooth_r(spc, width, height, alpha, neighbors),
      rust = tryCatch(
        graph_smooth_rust(spc, width, height, alpha, neighbors, solver),
        error = function(e) {
          graph_smooth_r(spc, width, height, alpha, neighbors)
        }
      )
    )

    res <- x
    res@data$spc <- smoothed
    res
  }
)
