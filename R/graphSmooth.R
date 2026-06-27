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
#'   ([graph_smooth_rust()]) built on `faer`. Requires the package's
#'   Rust extension to be compiled. Silently falls back to "r" if unavailable or fails.
#' * `backend = "r"`: the pure-R baseline
#'   ([graph_smooth_r()]) built on [Matrix] sparse routines. Used as the
#'   correctness reference in Stage 1.
#'
#' Pixel ordering is column-major: spectrum row `k` corresponds to grid
#' position `i = k %% width` (x), `j = k %/% width` (y), both 0-based.
#'
#' @param x [hyperSpec::hyperSpec] object whose `nrow()` equals
#'   `width * height`.
#' @param width,height integer image dimensions.
#' @param alpha non-negative numeric smoothing strength. Larger values
#'   produce smoother output; `alpha = 0` returns `x` unchanged.
#' @param neighbors integer, either `4` (von Neumann) or `8` (Moore)
#'   neighborhood connectivity.
#' @param backend one of `"r"` or `"rust"`; selects the solver
#'   implementation.
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
           backend = c("rust", "r")) {
    standardGeneric("graphSmooth")
  }
)

#' @rdname graphSmooth
#' @export
setMethod(
  "graphSmooth",
  signature(x = "hyperSpec"),
  function(x, width, height, alpha = 1.0, neighbors = 4L,
           backend = c("rust", "r")) {
    backend   <- match.arg(backend)
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
        graph_smooth_rust(spc, width, height, alpha, neighbors),
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
