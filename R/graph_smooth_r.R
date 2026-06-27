#' Pure-R baseline for graph-based spatial smoothing
#'
#' Stage 1 baseline implementation of the graph Laplacian smoother
#' used as a correctness reference for the Rust kernel
#' ([graph_smooth_rust()]). Builds the pixel adjacency, the combinatorial
#' graph Laplacian \eqn{L = D - W} on a `width` x `height` grid, and
#' solves the sparse linear system
#'
#' \deqn{(I + \alpha L)\, X = B}
#'
#' band-by-band using [Matrix::solve()]. With `neighbors = 4` or `8`,
#' `L` is symmetric positive semi-definite, so `I + alpha * L` is SPD for
#' any `alpha > 0` and the solve is well-posed.
#'
#' @param data numeric matrix of input spectra. Rows index pixels in
#'   column-major order (`k = j * width + i`, with `i` the x-coordinate
#'   and `j` the y-coordinate, both 0-based), columns index wavelength
#'   bands.
#' @param width,height integer image dimensions; `nrow(data)` must equal
#'   `width * height`.
#' @param alpha non-negative numeric smoothing strength. `alpha = 0`
#'   returns `data` unchanged.
#' @param neighbors integer, either `4` (von Neumann) or `8` (Moore)
#'   neighborhood connectivity.
#'
#' @return Dense numeric matrix of smoothed spectra with the same shape
#'   as `data`.
#'
#' @seealso [graphSmooth()], [graph_smooth_rust()].
#'
#' @examples
#' set.seed(1)
#' B <- matrix(rnorm(16 * 3), nrow = 16, ncol = 3)
#' graph_smooth_r(B, width = 4, height = 4, alpha = 0.5, neighbors = 4L)
#'
#' @export
graph_smooth_r <- function(data, width, height, alpha = 1.0, neighbors = 4L) {
  width  <- as.integer(width)
  height <- as.integer(height)
  alpha  <- as.numeric(alpha)
  neighbors <- as.integer(neighbors)

  if (length(width) != 1L  || is.na(width)  || width  < 1L) {
    stop("`width` must be a positive integer.")
  }
  if (length(height) != 1L || is.na(height) || height < 1L) {
    stop("`height` must be a positive integer.")
  }
  if (length(alpha) != 1L || is.na(alpha) || alpha < 0) {
    stop("`alpha` must be a non-negative scalar.")
  }
  if (!(neighbors %in% c(4L, 8L))) {
    stop("`neighbors` must be 4 or 8.")
  }

  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  storage.mode(data) <- "double"

  n <- width * height
  if (nrow(data) != n) {
    stop(sprintf(
      "nrow(data) (%d) must equal width * height (%d).",
      nrow(data), n
    ))
  }

  if (alpha == 0) {
    return(data)
  }

  L <- build_laplacian_r(width, height, neighbors)
  A <- Matrix::Diagonal(n) + alpha * L
  out <- Matrix::solve(A, data)
  as.matrix(out)
}

#' Build the pixel adjacency matrix for a `width` x `height` grid.
#'
#' Returns a symmetric, binary, sparse `dgCMatrix` of size
#' `(width * height) x (width * height)`.
#'
#' Pixel ordering matches the Rust kernel: `k = j * width + i` with
#' `i = k %% width` (x-coordinate) and `j = k %/% width` (y-coordinate),
#' both 0-based.
#'
#' @inheritParams graph_smooth_r
#' @return A symmetric sparse adjacency matrix.
#' @keywords internal
build_adjacency_r <- function(width, height, neighbors = 4L) {
  width  <- as.integer(width)
  height <- as.integer(height)
  neighbors <- as.integer(neighbors)

  offsets <- if (neighbors == 8L) {
    list(
      c(-1L, -1L), c( 0L, -1L), c( 1L, -1L),
      c(-1L,  0L),              c( 1L,  0L),
      c(-1L,  1L), c( 0L,  1L), c( 1L,  1L)
    )
  } else if (neighbors == 4L) {
    list(c(0L, -1L), c(-1L, 0L), c(1L, 0L), c(0L, 1L))
  } else {
    stop("`neighbors` must be 4 or 8.")
  }

  n  <- width * height
  ii <- rep(seq.int(0L, width  - 1L), times = height)  # x-coord
  jj <- rep(seq.int(0L, height - 1L), each  = width)   # y-coord
  k  <- jj * width + ii                                # node index

  rows <- integer(0)
  cols <- integer(0)
  for (off in offsets) {
    ni <- ii + off[1L]
    nj <- jj + off[2L]
    keep <- ni >= 0L & ni < width & nj >= 0L & nj < height
    rows <- c(rows, k[keep])
    cols <- c(cols, (nj[keep] * width + ni[keep]))
  }

  Matrix::sparseMatrix(
    i = rows + 1L,
    j = cols + 1L,
    x = 1,
    dims = c(n, n),
    symmetric = FALSE,
    check = FALSE
  )
}

#' Build the combinatorial graph Laplacian `L = D - W` for a grid.
#'
#' @inheritParams graph_smooth_r
#' @return A symmetric sparse Laplacian as a `dgCMatrix`.
#' @keywords internal
build_laplacian_r <- function(width, height, neighbors = 4L) {
  W <- build_adjacency_r(width, height, neighbors)
  d <- Matrix::rowSums(W)
  Matrix::Diagonal(x = d) - W
}
