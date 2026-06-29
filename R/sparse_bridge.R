#' R side of the dgCMatrix to faer FFI bridge
#'
#' Helpers that pull the three CSC slots (`@p`, `@i`, `@x`) out
#' of an R sparse matrix from the [Matrix] package and hand them to the
#' Rust extension. The Rust side reconstructs a `faer::SparseColMat`
#' through the `dgcmatrix-faer-bridge` crate, validates the structure,
#' and runs the requested computation.


#' Coerce a sparse Matrix into dgCMatrix and return its CSC slots
#'
#' Internal helper used by the FFI-bridge entry points. Any
#' [Matrix::CsparseMatrix-class] (including `dgCMatrix`, `dsCMatrix`,
#' `dtCMatrix`, `ddiMatrix`) is normalised to `dgCMatrix` so the
#' `@p`, `@i`, `@x` triplet has a uniform layout.
#'
#' @param A a sparse matrix coercible to `dgCMatrix`.
#' @return A list with components `p` (integer column pointers),
#'   `i` (integer row indices), `x` (numeric values), `nrow`, `ncol`.
#' @keywords internal
.dgc_slots <- function(A) {
 if (!methods::is(A, "Matrix") && !is.matrix(A)) {
   stop("`A` must be a matrix or sparse Matrix object.")
 }
 if (!methods::is(A, "dgCMatrix")) {
   A <- methods::as(A, "CsparseMatrix")
   A <- methods::as(A, "generalMatrix")
   A <- methods::as(A, "dgCMatrix")
 }
 list(
   p    = A@p,
   i    = A@i,
   x    = A@x,
   nrow = nrow(A),
   ncol = ncol(A)
 )
}


#' Row sums of a sparse matrix via the Rust dgCMatrix -> faer bridge
#'
#' Stage 2 end-to-end check of the FFI bridge: extract `@p`, `@i`, `@x`
#' from a `dgCMatrix`, ship them to Rust, reconstruct a faer
#' `SparseColMat` through the `dgcmatrix-faer-bridge` crate, and return
#' the row sums computed by walking the reconstructed sparse matrix.
#'
#' The result must agree with [Matrix::rowSums()] for every valid input;
#' this is the basis of the Stage 2 unit tests.
#'
#' @param A a sparse matrix coercible to [Matrix::dgCMatrix-class].
#' @return Numeric vector of length `nrow(A)`.
#' @seealso [dgc_row_sums_rust()], [Matrix::rowSums()].
#' @examples
#' \dontrun{
#'   A <- Matrix::sparseMatrix(
#'     i = c(1, 2, 3, 1, 4),
#'     j = c(1, 2, 3, 4, 4),
#'     x = c(2, 3, 4, 5, 6),
#'     dims = c(4, 4)
#'   )
#'   sparse_row_sums(A)        # via Rust bridge
#'   Matrix::rowSums(A)        # reference
#' }
#' @export
sparse_row_sums <- function(A) {
 s <- .dgc_slots(A)
 dgc_row_sums_rust(
   as.integer(s$p),
   as.integer(s$i),
   as.numeric(s$x),
   as.integer(s$nrow),
   as.integer(s$ncol)
 )
}
