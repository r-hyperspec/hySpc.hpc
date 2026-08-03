library(Matrix)

test_that("Laplacian matrix construction has correct properties", {
  # Small 4x3 grid
  w <- 4L
  h <- 3L
  n <- w * h
  
  # Fetch Laplacian representation from Rust bridge
  # L = D - W
  res <- laplacian_matrix_rust(w, h, 4L)
  
  # Reconstruct a sparse matrix in R using Matrix package for easy verification
  # Note: `i` and `p` returned from Rust are 0-indexed, so we use index1 = FALSE
  L <- Matrix::sparseMatrix(
    p = res$p,
    i = res$i,
    x = res$x,
    dims = c(res$n, res$n),
    index1 = FALSE
  )
  
  # Also get degrees to verify diagonal
  stats <- pixel_graph_stats_rust(w, h, 4L)
  
  # i. diagonal elements equal the node degree
  # We extract the diagonal directly from the CSC structure to avoid S4 method dispatch issues
  col_indices <- rep(0:(n-1), diff(L@p))
  diag_vals <- L@x[L@i == col_indices]
  expect_equal(diag_vals, as.numeric(stats$degrees))
  
  # ii. off-diagonal entries equal -1 for neighbors
  # We check that all off-diagonal non-zero entries are exactly -1
  # The diagonal entries are positive, all other non-zero entries must be -1
  off_diag_vals <- L@x[L@i != rep(0:(n-1), diff(L@p))]
  expect_true(all(off_diag_vals == -1))
  
  # Additionally, let's verify there are the correct number of off-diagonal -1s
  # In an undirected graph, number of off-diagonal entries = 2 * edges
  expect_equal(length(off_diag_vals), 2 * stats$edges)
  
  # iii. row sums equal zero for the Laplacian
  # A valid combinatorial Laplacian L = D - W has row sums exactly 0
  r_sums <- Matrix::rowSums(L)
  expect_true(all(abs(r_sums) < 1e-14))
})
