library(Matrix)

test_that("Laplacian matrix construction has correct properties", {
  # Small 4x3 grid
  w <- 4L
  h <- 3L
  n <- w * h
  res <- laplacian_matrix_rust(w, h, 4L)
  
  # Reconstruct a sparse matrix in R using Matrix package for easy verification
  L <- Matrix::sparseMatrix(
    p = res$p,
    i = res$i,
    x = res$x,
    dims = c(res$n, res$n),
    index1 = FALSE
  )

  stats <- pixel_graph_stats_rust(w, h, 4L)
  
  # i. diagonal elements equal the node degree
  col_indices <- rep(0:(n-1), diff(L@p))
  diag_vals <- L@x[L@i == col_indices]
  expect_equal(diag_vals, as.numeric(stats$degrees))
  
  # ii. off-diagonal entries equal -1 for neighbors
  off_diag_vals <- L@x[L@i != rep(0:(n-1), diff(L@p))]
  expect_true(all(off_diag_vals == -1))
  
  # In an undirected graph, number of off-diagonal entries = 2 * edges
  expect_equal(length(off_diag_vals), 2 * stats$edges)
  
  # iii. row sums equal zero for the Laplacian
  r_sums <- Matrix::rowSums(L)
  expect_true(all(abs(r_sums) < 1e-14))
})
