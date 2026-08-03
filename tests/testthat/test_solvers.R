library(Matrix)

test_that("CG solver agrees with exact dense solver", {
  w <- 5L
  h <- 4L
  n <- w * h
  alpha <- 1.5
  
  # Fetch combinatorial Laplacian L
  res <- laplacian_matrix_rust(w, h, 4L)
  L <- Matrix::sparseMatrix(
    p = res$p,
    i = res$i,
    x = res$x,
    dims = c(res$n, res$n),
    index1 = FALSE
  )
  
  # Construct the smoothing operator A = I + alpha * L
  A <- Matrix::Diagonal(n) + alpha * L
  # Generate some random spectra data (3 bands / columns)
  set.seed(42)
  b <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  
  x_exact <- as.matrix(solve(A, b))
  x_cg <- graph_smooth_rust(b, w, h, alpha, 4L, "cg")
  
  # Verify agreement within a tight tolerance (1e-6)
  expect_equal(x_cg, x_exact, tolerance = 1e-6)
})

test_that("BiCGSTAB solver agrees with exact dense solver", {
  w <- 5L
  h <- 4L
  n <- w * h
  alpha <- 0.8
  
  # Use 8-connectivity for variety
  res <- laplacian_matrix_rust(w, h, 8L)
  L <- Matrix::sparseMatrix(
    p = res$p,
    i = res$i,
    x = res$x,
    dims = c(res$n, res$n),
    index1 = FALSE
  )
  
  A <- Matrix::Diagonal(n) + alpha * L
  
  set.seed(123)
  b <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  x_exact <- as.matrix(solve(A, b))
  
  # Rust BiCGSTAB solver
  x_bicg <- graph_smooth_rust(b, w, h, alpha, 8L, "bicgstab")
  
  expect_equal(x_bicg, x_exact, tolerance = 1e-6)
})

test_that("Zero alpha returns exact input", {
  w <- 3L
  h <- 3L
  n <- w * h
  
  b <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  
  # Alpha = 0 means A = I, so output should be perfectly identical to input
  x_cg <- graph_smooth_rust(b, w, h, 0.0, 4L, "cg")
  x_bicg <- graph_smooth_rust(b, w, h, 0.0, 8L, "bicgstab")
  
  expect_equal(x_cg, b, tolerance = 1e-12)
  expect_equal(x_bicg, b, tolerance = 1e-12)
})
