test_that("graphSmooth works correctly with 4-connectivity", {
  width <- 4
  height <- 3
  n <- width * height
  n_lambda <- 2
  
  set.seed(42)
  spc <- matrix(rnorm(n * n_lambda), nrow = n)
  
  alpha <- 0.5
  L <- matrix(0, n, n)
  for (k in 1:n) {
    i <- (k - 1) %% width + 1
    j <- (k - 1) %/% width + 1
    nbrs <- c()
    if (i > 1) nbrs <- c(nbrs, k - 1)
    if (i < width) nbrs <- c(nbrs, k + 1)
    if (j > 1) nbrs <- c(nbrs, k - width)
    if (j < height) nbrs <- c(nbrs, k + width)
    
    L[k, k] <- length(nbrs)
    L[k, nbrs] <- -1
  }
  
  A <- diag(n) + alpha * L
  expected_spc <- solve(A, spc)
  
  # Ensure the package is loaded
  require(hySpc.hpc)
  
  obj <- new("hyperSpec", spc = spc)
  res <- graphSmooth(obj, width, height, alpha, neighbors = 4L)
  
  expect_equal(unname(res@data$spc), unname(expected_spc), tolerance = 1e-6)
})

test_that("graphSmooth works correctly with 8-connectivity", {
  width <- 3
  height <- 4
  n <- width * height
  n_lambda <- 2
  
  set.seed(43)
  spc <- matrix(rnorm(n * n_lambda), nrow = n)
  
  alpha <- 0.2
  L <- matrix(0, n, n)
  for (k in 1:n) {
    i <- (k - 1) %% width + 1
    j <- (k - 1) %/% width + 1
    nbrs <- c()
    
    for (dj in c(-1, 0, 1)) {
      for (di in c(-1, 0, 1)) {
        if (di == 0 && dj == 0) next
        ni <- i + di
        nj <- j + dj
        if (ni >= 1 && ni <= width && nj >= 1 && nj <= height) {
          nbrs <- c(nbrs, (nj - 1) * width + ni)
        }
      }
    }
    
    L[k, k] <- length(nbrs)
    L[k, nbrs] <- -1
  }
  
  A <- diag(n) + alpha * L
  expected_spc <- solve(A, spc)
  
  obj <- new("hyperSpec", spc = spc)
  res <- graphSmooth(obj, width, height, alpha, neighbors = 8L)
  
  expect_equal(unname(res@data$spc), unname(expected_spc), tolerance = 1e-6)
})
