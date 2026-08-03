test_that("extraction of dgCMatrix slots works correctly", {
  # Create a valid dgCMatrix
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 1, 4),
    j = c(1, 2, 3, 4, 4),
    x = c(2, 3, 4, 5, 6),
    dims = c(4, 4)
  )

  # .dgc_slots should extract p, i, x, nrow, ncol correctly
  slots <- .dgc_slots(A)
  expect_equal(slots$p, A@p)
  expect_equal(slots$i, A@i)
  expect_equal(slots$x, A@x)
  expect_equal(slots$nrow, 4L)
  expect_equal(slots$ncol, 4L)

  # The main wrapper should compute the correct result
  expect_equal(sparse_row_sums(A), Matrix::rowSums(A))
})

test_that("Rust bridge validates length(i) == length(x)", {
  expect_error(
    dgc_row_sums_rust(
      p = as.integer(c(0, 1, 2)),
      i = as.integer(c(0, 1)),
      x = as.numeric(c(1.0)), # mismatched length
      nrow = 2L,
      ncol = 2L
    ),
    "dgCMatrix -> faer conversion failed:"
  )
})

test_that("Rust bridge validates length(p) == ncol + 1", {
  expect_error(
    dgc_row_sums_rust(
      p = as.integer(c(0, 1)), # ncol = 2 requires length 3
      i = as.integer(c(0, 1)),
      x = as.numeric(c(1.0, 2.0)),
      nrow = 2L,
      ncol = 2L
    ),
    "dgCMatrix -> faer conversion failed:"
  )
})

test_that("Rust bridge validates column pointers are monotonic", {
  expect_error(
    dgc_row_sums_rust(
      p = as.integer(c(0, 2, 1)), # non-monotonic
      i = as.integer(c(0, 1)),
      x = as.numeric(c(1.0, 2.0)),
      nrow = 2L,
      ncol = 2L
    ),
    "dgCMatrix -> faer conversion failed:"
  )
})

test_that("Rust bridge validates row indices fall within valid bounds", {
  expect_error(
    dgc_row_sums_rust(
      p = as.integer(c(0, 1, 2)),
      i = as.integer(c(0, 3)), # row 3 is out of bounds for nrow=2
      x = as.numeric(c(1.0, 2.0)),
      nrow = 2L,
      ncol = 2L
    ),
    "Assertion failed"
  )
})

test_that("Rust bridge catches type errors and negative dimensions", {
  expect_error(
    dgc_row_sums_rust(
      p = as.integer(c(0, 1, 2)),
      i = as.integer(c(0, 1)),
      x = as.numeric(c(1.0, 2.0)),
      nrow = -1L, # invalid negative dimension
      ncol = 2L
    ),
    "nrow and ncol must be non-negative"
  )

  expect_error(
    dgc_row_sums_rust(
      p = c(0, 1, 2), # double instead of integer
      i = as.integer(c(0, 1)),
      x = as.numeric(c(1.0, 2.0)),
      nrow = 2L,
      ncol = 2L
    ),
    "`p` must be an integer vector"
  )
})
