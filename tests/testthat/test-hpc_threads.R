test_that("detect_hpc_cores respects scheduler environment variables", {
  # Slurm mock
  withr::with_envvar(c(SLURM_CPUS_PER_TASK = "8"), {
    expect_equal(detect_hpc_cores(), 8L)
  })

  # PBS mock
  withr::with_envvar(c(SLURM_CPUS_PER_TASK = "", PBS_NCPUS = "16"), {
    expect_equal(detect_hpc_cores(), 16L)
  })

  # SGE mock
  withr::with_envvar(c(SLURM_CPUS_PER_TASK = "", PBS_NCPUS = "", NSLOTS = "32"), {
    expect_equal(detect_hpc_cores(), 32L)
  })

  # Fallback to physical detection
  withr::with_envvar(c(SLURM_CPUS_PER_TASK = "", PBS_NCPUS = "", NSLOTS = "", OMP_NUM_THREADS = "", RAYON_NUM_THREADS = ""), {
    cores <- detect_hpc_cores()
    expect_true(is.integer(cores))
    expect_true(cores >= 1L)
  })
})

test_that("get_hpc_threads and set_hpc_threads work correctly", {
  init <- get_hpc_threads()
  expect_true(is.numeric(init))
  expect_true(init >= 1L)

  set_hpc_threads(2L)
  expect_equal(get_hpc_threads(), 2L)

  set_hpc_threads(4L)
  expect_equal(get_hpc_threads(), 4L)

  # Invalid inputs error
  expect_error(set_hpc_threads(0L))
  expect_error(set_hpc_threads(-2L))
  expect_error(set_hpc_threads("invalid"))

  # Restore initial threads
  set_hpc_threads(init)
  expect_equal(get_hpc_threads(), init)
})

test_that("with_hpc_threads scopes thread budget and restores upon exit", {
  set_hpc_threads(2L)
  expect_equal(get_hpc_threads(), 2L)

  res <- with_hpc_threads(4L, {
    expect_equal(get_hpc_threads(), 4L)
    "success"
  })
  expect_equal(res, "success")
  expect_equal(get_hpc_threads(), 2L)

  # Restores even on error
  expect_error(
    with_hpc_threads(1L, {
      expect_equal(get_hpc_threads(), 1L)
      stop("intentional error")
    })
  )
  expect_equal(get_hpc_threads(), 2L)
})

test_that("graphSmooth produces identical output across thread configurations", {
  width <- 4
  height <- 4
  n <- width * height
  n_lambda <- 10

  set.seed(123)
  spc <- matrix(rnorm(n * n_lambda), nrow = n)
  obj <- new("hyperSpec", spc = spc)

  # Solve with 1 thread
  res1 <- with_hpc_threads(1L, {
    graphSmooth(obj, width = width, height = height, alpha = 0.8, backend = "rust", solver = "cg")
  })

  # Solve with 2 threads
  res2 <- with_hpc_threads(2L, {
    graphSmooth(obj, width = width, height = height, alpha = 0.8, backend = "rust", solver = "cg")
  })

  # Solve with 4 threads
  res4 <- with_hpc_threads(4L, {
    graphSmooth(obj, width = width, height = height, alpha = 0.8, backend = "rust", solver = "cg")
  })

  expect_equal(res1@data$spc, res2@data$spc, tolerance = 1e-9)
  expect_equal(res1@data$spc, res4@data$spc, tolerance = 1e-9)
})
