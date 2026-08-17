#' HPC Thread Management and Scheduler Quota Detection
#'
#' Functions to dynamically query, configure, and scope the worker thread
#' allocation used by the Rust high-performance computing backend in `hySpc.hpc`.
#'
#' @details
#' In multi-tenant HPC cluster environments (Slurm, PBS, Grid Engine) or shared
#' multi-core cloud instances, multiple R worker processes frequently execute on
#' the same physical node. Setting environment variables like `RAYON_NUM_THREADS`
#' after Rayon has already initialized has no effect in pure Rust.
#'
#' `hySpc.hpc` solves this via a dedicated thread pool manager in Rust:
#' * [set_hpc_threads()] dynamically adjusts the active Rayon worker pool size
#'   at runtime without restarting the R process.
#' * [get_hpc_threads()] inspects the currently allocated thread budget.
#' * [detect_hpc_cores()] automatically parses scheduler environment variables
#'   (`SLURM_CPUS_PER_TASK`, `PBS_NCPUS`, `NSLOTS`) to determine the exact CPU
#'   allocation granted to the current job.
#' * [with_hpc_threads()] provides an RAII-style scoped context manager that
#'   temporarily executes code under a specific thread budget and restores the
#'   previous configuration upon exit.
#'
#' @name hpc_threads
NULL

#' Detect CPU core allocation from cluster schedulers or hardware
#'
#' Inspects environment variables set by HPC cluster schedulers (Slurm, PBS,
#' SGE/Grid Engine) as well as OpenMP and Rayon variables. If no scheduler quota
#' is found, falls back to physical core detection via [parallel::detectCores()].
#'
#' @return An integer representing the recommended number of worker threads.
#' @seealso [set_hpc_threads()], [get_hpc_threads()].
#' @examples
#' cores <- detect_hpc_cores()
#' message("Detected allocated cores: ", cores)
#' @export
detect_hpc_cores <- function() {
  # 1. Slurm job allocation
  slurm <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (nzchar(slurm)) {
    val <- suppressWarnings(as.integer(slurm))
    if (!is.na(val) && val >= 1L) return(val)
  }

  # 2. PBS Professional / OpenPBS
  pbs <- Sys.getenv("PBS_NCPUS", unset = "")
  if (nzchar(pbs)) {
    val <- suppressWarnings(as.integer(pbs))
    if (!is.na(val) && val >= 1L) return(val)
  }

  # 3. Sun Grid Engine / OGE / Son of Grid Engine
  sge <- Sys.getenv("NSLOTS", unset = "")
  if (nzchar(sge)) {
    val <- suppressWarnings(as.integer(sge))
    if (!is.na(val) && val >= 1L) return(val)
  }

  # 4. Standard OpenMP / Rayon environment variables
  for (var in c("OMP_NUM_THREADS", "RAYON_NUM_THREADS")) {
    env_val <- Sys.getenv(var, unset = "")
    if (nzchar(env_val)) {
      val <- suppressWarnings(as.integer(env_val))
      if (!is.na(val) && val >= 1L) return(val)
    }
  }

  # 5. Hardware fallback
  detected <- parallel::detectCores(logical = FALSE)
  if (is.na(detected) || detected < 1L) {
    detected <- parallel::detectCores(logical = TRUE)
  }
  if (is.na(detected) || detected < 1L) {
    detected <- 1L
  }
  as.integer(detected)
}

#' Query active worker threads in the Rust backend
#'
#' @return An integer representing the active number of Rayon worker threads.
#' @seealso [set_hpc_threads()], [with_hpc_threads()].
#' @examples
#' \dontrun{
#'   get_hpc_threads()
#' }
#' @export
get_hpc_threads <- function() {
  get_hpc_threads_rust()
}

#' Configure worker threads for the Rust backend
#'
#' Dynamically configures the Rayon thread pool budget in the native Rust
#' backend. If `n` is omitted or `NULL`, the thread count is automatically
#' determined via [detect_hpc_cores()].
#'
#' @param n integer number of worker threads (must be >= 1). If `NULL` (the default),
#'   automatically detects the allocated core quota.
#' @return An integer representing the newly configured thread count (invisibly).
#' @seealso [get_hpc_threads()], [with_hpc_threads()], [detect_hpc_cores()].
#' @examples
#' \dontrun{
#'   set_hpc_threads(4)
#'   get_hpc_threads()
#' }
#' @export
set_hpc_threads <- function(n = NULL) {
  if (is.null(n)) {
    n <- detect_hpc_cores()
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1L) {
    stop("`n` must be a single positive integer >= 1.")
  }
  res <- set_hpc_threads_rust(as.integer(n))
  invisible(res)
}

#' Execute code within a scoped thread budget
#'
#' Temporarily configures the Rust worker thread pool to `n` threads for the
#' duration of evaluating `code`, and automatically restores the previous thread
#' allocation when `code` finishes (even if an error occurs).
#'
#' @param n integer number of worker threads to use during execution.
#' @param code expression or code block to evaluate.
#' @return The result of evaluating `code`.
#' @seealso [set_hpc_threads()], [get_hpc_threads()].
#' @examples
#' \dontrun{
#'   res <- with_hpc_threads(2, {
#'     graphSmooth(spc_cube, width = 50, height = 50, alpha = 0.5)
#'   })
#' }
#' @export
with_hpc_threads <- function(n, code) {
  old_threads <- get_hpc_threads()
  on.exit(set_hpc_threads(old_threads), add = TRUE)
  set_hpc_threads(n)
  force(code)
}
