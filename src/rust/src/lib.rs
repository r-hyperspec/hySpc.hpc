use extendr_api::prelude::*;
use faer::sparse::SparseColMat;
use faer::Mat;
use faer::solvers::SpSolver;
use dgcmatrix_faer_bridge::{dgcmatrix_to_faer, DgCMatrixView};

/// Solve graph-based spatial smoothing: (I + alpha * L) x = b
/// @param data RMatrix of input spectra data.
/// @param width width of the image.
/// @param height height of the image.
/// @param alpha smoothing parameter.
/// @param neighbors connection type (4 or 8).
/// @export
#[extendr]
fn graph_smooth_rust(
  data: RMatrix<f64>,
  width: usize,
  height: usize,
  alpha: f64,
  neighbors: i32,
) -> extendr_api::Result<RMatrix<f64>> {
  let n = width * height;
  if data.nrows() != n {
      return Err(Error::Other(format!("Number of rows in data ({}) must equal width * height ({})", data.nrows(), n)));
  }
   let n_cols = data.ncols();
   // Construct the CSC components for (I + alpha * L)
  let mut col_ptrs = Vec::with_capacity(n + 1);
  let mut row_indices = Vec::new();
  let mut values = Vec::new();
   for k in 0..n {
      col_ptrs.push(row_indices.len());
      let i = k % width;
      let j = k / width;
    
      let mut nbrs = Vec::new();
    
      // Define neighbor offsets based on connectivity
      if neighbors == 8 {
          if j > 0 && i > 0 { nbrs.push(k - width - 1); }
          if j > 0 { nbrs.push(k - width); }
          if j > 0 && i + 1 < width { nbrs.push(k - width + 1); }
          if i > 0 { nbrs.push(k - 1); }
          nbrs.push(k); // diagonal
          if i + 1 < width { nbrs.push(k + 1); }
          if j + 1 < height && i > 0 { nbrs.push(k + width - 1); }
          if j + 1 < height { nbrs.push(k + width); }
          if j + 1 < height && i + 1 < width { nbrs.push(k + width + 1); }
      } else {
          // Default 4 connectivity
          if j > 0 { nbrs.push(k - width); }
          if i > 0 { nbrs.push(k - 1); }
          nbrs.push(k); // diagonal
          if i + 1 < width { nbrs.push(k + 1); }
          if j + 1 < height { nbrs.push(k + width); }
      }
    
      let degree = nbrs.len() as f64 - 1.0;
    
      for &m in &nbrs {
          row_indices.push(m);
          if m == k {
              values.push(1.0 + alpha * degree);
          } else {
              values.push(-alpha);
          }
      }
  }
  col_ptrs.push(row_indices.len());
   // Create the faer SparseColMat
  let symbolic = faer::sparse::SymbolicSparseColMat::new_checked(
      n,
      n,
      col_ptrs,
      None,
      row_indices,
  );
   let L_mat = SparseColMat::<usize, f64>::new(
      symbolic,
      values,
  );
   // Decompose using simplicial LLT (Cholesky)
  let cholesky = L_mat.sp_cholesky(faer::Side::Lower).map_err(|e| Error::Other(format!("Cholesky factorization failed: {:?}", e)))?;
   // Copy input data into faer Mat
  let mut b = Mat::<f64>::zeros(n, n_cols);
  let r_slice = data.as_real_slice().unwrap();
  for j in 0..n_cols {
      for i in 0..n {
          b.write(i, j, r_slice[j * n + i]);
      }
  }
   // Solve
  let x = cholesky.solve(&b);
   // Write back to RMatrix
  let mut out_data = vec![0.0; n * n_cols];
  for j in 0..n_cols {
      for i in 0..n {
          out_data[j * n + i] = x.read(i, j);
      }
  }
   Ok(RMatrix::new_matrix(n, n_cols, |r, c| out_data[c * n + r]))
}


/// Compute row sums of an R dgCMatrix using the Rust faer bridge.
///
/// Extract `p`, `i`, `x` slots of a `dgCMatrix` on
/// the R side, transfer them through extendr, reconstruct a faer
/// `SparseColMat` via the `dgcmatrix-faer-bridge` crate, and return the
/// row sums to R. This is the minimal round-trip that exercises the
/// dgCMatrix -> faer FFI bridge end-to-end.
///
/// @param p integer vector, dgCMatrix `@p` column pointers (length ncol + 1).
/// @param i integer vector, dgCMatrix `@i` row indices (length nnz).
/// @param x numeric vector, dgCMatrix `@x` values (length nnz).
/// @param nrow integer, number of rows of the matrix.
/// @param ncol integer, number of columns of the matrix.
/// @return numeric vector of length `nrow` containing the row sums.
/// @export
#[extendr]
fn dgc_row_sums_rust(
  p: Robj,
  i: Robj,
  x: Robj,
  nrow: i32,
  ncol: i32,
) -> extendr_api::Result<Doubles> {
  if nrow < 0 || ncol < 0 {
      return Err(Error::Other(format!(
          "nrow and ncol must be non-negative, got nrow = {}, ncol = {}",
          nrow, ncol
      )));
  }
  let nrow_u = nrow as usize;
  let ncol_u = ncol as usize;


  // Zero-copy borrows of the R vectors backing the dgCMatrix slots.
  // `as_integer_slice` / `as_real_slice` return `None` if the underlying
  // SEXP isn't INTSXP / REALSXP; we translate that to a clear error
  // (passing e.g. `dgCMatrix@p` already coerced to double would otherwise
  // break pointer arithmetic in Rust).
  let p_slice = p
      .as_integer_slice()
      .ok_or_else(|| Error::Other("`p` must be an integer vector".to_string()))?;
  let i_slice = i
      .as_integer_slice()
      .ok_or_else(|| Error::Other("`i` must be an integer vector".to_string()))?;
  let x_slice = x
      .as_real_slice()
      .ok_or_else(|| Error::Other("`x` must be a numeric (double) vector".to_string()))?;

  let view = DgCMatrixView::new(nrow_u, ncol_u, p_slice, i_slice, x_slice);
  let _mat: SparseColMat<usize, f64> = dgcmatrix_to_faer(&view)
      .map_err(|e| Error::Other(format!("dgCMatrix -> faer conversion failed: {}", e)))?;

  let mut out = vec![0.0f64; nrow_u];
  for j in 0..ncol_u {
      let start = p_slice[j] as usize;
      let end = p_slice[j + 1] as usize;
      for k in start..end {
          out[i_slice[k] as usize] += x_slice[k];
      }
  }
  Ok(Doubles::from_values(out))
}

extendr_module! {
  mod hy_spc_hpc;
  fn graph_smooth_rust;
  fn dgc_row_sums_rust;
}
