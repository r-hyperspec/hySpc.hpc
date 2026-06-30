use extendr_api::prelude::*;
use faer::sparse::SparseColMat;
use faer::Mat;
use faer::solvers::SpSolver;
use dgcmatrix_faer_bridge::{dgcmatrix_to_faer, DgCMatrixView};
use petgraph::graph::{NodeIndex, UnGraph};

/// Build the pixel-neighborhood graph for a `width` x `height` image grid.
/// Stage 3: each pixel is a node and edges connect spatially adjacent pixels. Pixel ordering is column-major, matching the R side:
/// node `k` sits at `i = k % width` (x), `j = k / width` (y), both 0-based. Connectivity is configurable:
/// * `neighbors = 4` (von Neumann): horizontal + vertical edges.
/// * `neighbors = 8` (Moore): also the two diagonals.
fn build_pixel_graph(width: usize, height: usize, neighbors: i32) -> UnGraph<(), ()> {
    let n = width * height;
    let mut g: UnGraph<(), ()> = UnGraph::with_capacity(n, n * 4);
    for _ in 0..n {
        g.add_node(());
    }

    for k in 0..n {
        let i = k % width;
        let j = k / width;

        // Horizontal + vertical (4-connectivity), forward direction only.
        if i + 1 < width {
            g.add_edge(NodeIndex::new(k), NodeIndex::new(k + 1), ());
        }
        if j + 1 < height {
            g.add_edge(NodeIndex::new(k), NodeIndex::new(k + width), ());
        }

        // Diagonals (8-connectivity), forward direction only.
        if neighbors == 8 {
            if i + 1 < width && j + 1 < height {
                g.add_edge(NodeIndex::new(k), NodeIndex::new(k + width + 1), ());
            }
            if i > 0 && j + 1 < height {
                g.add_edge(NodeIndex::new(k), NodeIndex::new(k + width - 1), ());
            }
        }
    }

    g
}

/// Solve graph-based spatial smoothing: (I + alpha * L) x = b
///
/// Stage 3: the pixel neighborhood graph is reconstructed in Rust with petgraph` from `width`/`height`/`neighbors` (no adjacency matrix
/// crosses the FFI boundary), and the combinatorial Laplacian L = D - W` is assembled directly as a `faer` `SparseColMat`. The
/// linear system `(I + alpha * L) x = b` is then solved band-by-band.
/// The solve currently uses a sparse Cholesky factorization; Stage 4 will replace this with the iterative CG / BiCGSTAB solvers.
#[extendr]
fn graph_smooth_rust(
  data: RMatrix<f64>,
  width: usize,
  height: usize,
  alpha: f64,
  neighbors: i32,
) -> extendr_api::Result<RMatrix<f64>> {
  if neighbors != 4 && neighbors != 8 {
      return Err(Error::Other(format!(
          "`neighbors` must be 4 or 8, got {}",
          neighbors
      )));
  }
  let n = width * height;
  if data.nrows() != n {
      return Err(Error::Other(format!("Number of rows in data ({}) must equal width * height ({})", data.nrows(), n)));
  }
   let n_cols = data.ncols();

   // Stage 3: reconstruct the pixel neighborhood graph in Rust.
  let graph = build_pixel_graph(width, height, neighbors);

   // Assemble (I + alpha * L) in CSC form straight from the graph, where
   // L = D - W. Column k holds the diagonal `1 + alpha * deg(k)` plus a
   // `-alpha` entry for every neighbor of k. faer requires the row
   // indices within each column to be sorted and unique, so we collect
   // and sort per column.
  let mut col_ptrs = Vec::with_capacity(n + 1);
  let mut row_indices = Vec::new();
  let mut values = Vec::new();
   for k in 0..n {
      col_ptrs.push(row_indices.len());

      let mut entries: Vec<(usize, f64)> = Vec::new();
      let mut degree = 0.0;
      for nb in graph.neighbors(NodeIndex::new(k)) {
          entries.push((nb.index(), -alpha));
          degree += 1.0;
      }
      entries.push((k, 1.0 + alpha * degree));
      entries.sort_by_key(|&(r, _)| r);

      for (r, v) in entries {
          row_indices.push(r);
          values.push(v);
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
/// Extract `p`, `i`, `x` slots of a `dgCMatrix` on the R side, transfer them through extendr, reconstruct a faer SparseColMat` via the `dgcmatrix-faer-bridge` crate,
///  and return the row sums to R. This is the minimal round-trip that exercises the dgCMatrix -> faer FFI bridge end-to-end.
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
