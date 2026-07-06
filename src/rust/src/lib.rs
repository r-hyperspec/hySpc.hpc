use extendr_api::prelude::*;
use faer::sparse::SparseColMat;
use dgcmatrix_faer_bridge::{dgcmatrix_to_faer, DgCMatrixView};
use petgraph::graph::{NodeIndex, UnGraph};

/// `y <- A * x` for a matrix `A` stored in compressed-sparse-column (CSC)
/// form (`col_ptrs`, `row_indices`, `values`).
///
/// This is the single mat-vec primitive shared by the iterative
/// solvers. Keeping it self-contained (rather than going through faer's
/// operator overloads) means the Krylov solvers depend only on plain
/// slices, which keeps the FFI kernel small and easy to reason about.
fn csc_matvec(
    n: usize,
    col_ptrs: &[usize],
    row_indices: &[usize],
    values: &[f64],
    x: &[f64],
    y: &mut [f64],
) {
    for yi in y.iter_mut() {
        *yi = 0.0;
    }
    for j in 0..n {
        let xj = x[j];
        for idx in col_ptrs[j]..col_ptrs[j + 1] {
            y[row_indices[idx]] += values[idx] * xj;
        }
    }
}

#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(u, v)| u * v).sum()
}

#[inline]
fn norm2(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Conjugate Gradient (symmetric positive-definite case).
///
/// Solves `A x = b` for a symmetric positive-definite `A` (here
/// `A = I + alpha * L`, which is SPD for `alpha > 0`). Starts from
/// `x0 = 0`, so the initial residual is simply `r0 = b`. Convergence is
/// declared when `||r|| <= tol * max(||b||, 1)`.
fn solve_cg(
    n: usize,
    col_ptrs: &[usize],
    row_indices: &[usize],
    values: &[f64],
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> std::result::Result<Vec<f64>, String> {
    let mut x = vec![0.0; n];
    let mut r = b.to_vec(); // r0 = b - A x0 = b
    let mut p = r.clone();
    let mut ap = vec![0.0; n];

    let threshold = tol * norm2(b).max(1.0);
    let mut rs_old = dot(&r, &r);
    if rs_old.sqrt() <= threshold {
        return Ok(x);
    }

    for _ in 0..max_iter {
        csc_matvec(n, col_ptrs, row_indices, values, &p, &mut ap);
        let pap = dot(&p, &ap);
        if !pap.is_finite() || pap <= 0.0 {
            return Err("CG breakdown: non-positive curvature (matrix not SPD?)".to_string());
        }
        let alpha = rs_old / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rs_new = dot(&r, &r);
        if rs_new.sqrt() <= threshold {
            return Ok(x);
        }
        let beta = rs_new / rs_old;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rs_old = rs_new;
    }

    if !x.iter().all(|v| v.is_finite()) {
        return Err("CG produced non-finite values".to_string());
    }
    Ok(x)
}

/// BiCGSTAB
///
/// Solves `A x = b` for a general square `A`. Provided for completeness
/// per the proposal: the smoothing operator `I + alpha * L` is symmetric
/// so CG is preferred, but BiCGSTAB handles the non-symmetric systems a
/// future weighted / directed graph Laplacian could produce. Starts from
/// `x0 = 0` with shadow residual `r_hat = r0 = b`.
fn solve_bicgstab(
    n: usize,
    col_ptrs: &[usize],
    row_indices: &[usize],
    values: &[f64],
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> std::result::Result<Vec<f64>, String> {
    let mut x = vec![0.0; n];
    let mut r = b.to_vec(); // r0 = b - A x0 = b
    let r_hat = r.clone();

    let threshold = tol * norm2(b).max(1.0);
    if norm2(&r) <= threshold {
        return Ok(x);
    }

    let mut rho = 1.0_f64;
    let mut alpha = 1.0_f64;
    let mut omega = 1.0_f64;
    let mut v = vec![0.0; n];
    let mut p = vec![0.0; n];
    let mut s = vec![0.0; n];
    let mut t = vec![0.0; n];

    for _ in 0..max_iter {
        let rho_new = dot(&r_hat, &r);
        if rho_new == 0.0 || !rho_new.is_finite() {
            return Err("BiCGSTAB breakdown: rho = 0".to_string());
        }
        let beta = (rho_new / rho) * (alpha / omega);
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        csc_matvec(n, col_ptrs, row_indices, values, &p, &mut v);
        let rhv = dot(&r_hat, &v);
        if rhv == 0.0 || !rhv.is_finite() {
            return Err("BiCGSTAB breakdown: <r_hat, v> = 0".to_string());
        }
        alpha = rho_new / rhv;
        for i in 0..n {
            s[i] = r[i] - alpha * v[i];
        }
        if norm2(&s) <= threshold {
            for i in 0..n {
                x[i] += alpha * p[i];
            }
            return Ok(x);
        }
        csc_matvec(n, col_ptrs, row_indices, values, &s, &mut t);
        let tt = dot(&t, &t);
        if tt == 0.0 || !tt.is_finite() {
            return Err("BiCGSTAB breakdown: <t, t> = 0".to_string());
        }
        omega = dot(&t, &s) / tt;
        for i in 0..n {
            x[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }
        if norm2(&r) <= threshold {
            return Ok(x);
        }
        if omega == 0.0 {
            return Err("BiCGSTAB breakdown: omega = 0".to_string());
        }
        rho = rho_new;
    }

    if !x.iter().all(|v| v.is_finite()) {
        return Err("BiCGSTAB produced non-finite values".to_string());
    }
    Ok(x)
}

/// Build the pixel-neighborhood graph for a `width` x `height` image grid.
/// Each pixel is a node and edges connect spatially adjacent pixels. Pixel ordering is column-major, matching the R side:
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

/// Assemble `A = I + alpha * L` in CSC form directly from the pixel
/// graph, where `L = D - W` is the combinatorial Laplacian.
///
/// Column `k` holds the diagonal `1 + alpha * deg(k)` plus a `-alpha`
/// entry for every neighbor of `k`. Row indices within each column are
/// sorted and unique. Returns `(n, col_ptrs, row_indices, values)`.
fn assemble_shifted_laplacian(
    width: usize,
    height: usize,
    neighbors: i32,
    alpha: f64,
) -> (usize, Vec<usize>, Vec<usize>, Vec<f64>) {
    let n = width * height;
    let graph = build_pixel_graph(width, height, neighbors);

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

    (n, col_ptrs, row_indices, values)
}

/// Solve graph-based spatial smoothing: (I + alpha * L) x = b
///
/// The pixel neighborhood graph is reconstructed in
/// Rust with `petgraph` from `width`/`height`/`neighbors` (no adjacency
/// matrix crosses the FFI boundary), the combinatorial Laplacian
/// `L = D - W` is assembled directly in CSC form, and the linear system
/// `(I + alpha * L) x = b` is solved band-by-band with an iterative
/// Krylov solver: Conjugate Gradient (`"cg"`, the default -- the
/// operator is symmetric positive-definite) or BiCGSTAB (`"bicgstab"`).
///
/// @param data RMatrix of input spectra data.
/// @param width width of the image.
/// @param height height of the image.
/// @param alpha smoothing parameter.
/// @param neighbors connection type (4 or 8).
/// @param solver iterative method, `"cg"` or `"bicgstab"`.
/// @export
#[extendr]
fn graph_smooth_rust(
  data: RMatrix<f64>,
  width: usize,
  height: usize,
  alpha: f64,
  neighbors: i32,
  solver: &str,
) -> extendr_api::Result<RMatrix<f64>> {
  if neighbors != 4 && neighbors != 8 {
      return Err(Error::Other(format!(
          "`neighbors` must be 4 or 8, got {}",
          neighbors
      )));
  }
  if solver != "cg" && solver != "bicgstab" {
      return Err(Error::Other(format!(
          "`solver` must be \"cg\" or \"bicgstab\", got {:?}",
          solver
      )));
  }
  let n = width * height;
  if data.nrows() != n {
      return Err(Error::Other(format!("Number of rows in data ({}) must equal width * height ({})", data.nrows(), n)));
  }
   let n_cols = data.ncols();

   // reconstruct the pixel neighborhood graph in Rust with
   // petgraph and assemble A = I + alpha * L (L = D - W) in CSC form.
  let (_, col_ptrs, row_indices, values) =
      assemble_shifted_laplacian(width, height, neighbors, alpha);

   // solve (I + alpha * L) x = b for every wavelength band with
   // the chosen iterative Krylov method. Each column of `data` is an
   // independent right-hand side (we will later parallelize this loop).
  let r_slice = data
      .as_real_slice()
      .ok_or_else(|| Error::Other("`data` must be a numeric (double) matrix".to_string()))?;

  let tol = 1e-10;
  // CG converges in at most n steps in exact arithmetic; the `+ 1000`
  // headroom covers finite-precision drift for tiny grids.
  let max_iter = n + 1000;

  let mut out_data = vec![0.0f64; n * n_cols];
  for jcol in 0..n_cols {
      let b = &r_slice[jcol * n..jcol * n + n];
      let x = match solver {
          "bicgstab" => solve_bicgstab(n, &col_ptrs, &row_indices, &values, b, tol, max_iter),
          _ => solve_cg(n, &col_ptrs, &row_indices, &values, b, tol, max_iter),
      }
      .map_err(|e| {
          Error::Other(format!("{} solver failed on band {}: {}", solver, jcol, e))
      })?;
      out_data[jcol * n..jcol * n + n].copy_from_slice(&x);
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

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graph::NodeIndex;

    #[test]
    fn graph_has_one_node_per_pixel() {
        let g = build_pixel_graph(4, 3, 4);
        assert_eq!(g.node_count(), 12);
    }

    #[test]
    fn edge_count_matches_connectivity() {
        let (w, h) = (4usize, 3usize);
        let g4 = build_pixel_graph(w, h, 4);
        // 4-connectivity: (w-1)*h horizontal + w*(h-1) vertical edges.
        let expected4 = (w - 1) * h + w * (h - 1);
        assert_eq!(g4.edge_count(), expected4);

        let g8 = build_pixel_graph(w, h, 8);
        // 8-connectivity adds 2*(w-1)*(h-1) diagonal edges.
        let diag = 2 * (w - 1) * (h - 1);
        assert_eq!(g8.edge_count(), expected4 + diag);
    }

    #[test]
    fn corner_edge_interior_degrees_4conn() {
        // width = 4, height = 3, column-major k = j*width + i.
        let g = build_pixel_graph(4, 3, 4);
        let deg = |k: usize| g.neighbors(NodeIndex::new(k)).count();
        // corners
        for k in [0usize, 3, 8, 11] {
            assert_eq!(deg(k), 2, "corner {k}");
        }
        // interior pixels (i in {1,2}, j = 1) -> k = 5, 6
        for k in [5usize, 6] {
            assert_eq!(deg(k), 4, "interior {k}");
        }
        // an edge (non-corner border) pixel, e.g. k = 1 (i=1,j=0)
        assert_eq!(deg(1), 3);
    }

    #[test]
    fn shifted_laplacian_rows_are_consistent() {
        // Row sum of L is 0, so row sum of (I + alpha*L) must be 1.
        let (n, col_ptrs, row_indices, values) = assemble_shifted_laplacian(5, 4, 8, 0.7);
        let mut row_sum = vec![0.0f64; n];
        for j in 0..n {
            for idx in col_ptrs[j]..col_ptrs[j + 1] {
                row_sum[row_indices[idx]] += values[idx];
            }
        }
        for (k, s) in row_sum.iter().enumerate() {
            assert!((s - 1.0).abs() < 1e-12, "row {k} sum = {s}");
        }
    }

    // ---- iterative solvers ---------------------------------------

    fn residual_norm(
        n: usize,
        col_ptrs: &[usize],
        row_indices: &[usize],
        values: &[f64],
        x: &[f64],
        b: &[f64],
    ) -> f64 {
        let mut ax = vec![0.0; n];
        csc_matvec(n, col_ptrs, row_indices, values, x, &mut ax);
        norm2(&ax.iter().zip(b).map(|(a, bb)| a - bb).collect::<Vec<_>>())
    }

    #[test]
    fn cg_solves_spd_system() {
        let (n, cp, ri, v) = assemble_shifted_laplacian(6, 5, 4, 1.3);
        let b: Vec<f64> = (0..n).map(|k| ((k * 7 + 1) % 11) as f64 - 5.0).collect();
        let x = solve_cg(n, &cp, &ri, &v, &b, 1e-12, n + 1000).unwrap();
        assert!(residual_norm(n, &cp, &ri, &v, &x, &b) < 1e-8);
    }

    #[test]
    fn bicgstab_solves_same_system_as_cg() {
        let (n, cp, ri, v) = assemble_shifted_laplacian(6, 5, 8, 0.9);
        let b: Vec<f64> = (0..n).map(|k| ((k * 3 + 2) % 13) as f64 - 6.0).collect();
        let xc = solve_cg(n, &cp, &ri, &v, &b, 1e-12, n + 1000).unwrap();
        let xb = solve_bicgstab(n, &cp, &ri, &v, &b, 1e-12, n + 1000).unwrap();
        let diff = norm2(&xc.iter().zip(&xb).map(|(a, c)| a - c).collect::<Vec<_>>());
        assert!(diff < 1e-6, "CG and BiCGSTAB disagree: {diff}");
        assert!(residual_norm(n, &cp, &ri, &v, &xb, &b) < 1e-8);
    }

    #[test]
    fn alpha_zero_is_identity() {
        // A = I, so the solution equals b for any solver.
        let (n, cp, ri, v) = assemble_shifted_laplacian(4, 4, 4, 0.0);
        let b: Vec<f64> = (0..n).map(|k| k as f64).collect();
        let x = solve_cg(n, &cp, &ri, &v, &b, 1e-12, n + 1000).unwrap();
        for i in 0..n {
            assert!((x[i] - b[i]).abs() < 1e-12);
        }
    }
}
