use extendr_api::prelude::*;
use faer::sparse::SparseColMat;
use faer::Mat;
use faer::solvers::SpSolver;

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

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod hy_spc_hpc;
    fn graph_smooth_rust;
}
