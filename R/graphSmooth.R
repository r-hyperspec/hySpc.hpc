#' Spatial Graph Smoothing
#'
#' Applies a spatial graph Laplacian smoother to a hyperSpec object representing an image.
#'
#' @param x hyperSpec object
#' @param width integer, width of the image
#' @param height integer, height of the image
#' @param alpha numeric, smoothing parameter (higher = more smoothing)
#' @param neighbors integer, 4 or 8 connectivity
#' @return smoothed hyperSpec object
#' @export
setGeneric("graphSmooth", function(x, width, height, alpha = 1.0, neighbors = 4L) {
  standardGeneric("graphSmooth")
})

#' @rdname graphSmooth
#' @export
setMethod("graphSmooth", signature(x = "hyperSpec"), function(x, width, height, alpha = 1.0, neighbors = 4L) {
  if (nrow(x) != width * height) {
    stop("The number of spectra in x must equal width * height")
  }
  
  if (neighbors != 4L && neighbors != 8L) {
    stop("neighbors must be 4 or 8")
  }
  
  # Call Rust kernel
  # The spectra are stored as x@data$spc, which is a matrix
  spc <- x@data$spc
  
  smoothed_spc <- graph_smooth_rust(spc, as.integer(width), as.integer(height), as.numeric(alpha), as.integer(neighbors))
  
  # Update hyperSpec object
  res <- x
  res@data$spc <- smoothed_spc
  
  return(res)
})
