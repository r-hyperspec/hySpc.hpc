test_that("pixel graph construction has correct number of nodes", {
  # 4x3 grid = 12 nodes
  stats <- pixel_graph_stats_rust(4L, 3L, 4L)
  expect_equal(stats$nodes, 12L)
})

test_that("pixel graph construction has correct number of edges", {
  # Grid: width = 4, height = 3
  w <- 4L
  h <- 3L
  
  # 4-connectivity: (w-1)*h + w*(h-1) edges
  stats_4 <- pixel_graph_stats_rust(w, h, 4L)
  expected_4 <- (w - 1L) * h + w * (h - 1L)
  expect_equal(stats_4$edges, expected_4)
  
  # 8-connectivity: 4-connectivity + 2*(w-1)*(h-1) edges
  stats_8 <- pixel_graph_stats_rust(w, h, 8L)
  expected_8 <- expected_4 + 2L * (w - 1L) * (h - 1L)
  expect_equal(stats_8$edges, expected_8)
})

test_that("interior, edge, and corner pixels have correct neighbor counts (4-connectivity)", {
  # 4x3 grid, column-major indexing k = j * width + i
  # (0,0) (1,0) (2,0) (3,0)  ->  0, 1, 2, 3
  # (0,1) (1,1) (2,1) (3,1)  ->  4, 5, 6, 7
  # (0,2) (1,2) (2,2) (3,2)  ->  8, 9, 10, 11
  
  stats <- pixel_graph_stats_rust(4L, 3L, 4L)
  degrees <- stats$degrees
  
  # R indexing is 1-based, so node k -> degrees[k+1]
  
  # Corners: (0,0)->0, (3,0)->3, (0,2)->8, (3,2)->11
  expect_equal(degrees[0 + 1], 2L)
  expect_equal(degrees[3 + 1], 2L)
  expect_equal(degrees[8 + 1], 2L)
  expect_equal(degrees[11 + 1], 2L)
  
  # Interior pixels: (1,1)->5, (2,1)->6
  expect_equal(degrees[5 + 1], 4L)
  expect_equal(degrees[6 + 1], 4L)
  
  # Edge pixels (non-corner): (1,0)->1
  expect_equal(degrees[1 + 1], 3L)
})

test_that("invalid connectivity throws error", {
  expect_error(
    pixel_graph_stats_rust(4L, 3L, 5L),
    "`neighbors` must be 4 or 8"
  )
})
