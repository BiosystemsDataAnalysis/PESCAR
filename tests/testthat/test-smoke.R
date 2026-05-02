test_that("PESCAR package smoke test", {
  expect_true(requireNamespace("PESCAR", quietly = TRUE))
  expect_gt(length(getNamespaceExports("PESCAR")), 0)
})
