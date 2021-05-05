test_that("cusumr works", {
  df <- c(0,0,0,0,1,1,1,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

  expect_type(cusumr(df), "list")
  expect_true(sum(is.na(cusumr(df)$cusum_score)) == 0)
})
