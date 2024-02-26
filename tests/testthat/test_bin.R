context("binary_ss_checks")

test_that("Test 4.1", {
  ss <- c(905, 4501, 4252, 4501)
  test <- pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpnormal = c(-5,2.5), oeciwidth = 1, simobs = 1000)
  test_ss <- test$results_table[,1]
  names(test_ss) <- c()
  expect_equal(test_ss,ss)
})

test_that("Test 4.2", {
  ss <- c(905, 1641, 4252, 4252)
  test <- pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpbeta = c(0.5,0.5), oeciwidth = 1, simobs = 1000)
  test_ss <- test$results_table[,1]
  names(test_ss) <- c()
  expect_equal(test_ss,ss)
})

test_that("Test 4.3", {
  ss <- c(905, 16175, 4252, 16175)
  test <- pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, lpcstat = -4.7, oeciwidth = 1, seed = 1234, simobs = 10000)
  test_ss <- test$results_table[,1]
  names(test_ss) <- c()
  expect_equal(test_ss,ss)
})


