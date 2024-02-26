context("syntax_checks")

# general checks
#======================================
test_that("Test 1.1.1", {
  errormessage <- 'argument "type" is missing, with no default'
  expect_error(pmvalsampsize(prevalence = 0.018, cstatistic = 0.8, lpnormal = c(-5,2.5), oeciwidth = 1), errormessage)
})

test_that("Test 1.1.2", {
  errormessage <- 'prevalence must be specified for binary sample size'
  expect_error(pmvalsampsize(type = "b", cstatistic = 0.8, lpnormal = c(-5,2.5)), errormessage)
})

test_that("Test 1.1.3", {
  errormessage <- 'cstatistic must be specified for binary outcome models'
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, lpnormal = c(-5,2.5)), errormessage)
})


# checks for binary sample size
#======================================
test_that("Test 1.3.1", {
  errormessage <- "An LP distribution must be specified"
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8), errormessage)
})

test_that("Test 1.3.2", {
  errormessage <- "Only one LP distribution option can be specified"
  expect_error(pmvalsampsize(type = "b", prevalence = 0.018, lpnormal = c(-5,2.5), lpbeta = c(0.5,0.5), cstatistic = 0.8), errormessage)
})




