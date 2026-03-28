test_that("rcp1armHazardRatio: formula approach returns correct structure", {
  result <- rcp1armHazardRatio(
    lambda  = log(2) / 10,
    lambda0 = log(2) / 5,
    Nj      = c(10, 90),
    t_a     = 3,
    t_f     = 10,
    PI      = 0.5,
    approach = "formula"
  )

  expect_s3_class(result, "rcp1armHazardRatio")
  expect_named(result, c("approach", "nsim", "lambda", "lambda0", "Nj",
                         "t_a", "t_f", "tau", "lambda_dropout", "PI",
                         "Method1_logHR", "Method1_linearHR", "Method2"))
  expect_equal(result$approach, "formula")
  expect_null(result$nsim)
  expect_equal(result$lambda,  log(2) / 10)
  expect_equal(result$lambda0, log(2) / 5)
  expect_equal(result$Nj,      c(10, 90))
  expect_equal(result$t_a,     3)
  expect_equal(result$t_f,     10)
  expect_equal(result$tau,     13)
  expect_true(is.na(result$lambda_dropout))
  expect_equal(result$PI,      0.5)
})

test_that("rcp1armHazardRatio: formula RCP values are in [0, 1]", {
  result <- rcp1armHazardRatio(
    lambda = log(2) / 10, lambda0 = log(2) / 5,
    Nj = c(10, 90), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )

  expect_gte(result$Method1_logHR,    0)
  expect_lte(result$Method1_logHR,    1)
  expect_gte(result$Method1_linearHR, 0)
  expect_lte(result$Method1_linearHR, 1)
  expect_gte(result$Method2,          0)
  expect_lte(result$Method2,          1)
})

test_that("rcp1armHazardRatio: tau is correctly computed", {
  result <- rcp1armHazardRatio(
    lambda = log(2) / 10, lambda0 = log(2) / 5,
    Nj = c(10, 90), t_a = 4, t_f = 8, approach = "formula"
  )
  expect_equal(result$tau, 12)
})

test_that("rcp1armHazardRatio: lambda_dropout is stored correctly", {
  result <- rcp1armHazardRatio(
    lambda = log(2) / 10, lambda0 = log(2) / 5,
    Nj = c(10, 90), t_a = 3, t_f = 10,
    lambda_dropout = 0.05, approach = "formula"
  )
  expect_equal(result$lambda_dropout, 0.05)
})

test_that("rcp1armHazardRatio: simulation approach returns correct structure", {
  result <- rcp1armHazardRatio(
    lambda  = log(2) / 10,
    lambda0 = log(2) / 5,
    Nj      = c(10, 90),
    t_a     = 3,
    t_f     = 10,
    PI      = 0.5,
    approach = "simulation",
    nsim     = 1000,
    seed     = 42
  )

  expect_s3_class(result, "rcp1armHazardRatio")
  expect_equal(result$approach, "simulation")
  expect_equal(result$nsim, 1000)
  expect_gte(result$Method1_logHR,    0)
  expect_lte(result$Method1_logHR,    1)
  expect_gte(result$Method1_linearHR, 0)
  expect_lte(result$Method1_linearHR, 1)
  expect_gte(result$Method2,          0)
  expect_lte(result$Method2,          1)
})

test_that("rcp1armHazardRatio: simulation is reproducible with same seed", {
  res1 <- rcp1armHazardRatio(
    lambda = log(2) / 10, lambda0 = log(2) / 5,
    Nj = c(10, 90), t_a = 3, t_f = 10,
    approach = "simulation", nsim = 1000, seed = 1
  )
  res2 <- rcp1armHazardRatio(
    lambda = log(2) / 10, lambda0 = log(2) / 5,
    Nj = c(10, 90), t_a = 3, t_f = 10,
    approach = "simulation", nsim = 1000, seed = 1
  )

  expect_equal(res1$Method1_logHR,    res2$Method1_logHR)
  expect_equal(res1$Method1_linearHR, res2$Method1_linearHR)
  expect_equal(res1$Method2,          res2$Method2)
})

test_that("rcp1armHazardRatio: formula and simulation agree within tolerance", {
  res_f <- rcp1armHazardRatio(
    lambda = log(2) / 10, lambda0 = log(2) / 5,
    Nj = c(30, 70), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )
  res_s <- rcp1armHazardRatio(
    lambda = log(2) / 10, lambda0 = log(2) / 5,
    Nj = c(30, 70), t_a = 3, t_f = 10,
    PI = 0.5, approach = "simulation", nsim = 50000, seed = 1
  )

  expect_equal(res_f$Method1_logHR,    res_s$Method1_logHR,    tolerance = 0.03)
  expect_equal(res_f$Method1_linearHR, res_s$Method1_linearHR, tolerance = 0.03)
  expect_equal(res_f$Method2,          res_s$Method2,          tolerance = 0.03)
})

test_that("rcp1armHazardRatio: J >= 3 regions work correctly", {
  result <- rcp1armHazardRatio(
    lambda = log(2) / 10, lambda0 = log(2) / 5,
    Nj = c(10, 45, 45), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )

  expect_s3_class(result, "rcp1armHazardRatio")
  expect_equal(length(result$Nj), 3)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armHazardRatio: HR >= 1 (no benefit) gives low Method2", {
  # lambda >= lambda0 means no benefit: Method2 should be near 0
  result <- rcp1armHazardRatio(
    lambda = log(2) / 3, lambda0 = log(2) / 5,
    Nj = c(30, 70), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )
  expect_lt(result$Method2, 0.5)
})

test_that("rcp1armHazardRatio: input validation errors", {
  expect_error(
    rcp1armHazardRatio(lambda = -0.1, lambda0 = log(2) / 5,
                       Nj = c(10, 90), t_a = 3, t_f = 10),
    "lambda must be a single positive number"
  )
  expect_error(
    rcp1armHazardRatio(lambda = log(2) / 10, lambda0 = 0,
                       Nj = c(10, 90), t_a = 3, t_f = 10),
    "lambda0 must be a single positive number"
  )
  expect_error(
    rcp1armHazardRatio(lambda = log(2) / 10, lambda0 = log(2) / 5,
                       Nj = c(10, 90), t_a = -1, t_f = 10),
    "t_a must be a single positive number"
  )
  expect_error(
    rcp1armHazardRatio(lambda = log(2) / 10, lambda0 = log(2) / 5,
                       Nj = c(10, 90), t_a = 3, t_f = -1),
    "t_f must be a single positive number"
  )
  expect_error(
    rcp1armHazardRatio(lambda = log(2) / 10, lambda0 = log(2) / 5,
                       Nj = c(10, 90), t_a = 3, t_f = 10,
                       lambda_dropout = -0.05),
    "lambda_dropout must be a single positive number or NULL"
  )
  expect_error(
    rcp1armHazardRatio(lambda = log(2) / 10, lambda0 = log(2) / 5,
                       Nj = c(10, 90), t_a = 3, t_f = 10, PI = 1.5),
    "PI must be a single numeric value in"
  )
  expect_error(
    rcp1armHazardRatio(lambda = log(2) / 10, lambda0 = log(2) / 5,
                       Nj = c(10, 90), t_a = 3, t_f = 10,
                       approach = "invalid"),
    'approach must be either'
  )
})

test_that("rcp1armHazardRatio: print method runs without error", {
  result <- rcp1armHazardRatio(
    lambda = log(2) / 10, lambda0 = log(2) / 5,
    Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
  )
  expect_output(print(result), "Regional Consistency Probability")
  expect_output(print(result), "Hazard Ratio")
  expect_output(print(result), "Log-HR")
  expect_output(print(result), "Linear-HR")
})
