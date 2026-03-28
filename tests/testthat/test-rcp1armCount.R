test_that("rcp1armCount: formula approach returns correct structure", {
  result <- rcp1armCount(
    lambda     = 2,
    lambda0    = 3,
    dispersion = 1,
    Nj         = c(10, 90),
    PI         = 0.5,
    approach   = "formula"
  )

  expect_s3_class(result, "rcp1armCount")
  expect_named(result, c("approach", "nsim", "lambda", "lambda0",
                         "dispersion", "Nj", "PI",
                         "Method1_logRR", "Method1_linearRR", "Method2"))
  expect_equal(result$approach, "formula")
  expect_null(result$nsim)
  expect_equal(result$lambda,     2)
  expect_equal(result$lambda0,    3)
  expect_equal(result$dispersion, 1)
  expect_equal(result$Nj,         c(10, 90))
  expect_equal(result$PI,         0.5)
})

test_that("rcp1armCount: formula RCP values are in [0, 1]", {
  result <- rcp1armCount(
    lambda = 2, lambda0 = 3, dispersion = 1,
    Nj = c(10, 90), PI = 0.5, approach = "formula"
  )

  expect_gte(result$Method1_logRR,    0)
  expect_lte(result$Method1_logRR,    1)
  expect_gte(result$Method1_linearRR, 0)
  expect_lte(result$Method1_linearRR, 1)
  expect_gte(result$Method2,          0)
  expect_lte(result$Method2,          1)
})

test_that("rcp1armCount: formula Method2 matches manual calculation", {
  lambda <- 2; lambda0 <- 3; dispersion <- 1; Nj <- c(10, 90)
  J <- length(Nj)
  probs <- numeric(J)
  for (j in seq_len(J)) {
    mu_j   <- Nj[j] * lambda
    size_j <- Nj[j] * dispersion
    y_thr  <- floor(Nj[j] * lambda0 - 1e-10)
    probs[j] <- if (y_thr < 0) 0 else pnbinom(y_thr, mu = mu_j, size = size_j)
  }
  expected_m2 <- prod(probs)

  result <- rcp1armCount(
    lambda = lambda, lambda0 = lambda0, dispersion = dispersion,
    Nj = Nj, PI = 0.5, approach = "formula"
  )

  expect_equal(result$Method2, expected_m2, tolerance = 1e-10)
})

test_that("rcp1armCount: simulation approach returns correct structure", {
  result <- rcp1armCount(
    lambda     = 2,
    lambda0    = 3,
    dispersion = 1,
    Nj         = c(10, 90),
    PI         = 0.5,
    approach   = "simulation",
    nsim       = 1000,
    seed       = 42
  )

  expect_s3_class(result, "rcp1armCount")
  expect_equal(result$approach, "simulation")
  expect_equal(result$nsim, 1000)
  expect_gte(result$Method1_logRR,    0)
  expect_lte(result$Method1_logRR,    1)
  expect_gte(result$Method1_linearRR, 0)
  expect_lte(result$Method1_linearRR, 1)
  expect_gte(result$Method2,          0)
  expect_lte(result$Method2,          1)
})

test_that("rcp1armCount: simulation is reproducible with same seed", {
  res1 <- rcp1armCount(
    lambda = 2, lambda0 = 3, dispersion = 1, Nj = c(10, 90),
    approach = "simulation", nsim = 1000, seed = 1
  )
  res2 <- rcp1armCount(
    lambda = 2, lambda0 = 3, dispersion = 1, Nj = c(10, 90),
    approach = "simulation", nsim = 1000, seed = 1
  )

  expect_equal(res1$Method1_logRR,    res2$Method1_logRR)
  expect_equal(res1$Method1_linearRR, res2$Method1_linearRR)
  expect_equal(res1$Method2,          res2$Method2)
})

test_that("rcp1armCount: formula and simulation agree within tolerance", {
  res_f <- rcp1armCount(
    lambda = 2, lambda0 = 3, dispersion = 1,
    Nj = c(30, 70), PI = 0.5, approach = "formula"
  )
  res_s <- rcp1armCount(
    lambda = 2, lambda0 = 3, dispersion = 1,
    Nj = c(30, 70), PI = 0.5, approach = "simulation",
    nsim = 50000, seed = 1
  )

  expect_equal(res_f$Method1_logRR,    res_s$Method1_logRR,    tolerance = 0.03)
  expect_equal(res_f$Method1_linearRR, res_s$Method1_linearRR, tolerance = 0.03)
  expect_equal(res_f$Method2,          res_s$Method2,          tolerance = 0.03)
})

test_that("rcp1armCount: J >= 3 regions work correctly", {
  result <- rcp1armCount(
    lambda = 2, lambda0 = 3, dispersion = 1,
    Nj = c(10, 45, 45), PI = 0.5, approach = "formula"
  )

  expect_s3_class(result, "rcp1armCount")
  expect_equal(length(result$Nj), 3)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armCount: lambda >= lambda0 gives low Method2", {
  # When lambda >= lambda0, RR >= 1 so benefit probability should be near 0
  result <- rcp1armCount(
    lambda = 5, lambda0 = 3, dispersion = 1,
    Nj = c(30, 70), PI = 0.5, approach = "formula"
  )
  expect_lt(result$Method2, 0.5)
})

test_that("rcp1armCount: input validation errors", {
  expect_error(
    rcp1armCount(lambda = -1, lambda0 = 3, dispersion = 1, Nj = c(10, 90)),
    "lambda must be a single positive number"
  )
  expect_error(
    rcp1armCount(lambda = 2, lambda0 = 0, dispersion = 1, Nj = c(10, 90)),
    "lambda0 must be a single positive number"
  )
  expect_error(
    rcp1armCount(lambda = 2, lambda0 = 3, dispersion = -1, Nj = c(10, 90)),
    "dispersion must be a single positive number"
  )
  expect_error(
    rcp1armCount(lambda = 2, lambda0 = 3, dispersion = 1, Nj = c(0, 90)),
    "Nj must be a vector of positive integers"
  )
  expect_error(
    rcp1armCount(lambda = 2, lambda0 = 3, dispersion = 1, Nj = c(10, 90),
                 PI = 2),
    "PI must be a single numeric value in"
  )
  expect_error(
    rcp1armCount(lambda = 2, lambda0 = 3, dispersion = 1, Nj = c(10, 90),
                 approach = "invalid"),
    'approach must be either'
  )
})

test_that("rcp1armCount: print method runs without error", {
  result <- rcp1armCount(
    lambda = 2, lambda0 = 3, dispersion = 1,
    Nj = c(10, 90), approach = "formula"
  )
  expect_output(print(result), "Regional Consistency Probability")
  expect_output(print(result), "Count")
  expect_output(print(result), "Log-RR")
  expect_output(print(result), "Linear-RR")
})
