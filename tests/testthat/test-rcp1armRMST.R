test_that("rcp1armRMST: formula approach returns correct structure", {
  lam0    <- log(2) / 5
  tstar   <- 8
  mu0_val <- (1 - exp(-lam0 * tstar)) / lam0

  result <- rcp1armRMST(
    lambda   = log(2) / 10,
    tau_star = tstar,
    mu0      = mu0_val,
    Nj       = c(10, 90),
    t_a      = 3,
    t_f      = 10,
    PI       = 0.5,
    approach = "formula"
  )

  expect_s3_class(result, "rcp1armRMST")
  expect_named(result, c("approach", "formula_type", "nsim", "lambda",
                         "Nj", "t_a", "t_f", "tau", "lambda_dropout",
                         "PI", "tau_star", "mu0", "mu_est",
                         "Method1", "Method2"))
  expect_equal(result$approach, "formula")
  expect_null(result$nsim)
  expect_equal(result$lambda,   log(2) / 10)
  expect_equal(result$t_a,      3)
  expect_equal(result$t_f,      10)
  expect_equal(result$tau,      13)
  expect_true(is.na(result$lambda_dropout))
  expect_equal(result$PI,       0.5)
  expect_equal(result$tau_star, tstar)
  expect_equal(result$mu0,      mu0_val, tolerance = 1e-10)
})

test_that("rcp1armRMST: formula_type is closed-form when tau_star <= t_f", {
  lam0  <- log(2) / 5; tstar <- 8
  result <- rcp1armRMST(
    lambda = log(2) / 10, tau_star = tstar,
    mu0 = (1 - exp(-lam0 * tstar)) / lam0,
    Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
  )
  # tau_star = 8 <= t_f = 10
  expect_equal(result$formula_type, "closed-form")
})

test_that("rcp1armRMST: formula_type is numerical-integration when tau_star > t_f", {
  lam0  <- log(2) / 5; tstar <- 12
  result <- rcp1armRMST(
    lambda = log(2) / 10, tau_star = tstar,
    mu0 = (1 - exp(-lam0 * tstar)) / lam0,
    Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
  )
  # tau_star = 12 > t_f = 10
  expect_equal(result$formula_type, "numerical-integration")
})

test_that("rcp1armRMST: mu_est matches exponential model", {
  lam   <- log(2) / 10; tstar <- 8
  lam0  <- log(2) / 5
  result <- rcp1armRMST(
    lambda = lam, tau_star = tstar,
    mu0 = (1 - exp(-lam0 * tstar)) / lam0,
    Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
  )
  expect_equal(result$mu_est, (1 - exp(-lam * tstar)) / lam, tolerance = 1e-10)
})

test_that("rcp1armRMST: formula RCP values are in [0, 1]", {
  lam0  <- log(2) / 5; tstar <- 8
  result <- rcp1armRMST(
    lambda = log(2) / 10, tau_star = tstar,
    mu0 = (1 - exp(-lam0 * tstar)) / lam0,
    Nj = c(10, 90), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )

  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armRMST: tau_star > t_a + t_f triggers error", {
  lam0 <- log(2) / 5
  expect_error(
    rcp1armRMST(
      lambda = log(2) / 10, tau_star = 15,
      mu0 = (1 - exp(-lam0 * 15)) / lam0,
      Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
    ),
    "tau_star must not exceed"
  )
})

test_that("rcp1armRMST: simulation approach returns correct structure", {
  lam0  <- log(2) / 5; tstar <- 8
  result <- rcp1armRMST(
    lambda   = log(2) / 10,
    tau_star = tstar,
    mu0      = (1 - exp(-lam0 * tstar)) / lam0,
    Nj       = c(10, 90),
    t_a      = 3,
    t_f      = 10,
    PI       = 0.5,
    approach = "simulation",
    nsim     = 1000,
    seed     = 42
  )

  expect_s3_class(result, "rcp1armRMST")
  expect_equal(result$approach, "simulation")
  expect_null(result$formula_type)
  expect_equal(result$nsim, 1000)
  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armRMST: simulation is reproducible with same seed", {
  lam0  <- log(2) / 5; tstar <- 8
  mu0_v <- (1 - exp(-lam0 * tstar)) / lam0
  res1 <- rcp1armRMST(
    lambda = log(2) / 10, tau_star = tstar, mu0 = mu0_v,
    Nj = c(10, 90), t_a = 3, t_f = 10,
    approach = "simulation", nsim = 1000, seed = 1
  )
  res2 <- rcp1armRMST(
    lambda = log(2) / 10, tau_star = tstar, mu0 = mu0_v,
    Nj = c(10, 90), t_a = 3, t_f = 10,
    approach = "simulation", nsim = 1000, seed = 1
  )

  expect_equal(res1$Method1, res2$Method1)
  expect_equal(res1$Method2, res2$Method2)
})

test_that("rcp1armRMST: Nj < 3 triggers error for simulation", {
  lam0  <- log(2) / 5; tstar <- 8
  expect_error(
    rcp1armRMST(
      lambda = log(2) / 10, tau_star = tstar,
      mu0 = (1 - exp(-lam0 * tstar)) / lam0,
      Nj = c(2, 90), t_a = 3, t_f = 10,
      approach = "simulation", nsim = 100, seed = 1
    ),
    "Nj must have all elements >= 3"
  )
})

test_that("rcp1armRMST: formula and simulation agree within tolerance", {
  lam0  <- log(2) / 5; tstar <- 8
  mu0_v <- (1 - exp(-lam0 * tstar)) / lam0
  res_f <- rcp1armRMST(
    lambda = log(2) / 10, tau_star = tstar, mu0 = mu0_v,
    Nj = c(30, 70), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )
  res_s <- rcp1armRMST(
    lambda = log(2) / 10, tau_star = tstar, mu0 = mu0_v,
    Nj = c(30, 70), t_a = 3, t_f = 10,
    PI = 0.5, approach = "simulation", nsim = 50000, seed = 1
  )

  expect_equal(res_f$Method1, res_s$Method1, tolerance = 0.03)
  expect_equal(res_f$Method2, res_s$Method2, tolerance = 0.03)
})

test_that("rcp1armRMST: closed-form equals numerical-integration at boundary", {
  lam   <- log(2) / 10; lam0 <- log(2) / 5
  # tau_star exactly at t_f
  res_cf <- rcp1armRMST(
    lambda = lam, tau_star = 10,
    mu0 = (1 - exp(-lam0 * 10)) / lam0,
    Nj = c(30, 70), t_a = 3, t_f = 10, approach = "formula"
  )
  # tau_star slightly above t_f
  res_ni <- rcp1armRMST(
    lambda = lam, tau_star = 10 + 1e-6,
    mu0 = (1 - exp(-lam0 * (10 + 1e-6))) / lam0,
    Nj = c(30, 70), t_a = 3, t_f = 10, approach = "formula"
  )

  expect_equal(res_cf$formula_type, "closed-form")
  expect_equal(res_ni$formula_type, "numerical-integration")
  expect_equal(res_cf$Method1, res_ni$Method1, tolerance = 1e-4)
  expect_equal(res_cf$Method2, res_ni$Method2, tolerance = 1e-4)
})

test_that("rcp1armRMST: J >= 3 regions work correctly", {
  lam0  <- log(2) / 5; tstar <- 8
  result <- rcp1armRMST(
    lambda = log(2) / 10, tau_star = tstar,
    mu0 = (1 - exp(-lam0 * tstar)) / lam0,
    Nj = c(10, 45, 45), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )

  expect_s3_class(result, "rcp1armRMST")
  expect_equal(length(result$Nj), 3)
  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
})

test_that("rcp1armRMST: input validation errors", {
  lam0  <- log(2) / 5; tstar <- 8
  mu0_v <- (1 - exp(-lam0 * tstar)) / lam0

  expect_error(
    rcp1armRMST(lambda = -0.1, tau_star = tstar, mu0 = mu0_v,
                Nj = c(10, 90), t_a = 3, t_f = 10),
    "lambda must be a single positive number"
  )
  expect_error(
    rcp1armRMST(lambda = log(2) / 10, tau_star = 0, mu0 = mu0_v,
                Nj = c(10, 90), t_a = 3, t_f = 10),
    "tau_star must be a single positive number"
  )
  expect_error(
    rcp1armRMST(lambda = log(2) / 10, tau_star = tstar, mu0 = -1,
                Nj = c(10, 90), t_a = 3, t_f = 10),
    "mu0 must be a single positive number"
  )
  expect_error(
    rcp1armRMST(lambda = log(2) / 10, tau_star = tstar, mu0 = mu0_v,
                Nj = c(10, 90), t_a = -1, t_f = 10),
    "t_a must be a single positive number"
  )
  expect_error(
    rcp1armRMST(lambda = log(2) / 10, tau_star = tstar, mu0 = mu0_v,
                Nj = c(10, 90), t_a = 3, t_f = 10,
                approach = "invalid"),
    'approach must be either'
  )
})

test_that("rcp1armRMST: print method runs without error", {
  lam0  <- log(2) / 5; tstar <- 8
  result <- rcp1armRMST(
    lambda = log(2) / 10, tau_star = tstar,
    mu0 = (1 - exp(-lam0 * tstar)) / lam0,
    Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
  )
  expect_output(print(result), "Regional Consistency Probability")
  expect_output(print(result), "RMST")
  expect_output(print(result), "Method 1")
  expect_output(print(result), "Method 2")
  expect_output(print(result), "Closed-Form")
})
