test_that("rcp1armMilestoneSurvival: formula approach returns correct structure", {
  result <- rcp1armMilestoneSurvival(
    lambda  = log(2) / 10,
    t_eval  = 8,
    S0      = exp(-log(2) * 8 / 5),
    Nj      = c(10, 90),
    t_a     = 3,
    t_f     = 10,
    PI      = 0.5,
    approach = "formula"
  )

  expect_s3_class(result, "rcp1armMilestoneSurvival")
  expect_named(result, c("approach", "formula_type", "nsim", "lambda",
                         "Nj", "t_a", "t_f", "tau", "lambda_dropout",
                         "PI", "eval_time", "S0", "S_est",
                         "Method1", "Method2"))
  expect_equal(result$approach, "formula")
  expect_null(result$nsim)
  expect_equal(result$lambda,    log(2) / 10)
  expect_equal(result$t_a,       3)
  expect_equal(result$t_f,       10)
  expect_equal(result$tau,       13)
  expect_true(is.na(result$lambda_dropout))
  expect_equal(result$eval_time, 8)
  expect_equal(result$S0,        exp(-log(2) * 8 / 5), tolerance = 1e-10)
})

test_that("rcp1armMilestoneSurvival: formula_type is closed-form when t_eval <= t_f", {
  result <- rcp1armMilestoneSurvival(
    lambda = log(2) / 10, t_eval = 8,
    S0 = exp(-log(2) * 8 / 5),
    Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
  )
  # t_eval = 8 <= t_f = 10
  expect_equal(result$formula_type, "closed-form")
})

test_that("rcp1armMilestoneSurvival: formula_type is numerical-integration when t_eval > t_f", {
  result <- rcp1armMilestoneSurvival(
    lambda = log(2) / 10, t_eval = 12,
    S0 = exp(-log(2) * 12 / 5),
    Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
  )
  # t_eval = 12 > t_f = 10
  expect_equal(result$formula_type, "numerical-integration")
})

test_that("rcp1armMilestoneSurvival: S_est matches exponential model", {
  lam    <- log(2) / 10
  t_eval <- 8
  result <- rcp1armMilestoneSurvival(
    lambda = lam, t_eval = t_eval,
    S0 = exp(-log(2) * 8 / 5),
    Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
  )
  expect_equal(result$S_est, exp(-lam * t_eval), tolerance = 1e-10)
})

test_that("rcp1armMilestoneSurvival: formula RCP values are in [0, 1]", {
  result <- rcp1armMilestoneSurvival(
    lambda = log(2) / 10, t_eval = 8,
    S0 = exp(-log(2) * 8 / 5),
    Nj = c(10, 90), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )

  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armMilestoneSurvival: simulation approach returns correct structure", {
  result <- rcp1armMilestoneSurvival(
    lambda  = log(2) / 10,
    t_eval  = 8,
    S0      = exp(-log(2) * 8 / 5),
    Nj      = c(10, 90),
    t_a     = 3,
    t_f     = 10,
    PI      = 0.5,
    approach = "simulation",
    nsim     = 1000,
    seed     = 42
  )

  expect_s3_class(result, "rcp1armMilestoneSurvival")
  expect_equal(result$approach, "simulation")
  expect_null(result$formula_type)
  expect_equal(result$nsim, 1000)
  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armMilestoneSurvival: simulation is reproducible with same seed", {
  res1 <- rcp1armMilestoneSurvival(
    lambda = log(2) / 10, t_eval = 8,
    S0 = exp(-log(2) * 8 / 5),
    Nj = c(10, 90), t_a = 3, t_f = 10,
    approach = "simulation", nsim = 1000, seed = 1
  )
  res2 <- rcp1armMilestoneSurvival(
    lambda = log(2) / 10, t_eval = 8,
    S0 = exp(-log(2) * 8 / 5),
    Nj = c(10, 90), t_a = 3, t_f = 10,
    approach = "simulation", nsim = 1000, seed = 1
  )

  expect_equal(res1$Method1, res2$Method1)
  expect_equal(res1$Method2, res2$Method2)
})

test_that("rcp1armMilestoneSurvival: formula and simulation agree within tolerance", {
  res_f <- rcp1armMilestoneSurvival(
    lambda = log(2) / 10, t_eval = 8,
    S0 = exp(-log(2) * 8 / 5),
    Nj = c(30, 70), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )
  res_s <- rcp1armMilestoneSurvival(
    lambda = log(2) / 10, t_eval = 8,
    S0 = exp(-log(2) * 8 / 5),
    Nj = c(30, 70), t_a = 3, t_f = 10,
    PI = 0.5, approach = "simulation", nsim = 50000, seed = 1
  )

  expect_equal(res_f$Method1, res_s$Method1, tolerance = 0.03)
  expect_equal(res_f$Method2, res_s$Method2, tolerance = 0.03)
})

test_that("rcp1armMilestoneSurvival: J >= 3 regions work correctly", {
  result <- rcp1armMilestoneSurvival(
    lambda = log(2) / 10, t_eval = 8,
    S0 = exp(-log(2) * 8 / 5),
    Nj = c(10, 45, 45), t_a = 3, t_f = 10,
    PI = 0.5, approach = "formula"
  )

  expect_s3_class(result, "rcp1armMilestoneSurvival")
  expect_equal(length(result$Nj), 3)
  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
})

test_that("rcp1armMilestoneSurvival: closed-form equals numerical-integration at boundary", {
  # When t_eval is just at t_f, both branches should give very similar results
  lam <- log(2) / 10
  S0  <- exp(-log(2) * 10 / 5)
  res_cf <- rcp1armMilestoneSurvival(
    lambda = lam, t_eval = 10, S0 = S0,
    Nj = c(30, 70), t_a = 3, t_f = 10, approach = "formula"
  )
  # t_eval slightly above t_f triggers numerical-integration
  res_ni <- rcp1armMilestoneSurvival(
    lambda = lam, t_eval = 10 + 1e-6, S0 = S0,
    Nj = c(30, 70), t_a = 3, t_f = 10, approach = "formula"
  )

  expect_equal(res_cf$formula_type, "closed-form")
  expect_equal(res_ni$formula_type, "numerical-integration")
  expect_equal(res_cf$Method1, res_ni$Method1, tolerance = 1e-4)
  expect_equal(res_cf$Method2, res_ni$Method2, tolerance = 1e-4)
})

test_that("rcp1armMilestoneSurvival: input validation errors", {
  expect_error(
    rcp1armMilestoneSurvival(lambda = -0.1, t_eval = 8,
                             S0 = 0.5, Nj = c(10, 90), t_a = 3, t_f = 10),
    "lambda must be a single positive number"
  )
  expect_error(
    rcp1armMilestoneSurvival(lambda = log(2) / 10, t_eval = 0,
                             S0 = 0.5, Nj = c(10, 90), t_a = 3, t_f = 10),
    "t_eval must be a single positive number"
  )
  expect_error(
    rcp1armMilestoneSurvival(lambda = log(2) / 10, t_eval = 8,
                             S0 = 0, Nj = c(10, 90), t_a = 3, t_f = 10),
    "S0 must be a single numeric value in"
  )
  expect_error(
    rcp1armMilestoneSurvival(lambda = log(2) / 10, t_eval = 8,
                             S0 = 1.1, Nj = c(10, 90), t_a = 3, t_f = 10),
    "S0 must be a single numeric value in"
  )
  expect_error(
    rcp1armMilestoneSurvival(lambda = log(2) / 10, t_eval = 8,
                             S0 = 0.5, Nj = c(10, 90), t_a = -1, t_f = 10),
    "t_a must be a single positive number"
  )
  expect_error(
    rcp1armMilestoneSurvival(lambda = log(2) / 10, t_eval = 8,
                             S0 = 0.5, Nj = c(10, 90), t_a = 3, t_f = 10,
                             approach = "invalid"),
    'approach must be either'
  )
})

test_that("rcp1armMilestoneSurvival: print method runs without error", {
  result <- rcp1armMilestoneSurvival(
    lambda = log(2) / 10, t_eval = 8,
    S0 = exp(-log(2) * 8 / 5),
    Nj = c(10, 90), t_a = 3, t_f = 10, approach = "formula"
  )
  expect_output(print(result), "Regional Consistency Probability")
  expect_output(print(result), "Milestone Survival")
  expect_output(print(result), "Method 1")
  expect_output(print(result), "Method 2")
})
