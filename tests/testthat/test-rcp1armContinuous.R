test_that("rcp1armContinuous: formula approach returns correct structure", {
  result <- rcp1armContinuous(
    mu  = 0.5,
    mu0 = 0.1,
    sd  = 1,
    Nj  = c(10, 90),
    PI  = 0.5,
    approach = "formula"
  )

  expect_s3_class(result, "rcp1armContinuous")
  expect_named(result, c("approach", "nsim", "mu", "mu0", "sd", "Nj",
                         "PI", "Method1", "Method2"))
  expect_equal(result$approach, "formula")
  expect_null(result$nsim)
  expect_equal(result$mu,  0.5)
  expect_equal(result$mu0, 0.1)
  expect_equal(result$sd,  1)
  expect_equal(result$Nj,  c(10, 90))
  expect_equal(result$PI,  0.5)
})

test_that("rcp1armContinuous: formula RCP values are in [0, 1]", {
  result <- rcp1armContinuous(
    mu  = 0.5,
    mu0 = 0.1,
    sd  = 1,
    Nj  = c(10, 90),
    PI  = 0.5,
    approach = "formula"
  )

  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armContinuous: formula Method1 matches manual calculation", {
  mu  <- 0.5; mu0 <- 0.1; sd <- 1
  Nj  <- c(10, 90); PI <- 0.5
  N   <- sum(Nj); f1 <- Nj[1] / N
  delta  <- mu - mu0
  mean_d <- (1 - PI) * delta
  var_d  <- (1 - PI * f1)^2 * (sd^2 / Nj[1]) +
            (PI * (1 - f1))^2 * (sd^2 / (N - Nj[1]))
  expected_m1 <- pnorm(mean_d / sqrt(var_d))

  result <- rcp1armContinuous(
    mu = mu, mu0 = mu0, sd = sd, Nj = Nj,
    PI = PI, approach = "formula"
  )

  expect_equal(result$Method1, expected_m1, tolerance = 1e-10)
})

test_that("rcp1armContinuous: formula Method2 matches manual calculation", {
  mu <- 0.5; mu0 <- 0.1; sd <- 1
  Nj <- c(10, 90); PI <- 0.5
  delta <- mu - mu0
  expected_m2 <- prod(pnorm(delta * sqrt(Nj) / sd))

  result <- rcp1armContinuous(
    mu = mu, mu0 = mu0, sd = sd, Nj = Nj,
    PI = PI, approach = "formula"
  )

  expect_equal(result$Method2, expected_m2, tolerance = 1e-10)
})

test_that("rcp1armContinuous: simulation approach returns correct structure", {
  result <- rcp1armContinuous(
    mu  = 0.5,
    mu0 = 0.1,
    sd  = 1,
    Nj  = c(10, 90),
    PI  = 0.5,
    approach = "simulation",
    nsim = 1000,
    seed = 42
  )

  expect_s3_class(result, "rcp1armContinuous")
  expect_equal(result$approach, "simulation")
  expect_equal(result$nsim, 1000)
  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armContinuous: simulation is reproducible with same seed", {
  res1 <- rcp1armContinuous(
    mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(10, 90),
    approach = "simulation", nsim = 1000, seed = 1
  )
  res2 <- rcp1armContinuous(
    mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(10, 90),
    approach = "simulation", nsim = 1000, seed = 1
  )

  expect_equal(res1$Method1, res2$Method1)
  expect_equal(res1$Method2, res2$Method2)
})

test_that("rcp1armContinuous: formula and simulation agree within tolerance", {
  res_f <- rcp1armContinuous(
    mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(30, 70),
    PI = 0.5, approach = "formula"
  )
  res_s <- rcp1armContinuous(
    mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(30, 70),
    PI = 0.5, approach = "simulation", nsim = 50000, seed = 1
  )

  expect_equal(res_f$Method1, res_s$Method1, tolerance = 0.02)
  expect_equal(res_f$Method2, res_s$Method2, tolerance = 0.02)
})

test_that("rcp1armContinuous: J >= 3 regions work correctly", {
  result <- rcp1armContinuous(
    mu  = 0.5,
    mu0 = 0.1,
    sd  = 1,
    Nj  = c(10, 45, 45),
    PI  = 0.5,
    approach = "formula"
  )

  expect_s3_class(result, "rcp1armContinuous")
  expect_equal(length(result$Nj), 3)
  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
})

test_that("rcp1armContinuous: PI = 0 gives Method1 = 1", {
  result <- rcp1armContinuous(
    mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(10, 90),
    PI = 0, approach = "formula"
  )
  # PI = 0: condition is always (hat_mu1 - mu0) > 0, which has probability
  # Phi(delta * sqrt(N1) / sd) > 0.5 when delta > 0.
  # mean_d = delta, sd_d = sd / sqrt(N1), so Method1 = Phi(delta*sqrt(N1)/sd)
  expect_gt(result$Method1, 0.5)
})

test_that("rcp1armContinuous: input validation errors", {
  expect_error(
    rcp1armContinuous(mu = c(0.5, 0.5), mu0 = 0.1, sd = 1, Nj = c(10, 90)),
    "mu must be a single numeric value"
  )
  expect_error(
    rcp1armContinuous(mu = 0.5, mu0 = 0.1, sd = -1, Nj = c(10, 90)),
    "sd must be a single positive number"
  )
  expect_error(
    rcp1armContinuous(mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(10.5, 89.5)),
    "Nj must be a vector of positive integers"
  )
  expect_error(
    rcp1armContinuous(mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(10, 90), PI = 1.5),
    "PI must be a single numeric value in"
  )
  expect_error(
    rcp1armContinuous(mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(10, 90),
                      approach = "invalid"),
    'approach must be either'
  )
  expect_error(
    rcp1armContinuous(mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(10, 90),
                      approach = "simulation", nsim = -100),
    "nsim must be a single positive integer"
  )
  expect_error(
    rcp1armContinuous(mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(10, 90),
                      approach = "simulation", nsim = 1000, seed = -1),
    "seed must be a single non-negative integer"
  )
})

test_that("rcp1armContinuous: print method runs without error", {
  result <- rcp1armContinuous(
    mu = 0.5, mu0 = 0.1, sd = 1, Nj = c(10, 90),
    approach = "formula"
  )
  expect_output(print(result), "Regional Consistency Probability")
  expect_output(print(result), "Continuous")
  expect_output(print(result), "Method 1")
  expect_output(print(result), "Method 2")
})
