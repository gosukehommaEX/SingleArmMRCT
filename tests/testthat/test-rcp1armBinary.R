test_that("rcp1armBinary: formula approach returns correct structure", {
  result <- rcp1armBinary(
    p  = 0.5,
    p0 = 0.2,
    Nj = c(10, 90),
    PI = 0.5,
    approach = "formula"
  )

  expect_s3_class(result, "rcp1armBinary")
  expect_named(result, c("approach", "nsim", "p", "p0", "Nj", "PI",
                         "Method1", "Method2"))
  expect_equal(result$approach, "formula")
  expect_null(result$nsim)
  expect_equal(result$p,  0.5)
  expect_equal(result$p0, 0.2)
  expect_equal(result$Nj, c(10, 90))
  expect_equal(result$PI, 0.5)
})

test_that("rcp1armBinary: formula RCP values are in [0, 1]", {
  result <- rcp1armBinary(
    p = 0.5, p0 = 0.2, Nj = c(10, 90),
    PI = 0.5, approach = "formula"
  )

  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armBinary: formula Method2 matches manual calculation", {
  p <- 0.5; p0 <- 0.2; Nj <- c(10, 90)
  y_min_vec <- ceiling(Nj * p0 + 1e-10)
  expected_m2 <- prod(ifelse(
    y_min_vec > Nj, 0,
    pbinom(y_min_vec - 1, Nj, p, lower.tail = FALSE)
  ))

  result <- rcp1armBinary(
    p = p, p0 = p0, Nj = Nj, PI = 0.5, approach = "formula"
  )

  expect_equal(result$Method2, expected_m2, tolerance = 1e-10)
})

test_that("rcp1armBinary: simulation approach returns correct structure", {
  result <- rcp1armBinary(
    p  = 0.5,
    p0 = 0.2,
    Nj = c(10, 90),
    PI = 0.5,
    approach = "simulation",
    nsim = 1000,
    seed = 42
  )

  expect_s3_class(result, "rcp1armBinary")
  expect_equal(result$approach, "simulation")
  expect_equal(result$nsim, 1000)
  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
  expect_gte(result$Method2, 0)
  expect_lte(result$Method2, 1)
})

test_that("rcp1armBinary: simulation is reproducible with same seed", {
  res1 <- rcp1armBinary(
    p = 0.5, p0 = 0.2, Nj = c(10, 90),
    approach = "simulation", nsim = 1000, seed = 1
  )
  res2 <- rcp1armBinary(
    p = 0.5, p0 = 0.2, Nj = c(10, 90),
    approach = "simulation", nsim = 1000, seed = 1
  )

  expect_equal(res1$Method1, res2$Method1)
  expect_equal(res1$Method2, res2$Method2)
})

test_that("rcp1armBinary: formula and simulation agree within tolerance", {
  res_f <- rcp1armBinary(
    p = 0.5, p0 = 0.2, Nj = c(30, 70),
    PI = 0.5, approach = "formula"
  )
  res_s <- rcp1armBinary(
    p = 0.5, p0 = 0.2, Nj = c(30, 70),
    PI = 0.5, approach = "simulation", nsim = 50000, seed = 1
  )

  expect_equal(res_f$Method1, res_s$Method1, tolerance = 0.02)
  expect_equal(res_f$Method2, res_s$Method2, tolerance = 0.02)
})

test_that("rcp1armBinary: J >= 3 regions work correctly", {
  result <- rcp1armBinary(
    p  = 0.5,
    p0 = 0.2,
    Nj = c(10, 45, 45),
    PI = 0.5,
    approach = "formula"
  )

  expect_s3_class(result, "rcp1armBinary")
  expect_equal(length(result$Nj), 3)
  expect_gte(result$Method1, 0)
  expect_lte(result$Method1, 1)
})

test_that("rcp1armBinary: p0 = 0 is accepted", {
  result <- rcp1armBinary(
    p = 0.5, p0 = 0, Nj = c(10, 90),
    PI = 0.5, approach = "formula"
  )
  # When p0 = 0, the condition hat_pj > 0 is equivalent to Yj >= 1.
  # Pr(Yj >= 1) = 1 - (1 - p)^Nj < 1 strictly, since Yj = 0 has positive
  # probability. Method2 = prod_j [1 - (1-p)^Nj].
  expected_m2 <- prod(1 - (1 - 0.5)^c(10, 90))
  expect_equal(result$Method2, expected_m2, tolerance = 1e-10)
})

test_that("rcp1armBinary: input validation errors", {
  expect_error(
    rcp1armBinary(p = 0, p0 = 0.2, Nj = c(10, 90)),
    "p must be a single numeric value in"
  )
  expect_error(
    rcp1armBinary(p = 1, p0 = 0.2, Nj = c(10, 90)),
    "p must be a single numeric value in"
  )
  expect_error(
    rcp1armBinary(p = 0.5, p0 = 1, Nj = c(10, 90)),
    "p0 must be a single numeric value in"
  )
  expect_error(
    rcp1armBinary(p = 0.5, p0 = 0.2, Nj = c(-10, 90)),
    "Nj must be a vector of positive integers"
  )
  expect_error(
    rcp1armBinary(p = 0.5, p0 = 0.2, Nj = c(10, 90), PI = -0.1),
    "PI must be a single numeric value in"
  )
  expect_error(
    rcp1armBinary(p = 0.5, p0 = 0.2, Nj = c(10, 90), approach = "invalid"),
    'approach must be either'
  )
})

test_that("rcp1armBinary: print method runs without error", {
  result <- rcp1armBinary(
    p = 0.5, p0 = 0.2, Nj = c(10, 90), approach = "formula"
  )
  expect_output(print(result), "Regional Consistency Probability")
  expect_output(print(result), "Binary")
  expect_output(print(result), "Method 1")
  expect_output(print(result), "Method 2")
})
