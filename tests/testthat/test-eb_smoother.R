test_that("eb_smoother with nbinom2() works", {
  set.seed(1234)

  test_eb_smoother_nbinom2 <- function(shape, rate,
                                       tolerance = 1e-2) {
    size <- length(population)
    observed <- rpois(size, rgamma(size, shape, rate) * population)

    eb_smoother_nbinom2 <- data.frame(observed = observed,
                                      population = population) |>
      eb_smoother(observed = "observed",
                  population = "population",
                  family = nbinom2())

    model <- eb_smoother_nbinom2$model
    family <- eb_smoother_nbinom2$family

    expect_equal(stats::sigma(model), shape,
                 tolerance = tolerance)
    expect_equal(family$linkinv(stats::coef(summary(model))$cond["(Intercept)", "Estimate"]), shape / rate,
                 tolerance = tolerance)
  }

  population <- sample(seq_len(1e2), 1e5,
                       replace = TRUE)

  test_eb_smoother_nbinom2(2, 3)
  test_eb_smoother_nbinom2(3, 2)
  test_eb_smoother_nbinom2(3, 5)
  test_eb_smoother_nbinom2(5, 3)
})

test_that("eb_smoother with betabinomial() works", {
  set.seed(1234)

  test_eb_smoother_betabinomial <- function(shape1, shape2,
                                            tolerance = 1e-2) {
    size <- length(population)
    observed <- rbinom(size, population, rbeta(size, shape1, shape2))

    eb_smoother_nbinom2 <- data.frame(observed = observed,
                                      population = population) |>
      eb_smoother(observed = "observed",
                  population = "population",
                  family = betabinomial())

    model <- eb_smoother_nbinom2$model
    family <- eb_smoother_nbinom2$family

    expect_equal(stats::sigma(model), shape1 + shape2,
                 tolerance = tolerance)
    expect_equal(family$linkinv(stats::coef(summary(model))$cond["(Intercept)", "Estimate"]), shape1 / (shape1 + shape2),
                 tolerance = tolerance)
  }

  population <- sample(seq_len(1e2), 1e5,
                       replace = TRUE)

  test_eb_smoother_betabinomial(2, 3)
  test_eb_smoother_betabinomial(3, 2)
  test_eb_smoother_betabinomial(3, 5)
  test_eb_smoother_betabinomial(5, 3)
})
