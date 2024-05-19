#' Empirical Bayes Smoother
#'
#' @param data A data frame.
#' @param observed A variable representing the observed counts specified in the
#' same way as [dplyr::pull()].
#' @param population A variable representing the population counts specified in
#' the same way as [dplyr::pull()].
#' @param family A family object created by [nbinom2()] or [betabinomial()].
#' @param formula A one-sided formula specifying the model.
#' @param ... Additional arguments passed to [glmmTMB::glmmTMB()].
#'
#' @return An object of class `eb_smoother`.
#'
#' @export
eb_smoother <- function(data, observed, population,
                        family = nbinom2(),
                        formula = ~ 1,
                        ...) {
  names_data <- names(data)
  observed <- rlang::sym(tidyselect::vars_pull(names_data, {{ observed }}))
  population <- rlang::sym(tidyselect::vars_pull(names_data, {{ population }}))

  if (!is.null(rlang::f_lhs(formula))) {
    cli::cli_abort("{.arg formula} must be a one-sided formula.")
  }

  if (!family$family %in% c("nbinom2", "betabinomial")) {
    cli::cli_abort("{.arg family} must be either {.code nbinom2()} or {.code betabinomial()}.")
  } else if (family$family == "nbinom2") {
    if (family$link != "log") {
      cli::cli_abort('Only {.code link = "log"} is supported for {.code nbinom2()}.')
    }

    rlang::f_lhs(formula) <- observed
    rlang::f_rhs(formula) <- rlang::expr(!!rlang::f_rhs(formula) + offset(log(!!population)))

    model <- glmmTMB::glmmTMB(formula = formula,
                              data = data,
                              family = family,
                              ...)
  } else if (family$family == "betabinomial") {
    rlang::f_lhs(formula) <- rlang::expr(cbind(!!observed, !!population - !!observed))

    model <- glmmTMB::glmmTMB(formula = formula,
                              data = data,
                              family = family,
                              ...)
  }
  structure(list(data = data,
                 observed = observed,
                 population = population,
                 family = family,
                 model = model),
            class = "eb_smoother")
}

#' @export
print.eb_smoother <- function(x, ...) {
  cli::cli_h2("Empirical Bayes Smoother")

  cli::cli_bullets(c("*" = "Observed: {.var {x$observed}}",
                     "*" = "Population: {.var {x$population}}"))

  cli::cli_h3("Family")
  print(x$family)

  cli::cli_h3("Model")
  print(x$model)
}

#' @export
predict.eb_smoother <- function(object,
                                newdata = NULL,
                                ...) {
  expected <- predict(object$model,
                      newdata = newdata,
                      type = "response",
                      ...)
  newdata <- newdata %||% object$data
  observed <- rlang::eval_tidy(object$observed, newdata)
  population <- rlang::eval_tidy(object$population, newdata)

  if (object$family$family == "nbinom2") {
    prior <- tibble::tibble(shape = stats::sigma(object$model),
                            rate = stats::sigma(object$model) * population / expected)

    posterior <- tibble::tibble(shape = observed + prior$shape,
                                rate = population + prior$rate)
    expected <- posterior$shape / posterior$rate
    variance <- posterior$shape / posterior$rate^2
  } else if (object$family$family == "betabinomial") {
    sigma_model <- stats::sigma(object$model)

    prior <- tibble::tibble(shape1 = expected * sigma_model,
                            shape2 = sigma_model - expected * sigma_model)
    posterior <- tibble::tibble(shape1 = observed + prior$shape1,
                                shape2 = population - observed + prior$shape2)
    expected <- posterior$shape1 / (posterior$shape1 + posterior$shape2)
    variance <- posterior$shape1 * posterior$shape2 / (posterior$shape1 + posterior$shape2)^2 / (posterior$shape1 + posterior$shape2 + 1)
  }
  tibble::tibble(prior = prior,
                 posterior = posterior,
                 expected = expected,
                 variance = variance)
}
