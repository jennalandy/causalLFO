#' Summary of a causalLFO_result object
#'
#' @param object causalLFO_result object
#'
#' @method summary causalLFO_result
#' @export
summary.causalLFO_result <- function(object) {
  print(data.frame(
    ATE = object$ATE
  ))
}

#' Summary of a causalLFO_bootstrap_result object
#'
#' @param object causalLFO_bootstrap_result object
#'
#' @method summary causalLFO_bootstrap_result
#' @export
summary.causalLFO_bootstrap_result <- function(object) {
  data.frame(
    mean = object$mean,
    lower = object$lower95,
    upper = object$upper95
  )
}

#' Plot a causalLFO_result object
#'
#' @param object causalLFO_result object
#'
#' @method plot causalLFO_result
#' @export
#' @import ggplot2
plot.causalLFO_result <- function(x) {
  K = length(x$ATE)
  data.frame(
    k = 1:K,
    ATE = x$ATE
  ) %>%
    ggplot2::ggplot(ggplot2::aes(x = k, y = ATE)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 1), size = 3) +
    ggplot2::labs(
      x = 'Latent Dimension',
      y = "ATE Estimate"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20)
    ) +
    ggplot2::geom_hline(yintercept = 0)
}

#' Plot multiple causalLFO_result objects
#'
#' @param res_list named list of causalLFO_result objects. Names will become colors
#' @param group_name string, name of grouping variable / legend title for colors
#'
#' @return ggplot2 object
#'
#' @export
#' @import ggplot2
plot_causalLFO_results <- function(res_list, group_name = "") {
  lapply(names(res_list), function(name) {
    x = res_list[[name]]
    K = length(x$ATE)
    data.frame(
      k = 1:K,
      ATE = x$ATE,
      name = name
    )
  }) %>%
    do.call(rbind, .) %>%
    ggplot2::ggplot(aes(x = k, y = ATE, color = name)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 1), size = 3) +
    ggplot2::labs(
      x = 'Latent Dimension',
      y = "ATE Estimate",
      color = group_name
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20)
    ) +
    ggplot2::geom_hline(yintercept = 0)
}

#' Plot a causalLFO_bootstrap_result object
#'
#' @param object causalLFO_result object
#'
#' @method plot causalLFO_bootstrap_result
#' @export
#' @import ggplot2
plot.causalLFO_bootstrap_result <- function(x) {
  K = length(x$mean)
  data.frame(
    k = 1:K,
    mean = x$mean,
    lower = x$lower95,
    upper = x$upper95
  ) %>%
    ggplot2::ggplot(aes(x = k, y = mean)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 1), size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), position = ggplot2::position_dodge(width = 1), width = 0.5, size = 1.5) +
    ggplot2::labs(
      x = 'Latent Dimension',
      y = "ATE Bootstrap Mean and 95% CI"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20)
    ) +
    ggplot2::geom_hline(yintercept = 0)
}


#' Plot multiple causalLFO_bootstrap_result objects
#'
#' @param res_list named list of causalLFO_bootstrap_result objects. Names will become colors
#' @param group_name string, name of grouping variable / legend title for colors
#'
#' @return ggplot2 object
#'
#' @export
#' @import ggplot2
plot_causalLFO_bootstrap_results <- function(res_list, group_name = "") {
  lapply(names(res_list), function(name) {
    x = res_list[[name]]
    K = length(x$mean)
    data.frame(
      k = 1:K,
      mean = x$mean,
      lower = x$lower95,
      upper = x$upper95,
      name = name
    )
  }) %>%
    do.call(rbind, .) %>%
    ggplot2::ggplot(aes(x = k, y = mean, color = name)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 1), size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), position = ggplot2::position_dodge(width = 1), width = 0.5, size = 1.5) +
    ggplot2::labs(
      x = 'Latent Dimension',
      y = "ATE Bootstrap Mean and 95% CI",
      color = group_name
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20)
    ) +
    ggplot2::geom_hline(yintercept = 0)
}
