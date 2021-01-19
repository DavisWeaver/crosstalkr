#' Identify proteins with a statistically significant relationship to user-provided seeds.
#'
#' @inheritParams bootstrap_null
#'
#' @inheritParams sparseRWR
#'
#' @inheritParams setup_init
#'
#'
#' @export

compute_crosstalk <- function(seed_proteins, ppi = "stringdb", n = 10000,
                              gamma=0.6, eps = 1e-10, tmax = 1000,
                              norm = TRUE, set_seed,
                              cache, seed_name = NULL,
                              ncores = 1)  {

}
