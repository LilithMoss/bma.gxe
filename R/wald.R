#' The P-value from a Wald Test
#'
#' @param est input A number.
#' @param se input A number.
#' @return The P-value.
#' @export

wald <- function(est,se){
  chi <- (est/se)^2
  pval <- pchisq(chi,1,lower.tail=F)
  return(pval)
}
