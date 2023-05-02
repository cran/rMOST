#' calc_out
#'
#' Helper function to calculate the expected objective outcome values based on
#' predictor weights solutions. Called by MOST().
#' @param x  Matrix of predictor weights solutions
#' @returns Expected objective outcomes
#' @keywords internal

calc_out <- function(x) {

  # package env variables:
  # optProb, Rx, Rxy1, Rxy2, Rxy3, sr, prop1, prop2, d1, d2

  # capture appropriate selection ratio and minority proportion(s)

  if (rMOSTenv$optProb == "3C") {

    Rx <- rMOSTenv$Rx_ParetoR
    Rxy1 <- rMOSTenv$Rxy1_ParetoR
    Rxy2 <- rMOSTenv$Rxy2_ParetoR
    Rxy3 <- rMOSTenv$Rxy3_ParetoR

    # predictor-objective relationships
    cmat <- t(cbind(as.matrix(Rxy1), as.matrix(Rxy2), as.matrix(Rxy3)))

  } else if (rMOSTenv$optProb == "2C_1AI") {

    Rx <- rMOSTenv$Rx_ParetoR
    Rxy1 <- rMOSTenv$Rxy1_ParetoR
    Rxy2 <- rMOSTenv$Rxy2_ParetoR
    sr <- rMOSTenv$sr_ParetoR
    prop1 <- rMOSTenv$prop1_ParetoR
    d1 <- rMOSTenv$d1_ParetoR

    # predictor-objective relationships
    cmat <- t(cbind(as.matrix(Rxy1), as.matrix(Rxy2), as.matrix(d1)))

  } else if (rMOSTenv$optProb == "1C_2AI") {

    Rx <- rMOSTenv$Rx_ParetoR
    Rxy1 <- rMOSTenv$Rxy1_ParetoR
    sr <- rMOSTenv$sr_ParetoR
    prop1 <- rMOSTenv$prop1_ParetoR
    prop2 <- rMOSTenv$prop2_ParetoR
    d1 <- rMOSTenv$d1_ParetoR
    d2 <- rMOSTenv$d2_ParetoR

    # predictor-objective relationships
    cmat <- t(cbind(as.matrix(Rxy1), as.matrix(d1), as.matrix(d2)))

  }

  # predictor correlations
  pmat <- as.matrix(Rx)

  # capture vector of predictor weights
  w <- as.vector(x)

  # calculate composite validities and/or mean subgroup differences
  cors <- as.vector(cmat %*% w / as.vector(sqrt(w %*% pmat %*% w)))

  # convert composite mean subgroup difference(s) to AI ratio(s)

  if (rMOSTenv$optProb == "2C_1AI"){

    cors[3] <- ai_ratio(d = cors[3], sr = sr, p = prop1)

  } else if (rMOSTenv$optProb == "1C_2AI"){

    cors[2] <- ai_ratio(d = cors[2], sr = sr, p = prop1)
    cors[3] <- ai_ratio(d = cors[3], sr = sr, p = prop2)

  }

  return(cors)

}

#' ai_ratio
#'
#' Helper function to convert mean subgroup differences to AI ratios (Newman et
#' al., 2007). Called by calc_out().
#' @param d  Mean subgroup difference of predictor(s)
#' @param sr Selection ratio in the full applicant pool
#' @param p Proportion of minority group in the full applicant pool
#' @keywords internal

ai_ratio <- function(d, sr, p){

  xcut <- (qnorm(1 - sr) * sqrt(1 + d^2 * (p * (1 - p)))) - (d * p)

  air <- (1.64 * xcut + sqrt(.76 * xcut ^ 2 + 4)) / (sqrt(exp(d^2 + 2 * xcut * d)) * (1.64 * (xcut + d) + sqrt(.76 * (xcut + d)^2 + 4)))

  return(air)

}
