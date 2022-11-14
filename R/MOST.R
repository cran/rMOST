# Multi-Objective Optimization via Normal Boundary Intersection for 3 Objectives
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com

#' MOST
#'
#' Optimizes 3 objectives with normal boundary intersection algorithm
#'
#' @details
#' # Inputs required by optimization problems
#' Different types of optimization problems require different input parameters:
#' * optProb = "3C": MOST(optProb, Rx, Rxy1, Rxy2, Rxy3)
#' * optProb = "2C_1AI": MOST(optProb, Rx, Rxy1, Rxy2, sr, prop1, d1)
#' * optProb = "1C_2AI": MOST(optProb, Rx, Rxy1, sr, prop1, d1, prop2, d2)
#'
#' # Notes regarding the inputs
#' * For personnel selection applications, all predictor-intercorrelations and criterion-related validity inputs should be corrected for range restriction and criterion unreliability to reflect the relations in the applicant sample.
#' * For optimization problems with 2 adverse impact objectives (i.e., optProb = "1C_2AI"), d1 and d2 should be the standardized mean difference between a minority group and the same reference group (e.g., Black-White and Hispanic-White, not Black-White and female-male)
#'
#' # Optimization
#' * Optimization may take several minutes to run.
#' * Optimization may fail in some applications due to non-convergence.
#'
#' For more details, please consult the vignette.
#'
#' @param optProb Optimization problem. "3C" = no adverse impact objectives and three non-adverse impact objectives; "2C_1AI" = one adverse impact objective and two non-adverse impact objectives; "1C_2AI" = two adverse impact objectives and one non-adverse impact objective.
#' @param Rx Predictor intercorrelation matrix
#' @param Rxy1 Needs to specify for all three types of optimization problems (optProb). Predictor criterion-related validity for non-adverse impact objective 1 (i.e., correlation between each predictor and non-adverse impact objective 1)
#' @param Rxy2 Only specify if optimization problem is "3C" or "2C_1AI". Predictor criterion-related validity for non-adverse impact objective 2 (i.e., correlation between each predictor and non-adverse impact objective 2)
#' @param Rxy3 Only specify if optimization problem is "3C". Predictor criterion-related validity for non-adverse impact objective 3 (i.e., correlation between each predictor and non-adverse impact objective 3)
#' @param sr Only specify if optimization problem is "2C_1AI" or "1C_2AI". Overall selection ratio.
#' @param prop1 Only specify if optimization problem is "2C_1AI" or "1C_2AI". Proportion of minority1 in the applicant pool; prop1 = (# of minority1 applicants)/(total # of applicants)
#' @param prop2 Only specify if optimization problem is "1C_2AI". Proportion of minority2 in the applicant pool; prop2 = (# of minority2 applicants)/(total # of applicants)
#' @param d1 Only specify if optimization problem is "2C_1AI" or "1C_2AI". Vector of standardized group-mean differences between majority and minority 1 for each predictor; d1 = avg_majority - avg_minority1
#' @param d2 Only specify if optimization problem is "1C_2AI". Vector of standardized group-mean differences between majority and minority 2 for each predictor; d2 = avg_majority - avg_minority2
#' @param Spac Determines the number of solutions.
#' @returns Pareto-Optimal solutions with objective values (e.g., C1, AI1) and the corresponding predictor weights (e.g., P1, P2)
#' @examples
#' # A sample optimization problem with 3 non-adverse impact objectives and 3 predictors
#' # For more examples, please consult the vignette.
#'
#' # Specify inputs
#' # Predictor inter-correlation matrix (Rx)
#' Rx <- matrix(c(1,  .50, .50,
#'                .50,  1, .50,
#'                .50, .50,  1), 3, 3)
#'
#' # Predictor-objective relation vectors (Rxy1, Rxy2, Rxy3)
#' # Criterion-related validities
#' ## Criterion 1
#' Rxy1 <- c(-.30, 0, .30)
#' ## Criterion 2
#' Rxy2 <- c(0, .30, -.30)
#' ## Criterion 3
#' Rxy3 <- c(.30, -.30, 0)
#'
#' # Get Pareto-optimal solutions
#' \donttest{
#' out <- MOST(optProb = "3C", Rx = Rx, Rxy1 = Rxy1, Rxy2 = Rxy2, Rxy3 = Rxy3, Spac = 10)
#' out
#' }
#' @export

MOST <- function(optProb, Rx, Rxy1, Rxy2, Rxy3, sr, prop1, prop2, d1, d2, Spac = 10){

  message('\n Estimating Multi-Objective Optimal Solution ... \n')

  # optProb <<- optProb
  assign("optProb", optProb, rMOSTenv)

  ## Obtain MOO Predictor Weights ##

  if (rMOSTenv$optProb == "3C"){

    # Rx_ParetoR <<- Rx
    # Rxy1_ParetoR <<- Rxy1
    # Rxy2_ParetoR <<- Rxy2
    # Rxy3_ParetoR <<- Rxy3
    assign("Rx_ParetoR", Rx, rMOSTenv)
    assign("Rxy1_ParetoR", Rxy1, rMOSTenv)
    assign("Rxy2_ParetoR", Rxy2, rMOSTenv)
    assign("Rxy3_ParetoR", Rxy3, rMOSTenv)

    out1 <- ParetoR_3C(Rx = Rx, Rxy1 = Rxy1, Rxy2 = Rxy2, Rxy3 = Rxy3, Spac = Spac)
    out2.1 <- ParetoR_2C(Rx = Rx, Rxy1 = Rxy1, Rxy2 = Rxy2, Spac = Spac, graph = FALSE)
    out2.2 <- ParetoR_2C(Rx = Rx, Rxy1 = Rxy1, Rxy2 = Rxy3, Spac = Spac, graph = FALSE)
    out2.3 <- ParetoR_2C(Rx = Rx, Rxy1 = Rxy2, Rxy2 = Rxy3, Spac = Spac, graph = FALSE)

  } else if (rMOSTenv$optProb == "2C_1AI"){

    # Rx_ParetoR <<- Rx
    # Rxy1_ParetoR <<- Rxy1
    # Rxy2_ParetoR <<- Rxy2
    # sr_ParetoR <<- sr
    # prop1_ParetoR <<- prop1
    # d1_ParetoR <<- d1
    assign("Rx_ParetoR", Rx, rMOSTenv)
    assign("Rxy1_ParetoR", Rxy1, rMOSTenv)
    assign("Rxy2_ParetoR", Rxy2, rMOSTenv)
    assign("sr_ParetoR", sr, rMOSTenv)
    assign("prop1_ParetoR", prop1, rMOSTenv)
    assign("d1_ParetoR", d1, rMOSTenv)

    out1 <- ParetoR_2C_1AIR(Rx = Rx, Rxy1 = Rxy1, Rxy2 = Rxy2, sr = sr, prop1 = prop1, d1 = d1, Spac = Spac)
    out2.1 <- ParetoR_2C(Rx = Rx, Rxy1 = Rxy1, Rxy2 = Rxy2, Spac = Spac, graph = FALSE)
    out2.2 <- ParetoR(Rx = Rx, Rxy1 = Rxy1, sr = sr, prop1 = prop1, d1 = d1, Spac = Spac, graph = FALSE)
    out2.3 <- ParetoR(Rx = Rx, Rxy1 = Rxy2, sr = sr, prop1 = prop1, d1 = d1, Spac = Spac, graph = FALSE)

  } else if (rMOSTenv$optProb == "1C_2AI"){

    # Rx_ParetoR <<- Rx
    # Rxy1_ParetoR <<- Rxy1
    # sr_ParetoR <<- sr
    # prop1_ParetoR <<- prop1
    # prop2_ParetoR <<- prop2
    # d1_ParetoR <<- d1
    # d2_ParetoR <<- d2
    assign("Rx_ParetoR", Rx, rMOSTenv)
    assign("Rxy1_ParetoR", Rxy1, rMOSTenv)
    assign("sr_ParetoR", sr, rMOSTenv)
    assign("prop1_ParetoR", prop1, rMOSTenv)
    assign("prop2_ParetoR", prop2, rMOSTenv)
    assign("d1_ParetoR", d1, rMOSTenv)
    assign("d2_ParetoR", d2, rMOSTenv)

    out1 <- ParetoR_1C_2AIR(Rx = Rx, Rxy1 = Rxy1, sr = sr, prop1 = prop1, prop2 = prop2, d1 = d1, d2 = d2, Spac = Spac)
    out2.1 <- ParetoR(Rx = Rx, Rxy1 = Rxy1, sr = sr, prop1 = prop1, d1 = d1, Spac = Spac)
    out2.2 <- ParetoR(Rx = Rx, Rxy1 = Rxy1, sr = sr, prop1 = prop2, d1 = d2, Spac = Spac)
    out2.3 <- NULL

  }

  ## Obtain Hiring Outcomes ##

  DF1a <- as.data.frame(round(out1$Pareto_Xmat,3))
  DF1b <-  as.data.frame(matrix(apply(DF1a, 1, calc_out), ncol=3, byrow=TRUE))
  colnames(DF1b) <- c("Ry1","Ry2","Ry3")
  rownames(DF1b) <- c(1:nrow(DF1b))
  tab1 <- as.data.frame(cbind(DF1b,DF1a))

  DF2.1a <- as.data.frame(round(out2.1$Pareto_Xmat,3))
  DF2.1b <-  as.data.frame(matrix(apply(DF2.1a, 1, calc_out), ncol=3, byrow=TRUE))
  colnames(DF2.1b) <- c("Ry1","Ry2","Ry3")
  rownames(DF2.1b) <- c(1:nrow(DF2.1b))
  tab2.1 <- as.data.frame(cbind(DF2.1b,DF2.1a))

  DF2.2a <- as.data.frame(round(out2.2$Pareto_Xmat,3))
  DF2.2b <-  as.data.frame(matrix(apply(DF2.2a, 1, calc_out), ncol=3, byrow=TRUE))
  colnames(DF2.2b) <- c("Ry1","Ry2","Ry3")
  rownames(DF2.2b) <- c(1:nrow(DF2.2b))
  tab2.2 <- as.data.frame(cbind(DF2.2b,DF2.2a))

  if (rMOSTenv$optProb == "3C" | rMOSTenv$optProb == "2C_1AI") {

    DF2.3a <- as.data.frame(round(out2.3$Pareto_Xmat,3))
    DF2.3b <-  as.data.frame(matrix(apply(DF2.3a, 1, calc_out), ncol=3, byrow=TRUE))
    colnames(DF2.3b) <- c("Ry1","Ry2","Ry3")
    rownames(DF2.3b) <- c(1:nrow(DF2.3b))
    tab2.3 <- as.data.frame(cbind(DF2.3b,DF2.3a))

  }

  ## Organize Outputs ##

  # Prepare output table

  if (rMOSTenv$optProb == "3C") {critNames <- c("C1","C2","C3")}
  else if (rMOSTenv$optProb == "2C_1AI") {critNames <- c("C1","C2","AI1")}
  else if (rMOSTenv$optProb == "1C_2AI") {critNames <- c("C1","AI1","AI2")}

  predNames <- c(paste("P", 1:dim(rMOSTenv$Rx_ParetoR)[1], sep = ""))

  names <- c("Solution", "Optimized", critNames, paste0("P",1:length(predNames)))

  # Fill in output table

  tabl1 <- cbind(c(1:nrow(tab1)), rep(paste(critNames[1:3], collapse=", "), nrow(tab1)), tab1)
  colnames(tabl1) <- names
  tabl2 <- cbind(c((nrow(tab1)+1):(nrow(tab1)+nrow(tab2.1))), rep(paste(critNames[1:2], collapse=", "), nrow(tab2.1)), tab2.1)
  colnames(tabl2) <- names
  tabl3 <- cbind(c((nrow(tab1)+nrow(tab2.1)+1):(nrow(tab1)+nrow(tab2.1)+nrow(tab2.2))), rep(paste(critNames[c(1,3)], collapse=", "), nrow(tab2.2)), tab2.2)
  colnames(tabl3) <- names

  if (rMOSTenv$optProb == "3C" | rMOSTenv$optProb == "2C_1AI") {

    tabl4 <- cbind(c((nrow(tab1)+nrow(tab2.1)+nrow(tab2.2)+1):(nrow(tab1)+nrow(tab2.1)+nrow(tab2.2)+nrow(tab2.3))), rep(paste(critNames[2:3], collapse=", "), nrow(tab2.3)), tab2.3)
    colnames(tabl4) <- names

    tabl <- rbind(tabl1,tabl2,tabl3,tabl4)

  } else {

    tabl <- rbind(tabl1,tabl2,tabl3)

  }

  tabl[,-2] <- round(tabl[,-2],3)

  message("\n Done. \n \n")

  return(tabl)
}
