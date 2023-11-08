# Pareto-Optimization via Normal Boundary Intersection: 3 continuous objectives
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com


#' ParetoR_3C
#'
#' Command function to optimize 3 non-adverse impact objectives via NBI
#' algorithm
#' @param Rx  Matrix with intercorrelations among predictors
#' @param Rxy1  Vector with correlation between each predictor and non-adverse impact objective 1
#' @param Rxy2  Vector with correlation between each predictor and non-adverse impact objective 2
#' @param Rxy3  Vector with correlation between each predictor and non-adverse impact objective 3
#' @param Spac Number of solutions
#' @return out Pareto-Optimal solution with objective outcome values (Criterion) and
#'   predictor weights (ParetoWeights)
#' @export
ParetoR_3C = function(Rx,
                      Rxy1,
                      Rxy2,
                      Rxy3,
                      # Spac = 20) { #SHINY: Spac = 20, has additional arguments graph = T and display_solution = T
                      Spac = 10) { #SHINY: Spac = 20, has additional arguments graph = T and display_solution = T
  assign("Rx_3C", Rx, rMOSTenv_3C)
  assign("Rxy1_3C", Rxy1, rMOSTenv_3C)
  assign("Rxy2_3C", Rxy2, rMOSTenv_3C)
  assign("Rxy3_3C", Rxy3, rMOSTenv_3C)

  # Tolerance Level for Algorithm
  TolCon = 1e-07
  TolF = 1e-15
  TolX = 1e-15

  # Do not change
  Fnum = 3
  Xnum = dim(Rx)[1]
  VLB = c(rep(0, Xnum))
  VUB = c(rep(1, Xnum))
  X0 = c(rep(1 / Xnum, Xnum))

  ###### Find Pareto-Optimal Solution ######
  out = suppressWarnings(NBI_3C(X0,
                                Spac, Fnum, VLB, VUB, TolX, TolF,
                                TolCon)) #SHINY: has additional arguments graph = T and display_solution = T
  return(out)
}


####################################### NBI Main Function ####################################

#' NBI Main Function
#'
#' Main function for obtaining pareto-optimal solution via normal boundary intersection.
#' @param X0 Initial input for predictor weight vector
#' @param Spac Number of Pareto spaces (i.e., number of Pareto points minus one)
#' @param Fnum Number of criterions
#' @param VLB Lower boundary for weight vector estimation
#' @param VUB Upper boundary for weight vector estimation
#' @param TolX Tolerance index for estimating weight vector, default is 1e-4
#' @param TolF Tolerance index for estimating criterion, default is 1e-4
#' @param TolCon Tolerance index for constraint conditions, default is 1e-7
#' @import nloptr
#' @return Pareto-Optimal solutions
#' @keywords internal

NBI_3C = function(X0,
                  Spac,
                  Fnum,
                  VLB = vector(),
                  VUB = vector(),
                  TolX = 1e-07,
                  TolF = 1e-07,
                  TolCon = 1e-07) { #SHINY: has additional parameters (graph = T, display_solution = TRUE)

  #------------------------------Initialize Options-------------------------------#
  X0 = assert_col_vec_3C(X0)
  VLB = assert_col_vec_3C(VLB)
  VUB = assert_col_vec_3C(VUB)

  # Number of variables
  nvars = length(X0)
  opts = list(
    maxeval = (nvars + 1) *
      1000,
    xtol_rel = TolX,
    ftol_rel = TolF
  )

  #Initialize PHI
  PHI = matrix(0, Fnum, Fnum)

  #------------------------------Shadow Point-------------------------------#

  # cat('\n ----Step 1: find shadow minimum---- \n')

  ShadowF = matrix(0, Fnum)
  ShadowX = matrix(0, nvars, Fnum)
  out = WeightsFun_3C(Fnum, Spac)
  Weight = out$Weights
  Near = out$Formers
  rm(out)
  Weight = Weight / Spac
  for (i in 1:Fnum) {
    temp = c(1, (Spac + 1), dim(Weight)[2])
    j = temp[i]
    assign("g_Weight_3C", Weight[, j],
           rMOSTenv_3C)
    out = suppressMessages(
      nloptr::slsqp(
        x0 = X0,
        fn = myLinCom_3C,
        lower = VLB,
        upper = VUB,
        hin = myCon_ineq_3C,
        heq = myCon_eq_3C,
        control = list(
          maxeval = (nvars +
                       1) * 1000,
          xtol_rel = TolX,
          ftol_rel = TolF
        )
      )
    )
    x = out$par
    x[1:nvars] = x[1:nvars] / sum(x[1:nvars])
    x[2:nvars] = x[2:nvars] / sum(x[2:nvars])
    ShadowX[, i] = x
    tempF = myFM_3C(x)
    ShadowF[i] = tempF[i]
  }

  #------------------------------Matrix PHI-------------------------------#

  # cat('\n ----Step 2: find PHI---- \n')

  for (i in 1:Fnum) {
    PHI[, i] = myFM_3C(ShadowX[, i]) -
      ShadowF
    PHI[i, i] = 0
  }
  assign("ShadowF_3C", ShadowF, rMOSTenv_3C)
  assign("PHI_3C", PHI, rMOSTenv_3C)

  #------------------------------Quasi-Normal Direction-------------------------------#

  # cat('\n ----Step 3: find Quasi-Normal---- \n')

  assign("g_Normal_3C", -rMOSTenv_3C$PHI_3C %*%
           matrix(1, Fnum, 1), rMOSTenv_3C)
  assign("g_Normal_3C", rMOSTenv_3C$g_Normal_3C / (norm(matrix(
    rMOSTenv_3C$g_Normal_3C
  ),
  "2") ^ 2), rMOSTenv_3C)

  #------------------------------weights-------------------------------#

  # cat('\n ----Step 4: create weights---- \n')

  out = WeightsFun_3C(Fnum, Spac)
  Weight = out$Weight
  Near = out$Formers
  Weight = Weight / Spac
  num_w = dimFun_3C(Weight)[2]

  #------------------------------NBI Subproblems-------------------------------#

  # Starting point for first NBI subproblem is the minimizer of f_1(x)

  xstart = c(ShadowX[, 1], 0) # changed based on try_NBI3D
  Pareto_Fmat = vector()
  Pareto_Xmat = vector()
  X_Near = vector()

  # solve NBI subproblems
  for (k in 1:num_w) {
    w = Weight[, k]
    if (Near[k] > 0) {
      xstart = X_Near[, Near[k]]
      # start X is the previous weight-order's X
    }
    assign("g_StartF_3C", rMOSTenv_3C$PHI_3C %*%
             w + rMOSTenv_3C$ShadowF_3C, rMOSTenv_3C)

    # SOLVE NBI SUBPROBLEM
    VLB_trials = c(VLB, -Inf)
    out = suppressMessages(
      nloptr::slsqp(
        x0 = xstart,
        fn = myT_3C,
        lower = VLB_trials,
        upper = c(VUB, Inf),
        hin = myCon_ineq_3C,
        heq = myTCon_eq_3C,
        control = opts
      )
    )
    x_trial = out$par
    rm(out)
    x_trial[1:nvars] = x_trial[1:nvars] / sum(x_trial[1:nvars])
    x_trial[2:nvars] = x_trial[2:nvars] / sum(x_trial[2:nvars])
    Pareto_Xmat = cbind(Pareto_Xmat,
                        x_trial[1:nvars]) # Pareto optima in X-space
    Pareto_Fmat = cbind(Pareto_Fmat,
                        -myFM_3C(x_trial[1:nvars])) # Pareto optima in F-space
    X_Near = cbind(X_Near, x_trial)
  }

  #------------------------------Plot Solutions-------------------------------#

  # if(graph==TRUE){plotPareto_3C(t(Pareto_Fmat))}

  #------------------------------Output Solutions-------------------------------#

  #   Output Solution

  Pareto_Fmat = t(Pareto_Fmat)
  Pareto_Xmat = t(Pareto_Xmat)
  colnames(Pareto_Fmat) = c("Ry1", "Ry2",
                            "Ry3")
  colnames(Pareto_Xmat) = c(paste0(rep("P",
                                       (nvars)), 1:(nvars)))

  # if(display_solution == TRUE){
  #
  #   solution = round(cbind(Pareto_Fmat,Pareto_Xmat),3)
  #   colnames(solution) = c("Ry1","Ry2","Ry3",paste0(rep("P",(nvars)),1:(nvars)))
  #   cat("\n Pareto-Optimal Solution \n \n")
  #   print(solution)
  #
  # }else{
  # cat("\n Done. \n \n")
  # }

  return(list(Pareto_Fmat = Pareto_Fmat,
              Pareto_Xmat = Pareto_Xmat))
}


########################### Supporting Functions ########################

# Supplementary Functions for NBI.r - Pareto-Optimization via Normal Boundary Intersection

###### combR_3C()######

#' combR_3C
#'
#' Support function to create predictor-criterion matrix
#' @param Rx Predictor inter-correlation matrix
#' @param Ry Predictor-criterion correlation (validity)
#' @return Rxy Predictor-criterion correlation matrix
#' @keywords internal
combR_3C <- function(Rx, Ry) {
  cbind(rbind(Rx, c(Ry)), c(Ry, 1))
}


###### myFM_3C() ######

#' myFM_3C
#'
#' Supporting function, defines criterion space
#' @param x Input predictor weight vector
#' @return f Criterion vector
#' @keywords internal
myFM_3C = function(x) {
  Rx <- rMOSTenv_3C$Rx_3C #SHINY: does not have this line, but likely needed for package
  Rxy1 <- rMOSTenv_3C$Rxy1_3C #SHINY: does not have this line, but likely needed for package
  Rxy2 <- rMOSTenv_3C$Rxy2_3C #SHINY: does not have this line, but likely needed for package
  Rxy3 <- rMOSTenv_3C$Rxy3_3C #SHINY: does not have this line, but likely needed for package
  b = x

  # Composite Validity R_cy1, Rcy2
  R_cy1 = t(c(t(b), 0) %*% combR_3C(Rx,
                                 Rxy1) %*% c(t(matrix(0, dimFun_3C(
                                   Rx
                                 )[1],
                                 1)), 1)) / sqrt(t(b) %*% Rx %*%
                                                   b)
  R_cy2 = t(c(t(b), 0) %*% combR_3C(Rx,
                                 Rxy2) %*% c(t(matrix(0, dimFun_3C(
                                   Rx
                                 )[1],
                                 1)), 1)) / sqrt(t(b) %*% Rx %*%
                                                   b)
  R_cy3 = t(c(t(b), 0) %*% combR_3C(Rx,
                                 Rxy3) %*% c(t(matrix(0, dimFun_3C(
                                   Rx
                                 )[1],
                                 1)), 1)) / sqrt(t(b) %*% Rx %*%
                                                   b)
  f = matrix(1, 3, 1)
  f[1, ] = -R_cy1
  f[2, ] = -R_cy2
  f[3, ] = -R_cy3
  return(f)
}


####### myCon_ineq_3C() ######

# Nonlinear inequalities at x

#' myCon_ineq_3C
#'
#' Support function, defines inequal constraint condition
#' @param x Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
#' @keywords internal

myCon_ineq_3C = function(x) {
  return(vector())
}


####### myCon_ineq_3C() ######

# Nonlinear inequalities at x

#' myCon_ineq_3C
#'
#' Support function, defines inequal constraint condition
#' @param x Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
#' @keywords internal

myCon_eq_3C = function(x) {
  Rx <- rMOSTenv_3C$Rx_3C #SHINY: does not have this line, but likely needed for package
  b = x
  ceq = (t(b) %*% Rx %*% b) - 1
  return(ceq)
}


####### myCon_ineq_3C() ######

# Nonlinear inequalities at x

#' myCon_ineq_3C
#'
#' Support function, defines inequal constraint condition
#' @param v Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
#' @keywords internal

assert_col_vec_3C = function(v) {
  if (is.null(dimFun_3C(v))) {
    v = v
  }
  else if (dimFun_3C(v)[1] < dimFun_3C(v)[2]) {
    v = t(t)
  }
  return(v)
} #SHINY: has an additional function assert_col_vec1a() which is not used


###### dimFun_3C() ######

#' dimFun_3C
#'
#' Support function, checks input predictor weight vector x
#' @param x Input predictor weight vector
#' @returns x Checked and refined input predictor weight vector
#' @keywords internal

dimFun_3C = function(x) {
  if (is.null(dim(x))) {
    return(c(0, 0))
  }
  else
    (return(dim(x)))
}


###### WeightsFun_3C() ######

#' WeightsFun_3C
#'
#' Support function, generates all possible weights for NBI subproblems
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weights All possible weights for NBI subproblem
#' @keywords internal

WeightsFun_3C = function(n, k) {

  # package environment variables
  # weight, Weights, Formers, Layer, lastone, currentone
  #
  # Generates all possible weights for NBI subproblems given:
  # n, the number of objectives
  # 1/k, the uniform spacing between two w_i (k integral)
  # This is essentially all the possible integral partitions
  # of integer k into n parts.

  assign("WeightSub_3C", matrix(0, 1, n),
         rMOSTenv_3C)
  assign("Weights_3C", vector(), rMOSTenv_3C)
  assign("Formers_3C", vector(), rMOSTenv_3C)
  assign("Layer_3C", n, rMOSTenv_3C)
  assign("lastone_3C", -1, rMOSTenv_3C)
  assign("currentone_3C", -1, rMOSTenv_3C)
  Weight_Generate_3C(1, k)
  return(list(Weights = rMOSTenv_3C$Weights_3C,
              Formers = rMOSTenv_3C$Formers_3C))
}


###### Weight_Generate_3C() ######

#' Weight_Generate_3C
#'
#' Function intended to test the weight generation scheme for NBI for > 2 objectives
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weight_Generate_3C
#' @keywords internal

Weight_Generate_3C = function(n, k) {

  # package env variables:
  # weight Weights Formers Layer lastone currentone

  # wtgener_test(n,k)
  #
  # Intended to test the weight generation scheme for NBI for > 2 objectives
  # n is the number of objectives
  # 1/k is the uniform spacing between two w_i (k integral)

  if (n == rMOSTenv_3C$Layer_3C) {
    if (rMOSTenv_3C$currentone_3C >= 0) {
      rMOSTenv_3C$Formers_3C <- c(rMOSTenv_3C$Formers_3C,
                            rMOSTenv_3C$lastone_3C)
      rMOSTenv_3C$lastone_3C <- rMOSTenv_3C$currentone_3C
      rMOSTenv_3C$currentone_3C <- -1
    }
    else {
      num = dimFun_3C(rMOSTenv_3C$Weights_3C)[2]
      rMOSTenv_3C$Formers_3C <- c(rMOSTenv_3C$Formers_3C,
                            num)
    }
    rMOSTenv_3C$WeightSub_3C[(rMOSTenv_3C$Layer_3C -
                          n + 1)] <- k
    rMOSTenv_3C$Weights_3C <- cbind(rMOSTenv_3C$Weights_3C,
                              t(rMOSTenv_3C$WeightSub_3C))
  }
  else {
    for (i in 0:k) {
      if (n == (rMOSTenv_3C$Layer_3C -
                2)) {
        num = dimFun_3C(rMOSTenv_3C$Weights_3C)[2]
        rMOSTenv_3C$currentone_3C <- num +
          1
      }
      rMOSTenv_3C$WeightSub_3C[(rMOSTenv_3C$Layer_3C -
                            n + 1)] <- i
      Weight_Generate_3C(n + 1, k -
                        i)
    }
  }
}


###### myLinCom_3C() ######

#' myLinCom_3C
#'
#' Support function
#' @param x Input predictor weight vector
#' @return f Criterion vector
#' @keywords internal

myLinCom_3C = function(x) {

  # package environment variable: g_Weight
  F = myFM_3C(x)
  f = t(rMOSTenv_3C$g_Weight_3C) %*% F
  return(f)
}


###### myT_3C() ######

#' myT_3C
#'
#' Support function, define criterion space for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return f Temporary criterion space
#' @keywords internal

myT_3C = function(x_t) {
  # f = x_t[length(x_t)]
  f = -x_t[length(x_t)]
  return(f)
}

###### myTCon_eq_3C() ######

#' myTCon_eq_3C
#'
#' Support function, define constraint condition for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return ceq Temporary constraint condition
#' @keywords internal

myTCon_eq_3C = function(x_t) {
  # package environment variables:
  # g_Normal g_StartF
  t = x_t[length(x_t)]
  x = x_t[1:(length(x_t) - 1)]
  fe = myFM_3C(x) - rMOSTenv_3C$g_StartF_3C -
    t * rMOSTenv_3C$g_Normal_3C
  ceq1 = myCon_eq_3C(x)
  ceq = c(ceq1, fe)
  return(ceq)
}

####### myTCon_ineq_3C() ######

# Nonlinear inequalities at x

#' myTCon_ineq_3C
#'
#' Support function, defines inequal constraint condition
#' @param x_t Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
#' @keywords internal

myTCon_ineq_3C = function(x_t) {
  return(vector())
}
