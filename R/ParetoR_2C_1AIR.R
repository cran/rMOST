# Command Function ParetoR_2C_1AIR: Optimize 2 criteria (continuous) and AI ratio of 1 minority group
## Developer: Q. Chelsea Song
## Contact: qianqisong@gmail.com

#' ParetoR_2C_1AIR
#'
#' Command function to optimize 2 non-adverse impact objectives and 1 adverse
#' impact objective via NBI algorithm
#' @param Rx  Matrix with intercorrelations among predictors
#' @param Rxy1  Vector with correlation between each predictor and non-adverse
#'   impact objective 1
#' @param Rxy2  Vector with correlation between each predictor and non-adverse
#'   impact objective 2
#' @param prop1 Proportion of minority applicants in full applicant pool
#' @param sr Selection ratio in full applicant pool
#' @param d1 Subgroup difference; standardized mean differences between minority
#'   and majority subgroups on each predictor in full applicant pool
#' @param Spac  Number of Pareto points
#' @return out Pareto-Optimal solution with objective outcome values (Criterion)
#'   and predictor weights (ParetoWeights)
#' @keywords internal

ParetoR_2C_1AIR = function(Rx, Rxy1, Rxy2, sr, prop1, d1, Spac = 10){

  # Rx_2C_1AIR <<- Rx
  # Rxy1_2C_1AIR <<- Rxy1
  # Rxy2_2C_1AIR <<- Rxy2
  # sr_2C_1AIR <<- sr
  # prop1_2C_1AIR <<- prop1
  # d1_2C_1AIR <<- d1
  assign("Rx_2C_1AIR", Rx, rMOSTenv_2C_1AIR)
  assign("Rxy1_2C_1AIR", Rxy1, rMOSTenv_2C_1AIR)
  assign("Rxy2_2C_1AIR", Rxy2, rMOSTenv_2C_1AIR)
  assign("sr_2C_1AIR", sr, rMOSTenv_2C_1AIR)
  assign("prop1_2C_1AIR", prop1, rMOSTenv_2C_1AIR)
  assign("d1_2C_1AIR", d1, rMOSTenv_2C_1AIR)

  # Tolerance Level for Algorithm
  TolCon 	= 1e-7 # tolerance of constraint
  TolF 	= 1e-15 # tolerance of objective function
  TolX 	= 1e-15 # tolerance of predictor variable

  # Do not change
  Fnum 	= 3
  Xnum = dim(Rx)[1]
  VLB = c(-(Xnum+1),rep(0,Xnum))
  VUB = c((Xnum+1),rep(1,Xnum))
  X0 = c(1,rep(1/Xnum,Xnum))

  ###### Find Pareto-Optimal Solution ######

  out = suppressWarnings(NBI_2C_1AIR(X0, Spac, Fnum, VLB, VUB, TolX, TolF, TolCon))
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

NBI_2C_1AIR = function(X0,Spac,Fnum,VLB=vector(),VUB=vector(),TolX=1e-4,TolF=1e-4,TolCon=1e-7){

  # cat('\n Estimating Pareto-Optimal Solution ... \n')

  #------------------------------Initialize Options-------------------------------#

  X0 = assert_col_vec_2C_1AIR(X0)
  VLB = assert_col_vec_2C_1AIR(VLB)
  VUB = assert_col_vec_2C_1AIR(VUB)

  # Number of variables
  nvars = length(X0)

  # Set options
  # algorithm: sequential (least-squares) quadratic programming algorithm
  # (SQP is algorithm for nonlinearly constrained, gradient-based optimization,
  # supporting both equality and inequality constraints.)
  # maxeval: Max Iterations
  # xtol_abs: Tolerance for X
  # ftol_abs: Tolerance for F
  # tol_constraints_ineq/eq: Tolerance for inequal/equal constraints
  # (for reference) MATLAB constraints:
  # options = optimoptions('fmincon','Algorithm','sqp','MaxIter',(nvars+1)*1000,'TolFun',TolF,'TolX',TolX,'TolCon',TolCon,'Display','off')

  suppressMessages(nloptr::nl.opts(optlist = list(
    maxeval = (nvars+1)*1000
    ,xtol_rel = TolX
    ,ftol_rel = TolF
  )))

  #Initialize PHI

  PHI = matrix(0,Fnum,Fnum)

  #------------------------------Shadow Point-------------------------------#

  # cat('\n ----Step 1: find shadow minimum---- \n')

  ShadowF = matrix(0,Fnum)
  ShadowX = matrix(0,nvars,Fnum)
  xstart  = X0
  out = WeightsFun_2C_1AIR(Fnum,Spac)
  Weight = out$Weights
  Near = out$Formers
  rm(out)
  Weight = Weight/Spac

  for(i in 1:Fnum){

    temp = c(1,(Spac+1),dim(Weight)[2])
    j = temp[i]
    # g_Weight_2C_1AIR <<- Weight[,j]
    assign("g_Weight_2C_1AIR", Weight[,j], rMOSTenv_2C_1AIR)

    out = suppressMessages(nloptr::slsqp(x0 = X0, fn = myLinCom_2C_1AIR
                                         ,lower = VLB, upper = VUB
                                         ,hin = myCon_ineq_2C_1AIR
                                         ,heq = myCon_eq_2C_1AIR,
                                         control = list("maxeval" = (nvars+1)*1000,
                                                        "xtol_rel" = TolX,
                                                        "ftol_rel" = TolF)
    ))

    x = out$par

    x[1:nvars] = x[1:nvars]/sum(x[1:nvars])
    x[2:nvars] = x[2:nvars]/sum(x[2:nvars])

    ShadowX[,i] = x

    tempF = myFM_2C_1AIR(x)
    ShadowF[i] = tempF[i]

  }


  #------------------------------Matrix PHI-------------------------------#

  # cat('\n ----Step 2: find PHI---- \n')

  for(i in 1:Fnum){

    PHI[,i] = myFM_2C_1AIR(ShadowX[,i]) + ShadowF
    PHI[i,i] = 0

  }

  # print(round(PHI,3))

  assign("ShadowF_2C_1AIR", ShadowF, rMOSTenv_2C_1AIR)
  assign("PHI_2C_1AIR", PHI, rMOSTenv_2C_1AIR)

  #Check to make sure that QPP is n-1 dimensional
  # if(rcond(PHI_2C_1AIR) < 1e-8){stop(' Phi matrix singular, aborting.')}
  if(rcond(rMOSTenv_2C_1AIR$PHI_2C_1AIR) < 1e-8){stop(' Phi matrix singular, aborting.')}

  #------------------------------Quasi-Normal Direction-------------------------------#

  # cat('\n ----Step 3: find Quasi-Normal---- \n')

  # g_Normal_2C_1AIR <<- -PHI_2C_1AIR%*%matrix(1,Fnum,1)
  # assign("g_Normal_2C_1AIR", -PHI_2C_1AIR%*%matrix(1,Fnum,1), rMOSTenv_2C_1AIR)
  assign("g_Normal_2C_1AIR", -rMOSTenv_2C_1AIR$PHI_2C_1AIR%*%matrix(1,Fnum,1), rMOSTenv_2C_1AIR)

  #------------------------------weights-------------------------------#

  # cat('\n ----Step 4: create weights---- \n')

  out = WeightsFun_2C_1AIR(Fnum,Spac)
  Weight = out$Weight
  Near = out$Formers
  Weight = Weight/Spac
  num_w = dimFun_2C_1AIR(Weight)[2]

  # cat('\n Weights in row: \n')
  # print(round(Weight,3))

  #------------------------------NBI Subproblems-------------------------------#

  # cat('\n ----Step 5: solve NBI sub-problems---- \n')

  # Starting point for first NBI subproblem is the minimizer of f_1(x)
  xstart = c(ShadowX[,1],0)

  Pareto_Fmat = vector()       # Pareto Optima in F-space
  Pareto_Xmat = vector()       # Pareto Optima in X-space
  X_Near      = vector()

  # solve NBI subproblems
  for(k in 1:num_w){

    w  = Weight[,k]

    # Solve problem only if it is not minimizing one of the individual objectives
    indiv_fn_index = which(w == 1)
    # the boundary solution which has been solved

    if(length(indiv_fn_index) != 0){

      # w has a 1 in indiv_fn_index th component, zero in rest
      # Just read in solution from shadow data
      # Pareto_Fmat = cbind(Pareto_Fmat, (-PHI_2C_1AIR[,indiv_fn_index] + ShadowF_2C_1AIR))
      Pareto_Fmat = cbind(Pareto_Fmat, (-rMOSTenv_2C_1AIR$PHI_2C_1AIR[,indiv_fn_index] + rMOSTenv_2C_1AIR$ShadowF_2C_1AIR))
      Pareto_Xmat = cbind(Pareto_Xmat, ShadowX[,indiv_fn_index])
      X_Near = cbind(X_Near, c(ShadowX[,indiv_fn_index],0))
      # print(Pareto_Fmat)

    }else{

      w = rev(w)

      if(Near[k] > 0){

        xstart = X_Near[,Near[k]]
        #start X is the previous weight-order's X

      }

      #start point in F-space
      # g_StartF_2C_1AIR <<- PHI_2C_1AIR%*%w + ShadowF_2C_1AIR
      assign("g_StartF_2C_1AIR", rMOSTenv_2C_1AIR$PHI_2C_1AIR%*%w + rMOSTenv_2C_1AIR$ShadowF_2C_1AIR, rMOSTenv_2C_1AIR)

      # SOLVE NBI SUBPROBLEM

      out = suppressMessages(nloptr::slsqp(x0 = xstart, fn = myT_2C_1AIR
                                           ,lower = c(VLB,-Inf)
                                           ,upper = c(VUB,Inf)
                                           ,hin = myCon_ineq_2C_1AIR
                                           ,heq = myTCon_eq_2C_1AIR))

      x_trial = out$par
      f = out$value
      rm(out)

      # success
      # if(fiasco >= 0){

      Pareto_Fmat = cbind(Pareto_Fmat, -myFM_2C_1AIR(x_trial[1:nvars]))  # Pareto optima in F-space
      som = 0

      for(k in 2:nvars){som = som + x_trial[k]}

      for(k in 2:nvars){x_trial[k] = x_trial[k]/som}

      Pareto_Xmat = cbind(Pareto_Xmat, x_trial[1:nvars])        # Pareto optima in X-space
      X_Near = cbind(X_Near,x_trial)

    }

  }

  #------------------------------Plot Solutions-------------------------------#

  #   cat('\n ----Step 6: plot---- \n')

  # if(graph==TRUE){plotPareto_2C_1AIR(Pareto_Fmat, Pareto_Xmat)}

  #------------------------------Output Solutions-------------------------------#

  #   Output Solution

  Pareto_Fmat = t(Pareto_Fmat)
  Pareto_Xmat = t(Pareto_Xmat[2:nrow(Pareto_Xmat),])
  colnames(Pareto_Fmat) = c("AIR","Ry1","Ry2")
  colnames(Pareto_Xmat) = c(paste0(rep("P",(nvars-1)),1:(nvars-1)))

  # if(display_solution == TRUE){
  #
  #   solution = round(cbind(Pareto_Fmat,Pareto_Xmat),3)
  #   colnames(solution) = c("AIR","Ry1","Ry2", paste0(rep("P",(nvars-1)),1:(nvars-1)))
  #   cat("\n Pareto-Optimal Solution \n \n")
  #   print(solution)
  #
  # }else{
  # cat("\n Done. \n \n")
  # }


  return(list(Pareto_Fmat = round(Pareto_Fmat, 3),
              Pareto_Xmat = round(Pareto_Xmat, 3)))

}

########################### Supporting Functions (A) ########################

###### combR_2C_1AIR()######

#' combR_2C_1AIR
#'
#' Support function to create predictor-criterion matrix
#' @param Rx Predictor inter-correlation matrix
#' @param Ry Predictor-criterion correlation (validity)
#' @return Rxy Predictor-criterion correlation matrix
#' @keywords internal

combR_2C_1AIR <- function(Rx, Ry){
  cbind(rbind(Rx,c(Ry)),c(Ry,1))
}

###### myFM_2C_1AIR() ######

#' myFM_2C_1AIR
#'
#' Supporting function, defines criterion space
#' @param x Input predictor weight vector
#' @importFrom stats pnorm qnorm runif
#' @return f Criterion vector
#' @keywords internal

myFM_2C_1AIR = function(x){

  d <- rMOSTenv_2C_1AIR$d1_2C_1AIR
  Rx <- rMOSTenv_2C_1AIR$Rx_2C_1AIR
  Rxy1 <- rMOSTenv_2C_1AIR$Rxy1_2C_1AIR
  Rxy2 <- rMOSTenv_2C_1AIR$Rxy2_2C_1AIR

  b = x[-1]

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%Rx%*%b)

  # mean of minority weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar = 0
  # mean of majority weighted predictor composite distribution (DeCorte, 1999)
  p_a_bar = d%*%x[-1]/sigma_p
  # minority group selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR_min = 1 - pnorm(x[1], p_i_bar)
  # majority group selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR_maj = 1 - pnorm(x[1], p_a_bar)

  # AIratio a_g (DeCorte et al., 2007)
  a_g = SR_min/SR_maj

  # Composite Validity R_xy

  R_cy1 = t(c(t(b),0)%*%combR_2C_1AIR(Rx,Rxy1)%*%c(t(matrix(0,dimFun_2C_1AIR(Rx)[1],1)),1))/sqrt(t(b)%*%Rx%*%b) # DeCorte et al., 2007
  R_cy2 = t(c(t(b),0)%*%combR_2C_1AIR(Rx,Rxy2)%*%c(t(matrix(0,dimFun_2C_1AIR(Rx)[1],1)),1))/sqrt(t(b)%*%Rx%*%b) # DeCorte et al., 2007

  f = matrix(1,3,1)
  f[1,] = -a_g
  f[2,] = -R_cy1
  f[3,] = -R_cy2

  return(f)

}

####### myCon_ineq_2C_1AIR() ######

# Nonlinear inequalities at x

#' myCon_ineq_2C_1AIR
#'
#' Support function, defines inequal constraint condition
#' @param x Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
#' @keywords internal

myCon_ineq_2C_1AIR = function(x){return(vector())}

####### myCon_eq_2C_1AIR() ######

# Nonlinear equalities at x

#' myCon_eq_2C_1AIR
#'
#' Support function, defines equal constraint condition
#' @param x Input predictor weight vector
#' @importFrom stats pnorm
#' @return Equal constraint condition for use in NBI()
#' @keywords internal

myCon_eq_2C_1AIR = function(x){

  # Obtain within-package 'global' variable from package env
  # prop <- prop1_2C_1AIR
  # sr <- sr_2C_1AIR
  # d <- d1_2C_1AIR
  # Rx <- Rx_2C_1AIR
  prop <- rMOSTenv_2C_1AIR$prop1_2C_1AIR
  sr <- rMOSTenv_2C_1AIR$sr_2C_1AIR
  d <- rMOSTenv_2C_1AIR$d1_2C_1AIR
  Rx <- rMOSTenv_2C_1AIR$Rx_2C_1AIR

  b = x[-1]

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%Rx%*%b)

  # mean of minority weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar = 0
  # mean of majority weighted predictor composite distribution (DeCorte, 1999)
  p_a_bar = d%*%x[-1]/sigma_p
  # p_a_bar = (x[2]*1.00+x[3]*0.23+x[4]*0.09+x[5]*0.33)/sigma_p
  # minority group selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR_min = 1 - pnorm(x[1], p_i_bar)
  # majority group selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR_maj = 1 - pnorm(x[1], p_a_bar)

  # Nonlinear equalities at x

  ceq = matrix(1,2,1)
  ceq[1,] = SR_min*prop + SR_maj*(1-prop) - sr # DeCorte et al. (2007)
  ceq[2,] = (t(b)%*%Rx%*%b) - 1

  return(ceq)

}

########################### Supporting Functions (B) ########################

# Supplementary Functions for NBI.r - Pareto-Optimization via Normal Boundary Intersection

# Function List
## assert_col_vec_2C_1AIR
## dimFun_2C_1AIR
## WeightsFun_2C_1AIR
## Weight_Generate_2C_1AIR
## myLinCom_2C_1AIR
## myT_2C_1AIR
## myTCon_eq_2C_1AIR
## plotPareto_2C_1AIR

###### assert_col_vec_2C_1AIR() ######

#' assert_col_vec_2C_1AIR
#'
#' Support function, refines intermediate variable for use in NBI()
#' @param v Intermediate variable v
#' @return Refined variable v
#' @keywords internal

assert_col_vec_2C_1AIR = function(v){
  if(is.null(dimFun_2C_1AIR(v))){
    v=v
  }else if(dimFun_2C_1AIR(v)[1] < dimFun_2C_1AIR(v)[2]){v = t(t)}
  return(v)}

###### dimFun_2C_1AIR() ######

#' dimFun_2C_1AIR
#'
#' Support function, checks input predictor weight vector x
#' @param x Input predictor weight vector
#' @return x Checked and refined input predictor weight vector
#' @keywords internal

dimFun_2C_1AIR = function(x){
  if(is.null(dim(x))){
    return(c(0,0))
  }else(return(dim(x)))
}

###### WeightsFun_2C_1AIR() ######

#' WeightsFun_2C_1AIR
#'
#' Support function, generates all possible weights for NBI subproblems
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weights All possible weights for NBI subproblem
#' @keywords internal

WeightsFun_2C_1AIR = function(n, k){

  # package env variables
  # weight, Weights, Formers, Layer, lastone, currentone
  #
  # Generates all possible weights for NBI subproblems given:
  # n, the number of objectives
  # 1/k, the uniform spacing between two w_i (k integral)
  # This is essentially all the possible integral partitions
  # of integer k into n parts.

  assign("WeightSub_2C_1AIR", matrix(0,1,n), rMOSTenv_2C_1AIR)
  # WeightSub <<- matrix(0,1,n)
  assign("Weights_2C_1AIR", vector(), rMOSTenv_2C_1AIR)
  # Weights <<- vector()
  assign("Formers_2C_1AIR", vector(), rMOSTenv_2C_1AIR)
  # Formers <<- vector()
  assign("Layer_2C_1AIR", n, rMOSTenv_2C_1AIR)
  # Layer <<- n
  assign("lastone_2C_1AIR", -1, rMOSTenv_2C_1AIR)
  # lastone <<- -1
  assign("currentone_2C_1AIR", -1, rMOSTenv_2C_1AIR)
  # currentone <<- -1

  Weight_Generate_2C_1AIR(1, k)

  return(list(Weights = rMOSTenv_2C_1AIR$Weights_2C_1AIR, Formers = rMOSTenv_2C_1AIR$Formers_2C_1AIR))

}

###### Weight_Generate_2C_1AIR() ######

#' Weight_Generate_2C_1AIR
#'
#' Function intended to test the weight generation scheme for NBI for > 2 objectives
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weight_Generate_2C_1AIR
#' @keywords internal

Weight_Generate_2C_1AIR = function(n, k){

  # package env variables:
  # weight Weights Formers Layer lastone currentone

  # wtgener_test(n,k)
  #
  # Intended to test the weight generation scheme for NBI for > 2 objectives
  # n is the number of objectives
  # 1/k is the uniform spacing between two w_i (k integral)

  # if(n == Layer){
  #
  #   if(currentone >= 0){
  #     Formers <<- c(Formers,lastone)
  #     lastone <<- currentone
  #     currentone <<- -1
  #   }else{
  #     num = dimFun_2C_1AIR(Weights)[2]
  #     Formers <<- c(Formers,num)
  #   }
  #
  #   WeightSub[(Layer - n + 1)] <<- k
  #   Weights <<- cbind(Weights,t(WeightSub))
  #
  # }else{
  #
  #   for(i in 0:k){
  #     if(n == (Layer - 2)){
  #       num = dimFun_2C_1AIR(Weights)[2]
  #       currentone <<- num+1
  #     }
  #
  #     WeightSub[(Layer - n + 1)] <<- i
  #     Weight_Generate_2C_1AIR(n+1, k-i)
  #   }
  #
  # }

  if(n == rMOSTenv_2C_1AIR$Layer_2C_1AIR){

    if(rMOSTenv_2C_1AIR$currentone_2C_1AIR >= 0){
      rMOSTenv_2C_1AIR$Formers_2C_1AIR <- c(rMOSTenv_2C_1AIR$Formers_2C_1AIR,rMOSTenv_2C_1AIR$lastone_2C_1AIR)
      rMOSTenv_2C_1AIR$lastone_2C_1AIR <- rMOSTenv_2C_1AIR$currentone_2C_1AIR
      rMOSTenv_2C_1AIR$currentone_2C_1AIR <- -1

    }else{
      num = dimFun_2C_1AIR(rMOSTenv_2C_1AIR$Weights_2C_1AIR)[2]
      rMOSTenv_2C_1AIR$Formers_2C_1AIR <- c(rMOSTenv_2C_1AIR$Formers_2C_1AIR,num)
    }

    rMOSTenv_2C_1AIR$WeightSub_2C_1AIR[(rMOSTenv_2C_1AIR$Layer_2C_1AIR - n + 1)] <- k
    rMOSTenv_2C_1AIR$Weights_2C_1AIR <- cbind(rMOSTenv_2C_1AIR$Weights_2C_1AIR,t(rMOSTenv_2C_1AIR$WeightSub_2C_1AIR))

  }else{

    for(i in 0:k){
      if(n == (rMOSTenv_2C_1AIR$Layer_2C_1AIR - 2)){
        num = dimFun_2C_1AIR(rMOSTenv_2C_1AIR$Weights_2C_1AIR)[2]
        rMOSTenv_2C_1AIR$currentone_2C_1AIR <- num+1
      }

      rMOSTenv_2C_1AIR$WeightSub_2C_1AIR[(rMOSTenv_2C_1AIR$Layer_2C_1AIR - n + 1)] <- i
      Weight_Generate_2C_1AIR(n+1, k-i)
    }

  }

}

###### myLinCom_2C_1AIR() ######

#' myLinCom_2C_1AIR
#'
#' Support function
#' @param x Input predictor weight vector
#' @return f Criterion vector
#' @keywords internal

myLinCom_2C_1AIR = function(x){

  # package env variable: g_Weight_2C_1AIR
  F = myFM_2C_1AIR(x)
  f = t(rMOSTenv_2C_1AIR$g_Weight_2C_1AIR)%*%F

  return(f)

}

###### myT_2C_1AIR() ######

#' myT_2C_1AIR
#'
#' Support function, define criterion space for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return f Temporary criterion space
#' @keywords internal

myT_2C_1AIR = function(x_t){

  f = x_t[length(x_t)]
  return(f)

}

###### myTCon_eq_2C_1AIR() ######

#' myTCon_eq_2C_1AIR
#'
#' Support function, define constraint condition for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return ceq Temporary constraint condition
#' @keywords internal

myTCon_eq_2C_1AIR = function(x_t){

  # package env variables:
  # g_Normal_2C_1AIR g_StartF_2C_1AIR

  t = x_t[length(x_t)]
  x = x_t[1:(length(x_t)-1)]

  fe  = -myFM_2C_1AIR(x) - rMOSTenv_2C_1AIR$g_StartF_2C_1AIR - t * rMOSTenv_2C_1AIR$g_Normal_2C_1AIR

  ceq1 = myCon_eq_2C_1AIR(x)
  ceq = c(ceq1,fe)

  return(ceq)

}
