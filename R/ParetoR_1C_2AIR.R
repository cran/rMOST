# Command Function (ParetoR) - Pareto-Optimization via Normal Boundary Intersection for 1 continuous objective and 2 adverse impact objectives
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com


#' ParetoR_1C_2AIR
#'
#' Command function to optimize 1 non-adverse impact objective and 2 adverse
#' impact objectives via NBI algorithm
#' @param sr Selection ratio in the full applicant pool
#' @param prop1 Proportion of minority1 applicants in the full applicant pool
#' @param prop2 Proportion of minority2 applicants in the full applicant pool
#' @param Rx Matrix with intercorrelations among predictors
#' @param Rxy1 Vector with correlation between each predictor and the non-adverse
#'   impact objective
#' @param d1 Subgroup difference 1; standardized mean differences between minority1
#'   and majority subgroups on each predictor in full applicant pool
#' @param d2 Subgroup difference 2; standardized mean differences between minority2
#'   and majority subgroups on each predictor in full applicant pool
#' @param Spac Number of solutions
#' @return out Pareto-Optimal solution with objective outcome values (Criterion) and
#'   predictor weights (ParetoWeights)
#' @export

ParetoR_1C_2AIR = function(sr, prop1, prop2, Rx, Rxy1, d1, d2,
                           # Spac = 20){ # SHINY app says 20
                           Spac = 10){

  # sr_1C_2AIR <<- sr
  # prop1_1C_2AIR <<- prop1
  # prop2_1C_2AIR <<- prop2
  # Rx_1C_2AIR <<- Rx
  # Rxy1_1C_2AIR <<- Rxy1
  # d1_1C_2AIR <<- d1
  # d2_1C_2AIR <<- d2
  assign("sr_1C_2AIR", sr, rMOSTenv_1C_2AIR)
  assign("prop1_1C_2AIR", prop1, rMOSTenv_1C_2AIR)
  assign("prop2_1C_2AIR", prop2, rMOSTenv_1C_2AIR)
  assign("Rx_1C_2AIR", Rx, rMOSTenv_1C_2AIR)
  assign("Rxy1_1C_2AIR", Rxy1, rMOSTenv_1C_2AIR)
  assign("d1_1C_2AIR", d1, rMOSTenv_1C_2AIR)
  assign("d2_1C_2AIR", d2, rMOSTenv_1C_2AIR)

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

  out = suppressWarnings(NBI_1C_2AIR(X0, Spac, Fnum, VLB, VUB, TolX, TolF, TolCon))

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

NBI_1C_2AIR = function(X0,Spac,Fnum,VLB=vector(),VUB=vector(),TolX=1e-7,TolF=1e-7,TolCon=1e-7){

  # cat('\n Estimating Pareto-Optimal Solution ... \n')

  #------------------------------Initialize Options-------------------------------#

  X0 = assert_col_vec_1C_2AIR(X0)
  VLB = assert_col_vec_1C_2AIR(VLB)
  VUB = assert_col_vec_1C_2AIR(VUB)

  # Number of variables
  nvars = length(X0)

  opts = list("maxeval" = (nvars+1)*1000
              ,"xtol_rel" = TolX
              ,"ftol_rel" = TolF
  )

  #Initialize PHI

  PHI = matrix(0,Fnum,Fnum)

  #------------------------------Shadow Point-------------------------------#

  # cat('\n ----Step 1: find shadow minimum---- \n')

  ShadowF = matrix(0,Fnum)
  ShadowX = matrix(0,nvars,Fnum)
  out = WeightsFun_1C_2AIR(Fnum,Spac)
  Weight = out$Weights
  Near = out$Formers
  rm(out)
  Weight = Weight/Spac

  for(i in 1:Fnum){

    temp = c(1,(Spac+1),dim(Weight)[2])
    j = temp[i]
    # g_Weight_1C_2AIR <<- Weight[,j]
    assign("g_Weight_1C_2AIR", Weight[,j], rMOSTenv_1C_2AIR)

    out = suppressMessages(nloptr::slsqp(x0 = X0, fn = myLinCom_1C_2AIR
                                         ,lower = VLB, upper = VUB
                                         ,hin = myCon_ineq_1C_2AIR
                                         ,heq = myCon_eq_1C_2AIR,
                                         control = list("maxeval" = (nvars+1)*1000,
                                                        "xtol_rel" = TolX,
                                                        "ftol_rel" = TolF)
    ))

    x = out$par

    # ########
    # print(paste0("i = ", i))
    # print("out$par:")
    # print(round(out$par, 3))
    # print("out$value: ")
    # print(round(out$value, 3))
    # ########

    x[1:nvars] = x[1:nvars]/sum(x[1:nvars])
    x[2:nvars] = x[2:nvars]/sum(x[2:nvars])

    ShadowX[,i] = x

    tempF = myFM_1C_2AIR(x)
    ShadowF[i] = tempF[i]


  }

  #------------------------------Matrix PHI-------------------------------#

  # cat('\n ----Step 2: find PHI---- \n')

  for(i in 1:Fnum){

    PHI[,i] = myFM_1C_2AIR(ShadowX[,i]) - ShadowF
    PHI[i,i] = 0

  }

  # ShadowF_1C_2AIR <<- ShadowF
  # PHI_1C_2AIR <<- PHI
  assign("ShadowF_1C_2AIR", ShadowF, rMOSTenv_1C_2AIR)
  assign("PHI_1C_2AIR", PHI, rMOSTenv_1C_2AIR)

  #------------------------------Quasi-Normal Direction-------------------------------#

  # cat('\n ----Step 3: find Quasi-Normal---- \n')

  assign("g_Normal_1C_2AIR", -rMOSTenv_1C_2AIR$PHI_1C_2AIR%*%matrix(1,Fnum,1), rMOSTenv_1C_2AIR)
  assign("g_Normal_1C_2AIR", rMOSTenv_1C_2AIR$g_Normal_1C_2AIR/(norm(matrix(rMOSTenv_1C_2AIR$g_Normal_1C_2AIR),"2")^2), rMOSTenv_1C_2AIR) # added based on try_NBI3D

  #------------------------------weights-------------------------------#

  # cat('\n ----Step 4: create weights---- \n')

  out = WeightsFun_1C_2AIR(Fnum,Spac)
  Weight = out$Weight
  Near = out$Formers
  Weight = Weight/Spac
  num_w = dimFun_1C_2AIR(Weight)[2]

  #------------------------------NBI Subproblems-------------------------------#


  # Starting point for first NBI subproblem is the minimizer of f_1(x)
  xstart = c(ShadowX[,1],0) # changed based on try_NBI3D

  Pareto_Fmat = vector()       # Pareto Optima in F-space
  Pareto_Xmat = vector()       # Pareto Optima in X-space
  X_Near      = vector()

  # solve NBI subproblems
  for(k in 1:num_w){

    w  = Weight[,k]

    if(Near[k] > 0){

      xstart = X_Near[,Near[k]]
      # start X is the previous weight-order's X

    }

    # g_StartF_1C_2AIR <<- PHI_1C_2AIR%*%w + ShadowF_1C_2AIR
    assign("g_StartF_1C_2AIR", rMOSTenv_1C_2AIR$PHI_1C_2AIR%*%w + rMOSTenv_1C_2AIR$ShadowF_1C_2AIR, rMOSTenv_1C_2AIR)

    # SOLVE NBI SUBPROBLEM

    # ##########
    # print("")
    # print("SOLVE NBI SUBPROBLEM")
    # print(paste0("k = ", k))
    # ##########

    # if(k %in% c(1,(Spac+1),num_w)){
    if(k %in% c(1,(Spac+1))){

      VLB_trials = c(VLB,0)

      out = suppressMessages(nloptr::slsqp(x0 = xstart, fn = myT_1C_2AIR
                                           ,lower = VLB_trials
                                           ,upper = c(VUB,Inf)
                                           ,hin = myCon_ineq_1C_2AIR
                                           ,heq = myTCon_eq_1_1C_2AIR
                                           ,control = opts
      ))

      # ########
      # print("First IF statement")
      # print(round(out$par, 3))
      # print(round(out$value, 3))
      # ########

    }

    if(k %in% num_w){

      VLB_trials = c(VLB,0)

      out = suppressMessages(nloptr::slsqp(x0 = xstart, fn = myT_1C_2AIR
                                           ,lower = VLB_trials
                                           ,upper = c(VUB,Inf)
                                           ,hin = myCon_ineq_1C_2AIR
                                           ,heq = myTCon_eq_2_1C_2AIR
                                           ,control = opts
      ))

      # ########
      # print("Second IF statement")
      # print(round(out$par, 3))
      # print(round(out$value, 3))
      # ########

    }

    if(!(k %in% c(1,(Spac+1),num_w))){

      VLB_trials = c(VLB,-Inf)

      out = suppressMessages(nloptr::slsqp(x0 = xstart, fn = myT_1C_2AIR
                                           ,lower = VLB_trials
                                           ,upper = c(VUB,Inf)
                                           ,hin = myCon_ineq_1C_2AIR
                                           ,heq = myTCon_eq_1C_2AIR
                                           ,control = opts
      ))

      # ########
      # print("Third IF statement")
      # print(round(out$par, 3))
      # print(round(out$value, 3))
      # ########

    }

    # ########
    # print("x_trial")
    # print(round(out$par, 3))
    # print(paste0("nvars = ", nvars))
    # ########

    x_trial = out$par
    rm(out)

    x_trial[1:nvars] = x_trial[1:nvars]/sum(x_trial[1:nvars])
    # ########
    # print(round(x_trial, 3))
    # ########
    x_trial[2:nvars] = x_trial[2:nvars]/sum(x_trial[2:nvars])
    # ########
    # print(round(x_trial, 3))
    # ########

    Pareto_Xmat = cbind(Pareto_Xmat, x_trial[1:nvars])        # Pareto optima in X-space
    Pareto_Fmat = cbind(Pareto_Fmat, -myFM_1C_2AIR(x_trial[1:nvars]))  # Pareto optima in F-space
    X_Near = cbind(X_Near,x_trial)
    # ###########
    # print(round(Pareto_Xmat, 3))
    # print(round(Pareto_Fmat, 3))
    # ###########

  }

  # ###########
  # print(round(X_Near, 3))
  # ###########

  #------------------------------Plot Solutions-------------------------------#


  # if(graph==TRUE){plotPareto(t(Pareto_Fmat))}

  #------------------------------Output Solutions-------------------------------#

  #   Output Solution

  Pareto_Fmat = t(Pareto_Fmat)
  Pareto_Xmat = t(Pareto_Xmat)[,-1]
  # ###########
  # print(round(Pareto_Xmat, 3))
  # print(round(Pareto_Fmat, 3))
  # ###########
  colnames(Pareto_Fmat) = c("Ry","AIR1","AIR2")
  colnames(Pareto_Xmat) = c(paste0(rep("P",(nvars-1)),1:(nvars-1)))

  # if(display_solution == TRUE){
  #
  #   solution = round(cbind(Pareto_Fmat,Pareto_Xmat),3)
  #   colnames(solution) = c("Ry","AIR1","AIR2",paste0(rep("P",(nvars-1)),1:(nvars-1)))
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

# Function List

## myFM_1C_2AIR
## myCon_ineq_1C_2AIR
## myCon_eq_1C_2AIR
## assert_col_vec_1C_2AIR
## dimFun_1C_2AIR
## WeightsFun_1C_2AIR
## Weight_Generate_1C_2AIR
## myLinCom_1C_2AIR
## myT_1C_2AIR
## myTCon_eq_1C_2AIR
## plotPareto_1C_2AIR

###### combR_1C_2AIR()######

#' combR_1C_2AIR
#'
#' Support function to create predictor-criterion matrix
#' @param Rx Predictor inter-correlation matrix
#' @param Ry Predictor-criterion correlation (validity)
#' @return Rxy Predictor-criterion correlation matrix
#' @keywords internal

combR_1C_2AIR <- function(Rx, Ry){
  cbind(rbind(Rx, c(Ry)), c(Ry, 1))
}

# ###### myFS() ######
#
# #' myFS
# #'
# #' Supporting function, defines criterion space
# #' @param x Input predictor weight vector
# #' @importFrom stats pnorm
# #' @return f Criterion vector
#
# myFS = function(x){
#
#   # b = x
#   b = x[-1] # Predictor weight
#   p_c = x[1] # Cut-off score
#
#   # Obtain within-package 'global' variables from package env
#   # d1 <- d1_1C_2AIR
#   # d2 <- d2_1C_2AIR
#   # R <- combR_1C_2AIR(Rx_1C_2AIR, Rxy1_1C_2AIR)
#   # R_u <- Rx_1C_2AIR
#   d1 <- rMOSTenv_1C_2AIR$d1_1C_2AIR
#   d2 <- rMOSTenv_1C_2AIR$d2_1C_2AIR
#   R <- combR_1C_2AIR(rMOSTenv_1C_2AIR$Rx_1C_2AIR, rMOSTenv_1C_2AIR$Rxy1_1C_2AIR)
#   R_u <- rMOSTenv_1C_2AIR$Rx_1C_2AIR
#
#   # variance of minority and majority applicant weighted predictor
#   # composite (P) distribution (DeCorte, 1999)
#   sigma_p = sqrt(t(b)%*%R_u%*%b)
#
#   # mean of Minority_1 weighted predictor composite distribution (DeCorte, 1999)
#   p_i_bar1 = d2%*%b/sigma_p
#   # mean of Minority_2 weighted predictor composite distribution (DeCorte, 1999)
#   p_i_bar2 = d1%*%b/sigma_p
#   # mean of Majority weighted predictor composite distribution (DeCorte, 1999)
#   p_a_bar0 = (d1+d2)%*%b/sigma_p
#
#   # Minority_1 selection ratio (denoted as h_i in DeCorte et al., 1999)
#   SR1 = 1 - pnorm(p_c, p_i_bar1)
#   # Minority_2 selection ratio (denoted as h_i in DeCorte et al., 1999)
#   SR2 = 1 - pnorm(p_c, p_i_bar2)
#   # Majority group selection ratio (denoted as h_i in DeCorte et al., 1999)
#   SR0 = 1 - pnorm(p_c, p_a_bar0)
#
#   # AIratio a_g (DeCorte et al., 2007)
#   a_g1 = SR1/SR0 # AI ratio of Minority_1
#   a_g2 = SR2/SR0 # AI ratio of Minority_2
#
#   # Composite Validity R_xy
#   Ry = t(c(t(b),0)%*%R%*%c(t(matrix(0,dimFun_1C_2AIR(R_u)[1],1)),1))/sqrt(t(b)%*%R_u%*%b) # DeCorte et al., 2007
#
#   if(g_Index == 1){f = -Ry}
#   if(g_Index == 2){f = -a_g1}
#   if(g_Index == 3){f = -a_g2}
#
#   return(f)
#
# }

###### myFM_1C_2AIR() ######

#' myFM_1C_2AIR
#'
#' Supporting function, defines criterion space
#' @param x Input predictor weight vector
#' @importFrom stats pnorm
#' @return f Criterion vector
#' @keywords internal

myFM_1C_2AIR = function(x){

  # b = x
  b = x[-1] # Predictor weight
  p_c = x[1] # Cutoff score

  # Obtain within-package 'global' variables from package env
  # d1 <- d1_1C_2AIR
  # d2 <- d2_1C_2AIR
  # R <- combR_1C_2AIR(Rx_1C_2AIR, Rxy1_1C_2AIR)
  # R_u <- Rx_1C_2AIR
  d1 <- rMOSTenv_1C_2AIR$d1_1C_2AIR
  d2 <- rMOSTenv_1C_2AIR$d2_1C_2AIR
  R <- combR_1C_2AIR(rMOSTenv_1C_2AIR$Rx_1C_2AIR, rMOSTenv_1C_2AIR$Rxy1_1C_2AIR)
  R_u <- rMOSTenv_1C_2AIR$Rx_1C_2AIR

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%R_u%*%b)

  # mean of Minority_1 weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar1 = d2%*%b/sigma_p
  # mean of Minority_2 weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar2 = d1%*%b/sigma_p
  # mean of Majority weighted predictor composite distribution (DeCorte, 1999)
  p_a_bar0 = (d1+d2)%*%b/sigma_p

  # Minority_1 selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR1 = 1 - pnorm(p_c, p_i_bar1)
  # Minority_2 selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR2 = 1 - pnorm(p_c, p_i_bar2)
  # Majority group selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR0 = 1 - pnorm(p_c, p_a_bar0)

  # AIratio a_g (DeCorte et al., 2007)
  a_g1 = SR1/SR0 # AI ratio of Minority_1
  a_g2 = SR2/SR0 # AI ratio of Minority_2

  # Composite Validity R_xy
  Ry = t(c(t(b),0)%*%R%*%c(t(matrix(0,dimFun_1C_2AIR(R_u)[1],1)),1))/sqrt(t(b)%*%R_u%*%b) # DeCorte et al., 2007

  f = matrix(1,3,1)
  f[1,] = -Ry
  f[2,] = -a_g1
  f[3,] = -a_g2

  return(f)

}

####### myCon_ineq_1C_2AIR() ######

# Nonlinear inequalities at x

#' myCon_ineq_1C_2AIR
#'
#' Support function, defines inequal constraint condition
#' @param x Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
#' @keywords internal

myCon_ineq_1C_2AIR = function(x){return(vector())}

####### myCon_eq_1C_2AIR() ######

# Nonlinear equalities at x

#' myCon_eq_1C_2AIR
#'
#' Support function, defines equal constraint condition
#' @param x Input predictor weight vector
#' @importFrom stats pnorm
#' @return Equal constraint condition for use in NBI()
#' @keywords internal

myCon_eq_1C_2AIR = function(x){

  # Obtain within-package 'global' variable from package env
  # sr <- sr_1C_2AIR
  # prop1 <- prop1_1C_2AIR
  # prop2 <- prop2_1C_2AIR
  # d1 <- d1_1C_2AIR
  # d2 <- d2_1C_2AIR
  # R_u <- Rx_1C_2AIR
  sr <- rMOSTenv_1C_2AIR$sr_1C_2AIR
  prop1 <- rMOSTenv_1C_2AIR$prop1_1C_2AIR
  prop2 <- rMOSTenv_1C_2AIR$prop2_1C_2AIR
  d1 <- rMOSTenv_1C_2AIR$d1_1C_2AIR
  d2 <- rMOSTenv_1C_2AIR$d2_1C_2AIR
  R_u <- rMOSTenv_1C_2AIR$Rx_1C_2AIR

  b <- x[-1] # Predictor weight
  p_c <- x[1] # Cutoff score

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%R_u%*%b)

  # mean of Minority_1 weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar1 = d2%*%b/sigma_p
  # mean of Minority_2 weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar2 = d1%*%b/sigma_p
  # mean of Majority weighted predictor composite distribution (DeCorte, 1999)
  p_a_bar0 = (d1+d2)%*%b/sigma_p

  # Black selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR1 = 1 - pnorm(p_c, p_i_bar1)
  # Hispanic selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR2 = 1 - pnorm(p_c, p_i_bar2)
  # White group selection ratio (denoted as h_a in DeCorte et al., 1999)
  SR0 = 1 - pnorm(p_c, p_a_bar0)

  # Nonlinear equalities at x

  ceq = matrix(1,2,1)
  ceq[1,] = SR1*prop1 + SR2*prop2 + SR0*(1 - prop1 - prop2) - sr # DeCorte et al. (2007)
  ceq[2,] = (t(b)%*%R_u%*%b) - 1

  return(ceq)

}

####### myCon_eq_1_1C_2AIR() ######

# Equality constraint for job performance validity
# Equality constraint for AI ratio of minority group 1

# Nonlinear equalities at x

#' myCon_eq_1_1C_2AIR
#'
#' Support function, defines equal constraint condition
#' @param x Input predictor weight vector
#' @importFrom stats pnorm
#' @return Equal constraint condition for use in NBI()
#' @keywords internal

myCon_eq_1_1C_2AIR = function(x){

  # Obtain within-package 'global' variable from package env
  # sr <- sr_1C_2AIR
  # prop1 <- prop1_1C_2AIR
  # d1 <- d1_1C_2AIR
  # R_u <- Rx_1C_2AIR
  sr <- rMOSTenv_1C_2AIR$sr_1C_2AIR
  prop1 <- rMOSTenv_1C_2AIR$prop1_1C_2AIR
  d1 <- rMOSTenv_1C_2AIR$d1_1C_2AIR
  R_u <- rMOSTenv_1C_2AIR$Rx_1C_2AIR

  b <- x[-1] # Predictor weight
  p_c <- x[1] # Cutoff score

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%R_u%*%b)

  # mean of minority weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar1 = 0
  # mean of Majority weighted predictor composite distribution (DeCorte, 1999)
  p_a_bar0 = d1%*%b/sigma_p

  # Black selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR1 = 1 - pnorm(p_c, p_i_bar1)
  # White group selection ratio (denoted as h_a in DeCorte et al., 1999)
  SR0 = 1 - pnorm(p_c, p_a_bar0)

  # Nonlinear equalities at x

  ceq = matrix(1,2,1)
  ceq[1,] = SR1*prop1 + SR0*(1 - prop1) - sr # DeCorte et al. (2007)
  ceq[2,] = (t(b)%*%R_u%*%b) - 1

  return(ceq)

}

####### myCon_eq_2_1C_2AIR() ######

# Equality constraint for AI ratio of minority group 2

# Nonlinear equalities at x

#' myCon_eq_2_1C_2AIR
#'
#' Support function, defines equal constraint condition
#' @param x Input predictor weight vector
#' @importFrom stats pnorm
#' @return Equal constraint condition for use in NBI()
#' @keywords internal

myCon_eq_2_1C_2AIR = function(x){

  # Obtain within-package 'global' variable from package env
  # sr <- sr_1C_2AIR
  # prop2 <- prop1_1C_2AIR
  # d2 <- d1_1C_2AIR
  # R_u <- Rx_1C_2AIR
  sr <- rMOSTenv_1C_2AIR$sr_1C_2AIR
  prop2 <- rMOSTenv_1C_2AIR$prop2_1C_2AIR ####REVISED#####
  d2 <- rMOSTenv_1C_2AIR$d2_1C_2AIR ####REVISED#####
  R_u <- rMOSTenv_1C_2AIR$Rx_1C_2AIR

  b <- x[-1] # Predictor weight
  p_c <- x[1] # Cutoff score

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%R_u%*%b)

  # mean of minority weighted predictor composite distribution (DeCorte, 1999)
  p_i_bar2 = 0
  # mean of Majority weighted predictor composite distribution (DeCorte, 1999)
  p_a_bar0 = d2%*%b/sigma_p

  # Black selection ratio (denoted as h_i in DeCorte et al., 1999)
  SR2 = 1 - pnorm(p_c, p_i_bar2)
  # White group selection ratio (denoted as h_a in DeCorte et al., 1999)
  SR0 = 1 - pnorm(p_c, p_a_bar0)

  # Nonlinear equalities at x

  ceq = matrix(1,2,1)
  ceq[1,] = SR2*prop2 + SR0*(1 - prop2) - sr # DeCorte et al. (2007)
  ceq[2,] = (t(b)%*%R_u%*%b) - 1

  return(ceq)

}

###### assert_col_vec_1C_2AIR() ######

#' assert_col_vec_1C_2AIR
#'
#' Support function, refines intermediate variable for use in NBI()
#' @param v Intermediate variable v
#' @return Refined variable v
#' @keywords internal

assert_col_vec_1C_2AIR = function(v){
  if(is.null(dimFun_1C_2AIR(v))){
    v=v
  }else if(dimFun_1C_2AIR(v)[1] < dimFun_1C_2AIR(v)[2]){v = t(t)}
  return(v)}

###### dimFun_1C_2AIR() ######

#' dimFun_1C_2AIR
#'
#' Support function, checks input predictor weight vector x
#' @param x Input predictor weight vector
#' @return x Checked and refined input predictor weight vector
#' @keywords internal

dimFun_1C_2AIR = function(x){
  if(is.null(dim(x))){
    return(c(0,0))
  }else(return(dim(x)))
}

###### WeightsFun_1C_2AIR() ######

#' WeightsFun_1C_2AIR
#'
#' Support function, generates all possible weights for NBI subproblems
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weights All possible weights for NBI subproblem
#' @keywords internal

WeightsFun_1C_2AIR = function(n, k){

  # package env variables
  # weight, Weights, Formers, Layer, lastone, currentone
  #
  # Generates all possible weights for NBI subproblems given:
  # n, the number of objectives
  # 1/k, the uniform spacing between two w_i (k integral)
  # This is essentially all the possible integral partitions
  # of integer k into n parts.

  assign("WeightSub_1C_2AIR", matrix(0,1,n), rMOSTenv_1C_2AIR)
  # WeightSub <<- matrix(0,1,n)
  assign("Weights_1C_2AIR", vector(), rMOSTenv_1C_2AIR)
  # Weights <<- vector()
  assign("Formers_1C_2AIR", vector(), rMOSTenv_1C_2AIR)
  # Formers <<- vector()
  assign("Layer_1C_2AIR", n, rMOSTenv_1C_2AIR)
  # Layer <<- n
  assign("lastone_1C_2AIR", -1, rMOSTenv_1C_2AIR)
  # lastone <<- -1
  assign("currentone_1C_2AIR", -1, rMOSTenv_1C_2AIR)
  # currentone <<- -1

  Weight_Generate_1C_2AIR(1, k)

  return(list(Weights = rMOSTenv_1C_2AIR$Weights_1C_2AIR, Formers = rMOSTenv_1C_2AIR$Formers_1C_2AIR))

}

###### Weight_Generate_1C_2AIR() ######

#' Weight_Generate_1C_2AIR
#'
#' Function intended to test the weight generation scheme for NBI for > 2 objectives
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weight_Generate_1C_2AIR
#' @keywords internal

Weight_Generate_1C_2AIR = function(n, k){

  # package env variables:
  # weight_1C_2AIR Weights_1C_2AIR Formers_1C_2AIR Layer_1C_2AIR lastone_1C_2AIR currentone_1C_2AIR

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
  #     num = dimFun_1C_2AIR(Weights)[2]
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
  #       num = dimFun_1C_2AIR(Weights)[2]
  #       currentone <<- num+1
  #     }
  #
  #     WeightSub[(Layer - n + 1)] <<- i
  #     Weight_Generate_1C_2AIR(n+1, k-i)
  #   }
  #
  # }
  if(n == rMOSTenv_1C_2AIR$Layer_1C_2AIR){

    if(rMOSTenv_1C_2AIR$currentone_1C_2AIR >= 0){
      rMOSTenv_1C_2AIR$Formers_1C_2AIR <- c(rMOSTenv_1C_2AIR$Formers_1C_2AIR,rMOSTenv_1C_2AIR$lastone_1C_2AIR)
      rMOSTenv_1C_2AIR$lastone_1C_2AIR <- rMOSTenv_1C_2AIR$currentone_1C_2AIR
      rMOSTenv_1C_2AIR$currentone_1C_2AIR <- -1
    }else{
      num = dimFun_1C_2AIR(rMOSTenv_1C_2AIR$Weights_1C_2AIR)[2]
      rMOSTenv_1C_2AIR$Formers_1C_2AIR <- c(rMOSTenv_1C_2AIR$Formers_1C_2AIR,num)
    }

    rMOSTenv_1C_2AIR$WeightSub_1C_2AIR[(rMOSTenv_1C_2AIR$Layer_1C_2AIR - n + 1)] <- k
    rMOSTenv_1C_2AIR$Weights_1C_2AIR <- cbind(rMOSTenv_1C_2AIR$Weights_1C_2AIR,t(rMOSTenv_1C_2AIR$WeightSub_1C_2AIR))

  }else{

    for(i in 0:k){
      if(n == (rMOSTenv_1C_2AIR$Layer_1C_2AIR - 2)){
        num = dimFun_1C_2AIR(rMOSTenv_1C_2AIR$Weights_1C_2AIR)[2]
        rMOSTenv_1C_2AIR$currentone_1C_2AIR <- num+1
      }

      rMOSTenv_1C_2AIR$WeightSub_1C_2AIR[(rMOSTenv_1C_2AIR$Layer_1C_2AIR - n + 1)] <- i
      Weight_Generate_1C_2AIR(n+1, k-i)
    }

  }

}

###### myLinCom_1C_2AIR() ######

#' myLinCom_1C_2AIR
#'
#' Support function
#' @param x Input predictor weight vector
#' @return f Criterion vector
#' @keywords internal

myLinCom_1C_2AIR = function(x){

  # package env variable: g_Weight_1C_2AIR
  F   = myFM_1C_2AIR(x)
  f = t(rMOSTenv_1C_2AIR$g_Weight_1C_2AIR)%*%F

  return(f)

}

###### myT_1C_2AIR() ######

#' myT_1C_2AIR
#'
#' Support function, define criterion space for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return f Temporary criterion space
#' @keywords internal

myT_1C_2AIR = function(x_t){

  # f = x_t[length(x_t)]
  f = -x_t[length(x_t)]

  return(f)

}

###### myTCon_eq_1C_2AIR() ######

#' myTCon_eq_1C_2AIR
#'
#' Support function, define constraint condition for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return ceq Temporary constraint condition
#' @keywords internal

myTCon_eq_1C_2AIR = function(x_t){

  # package env variables:
  # g_Normal_1C_2AIR g_StartF_1C_2AIR

  t = x_t[length(x_t)]
  x = x_t[1:(length(x_t)-1)]

  fe  = myFM_1C_2AIR(x) - rMOSTenv_1C_2AIR$g_StartF_1C_2AIR - t * rMOSTenv_1C_2AIR$g_Normal_1C_2AIR

  ceq1 = myCon_eq_1C_2AIR(x)
  ceq = c(ceq1,fe)

  return(ceq)

}

###### myTCon_eq_1_1C_2AIR() ######

# myTCon_eq_1C_2AIR() constraint for:
# - the first criterion: job performance validity;
# - second criterion: AI ratio of minority group 1

#' myTCon_eq_1_1C_2AIR
#'
#' Support function, define constraint condition for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return ceq Temporary constraint condition
#' @keywords internal

myTCon_eq_1_1C_2AIR = function(x_t){

  # package env variables:
  # g_Normal_1C_2AIR g_StartF_1C_2AIR

  t = x_t[length(x_t)]
  x = x_t[1:(length(x_t)-1)]

  fe  = myFM_1C_2AIR(x) - rMOSTenv_1C_2AIR$g_StartF_1C_2AIR - t * rMOSTenv_1C_2AIR$g_Normal_1C_2AIR

  ceq1 = myCon_eq_1_1C_2AIR(x)
  ceq = c(ceq1,fe)

  return(ceq)

}

###### myTCon_eq_2_1C_2AIR() ######

# myTCon_eq() constraint for:
# - second criterion: AI ratio of minority group 2

#' myTCon_eq_2_1C_2AIR
#'
#' Support function, define constraint condition for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return ceq Temporary constraint condition
#' @keywords internal

myTCon_eq_2_1C_2AIR = function(x_t){

  # package env variables:
  # g_Normal_1C_2AIR g_StartF_1C_2AIR

  t = x_t[length(x_t)]
  x = x_t[1:(length(x_t)-1)]

  fe  = myFM_1C_2AIR(x) - rMOSTenv_1C_2AIR$g_StartF_1C_2AIR - t * rMOSTenv_1C_2AIR$g_Normal_1C_2AIR

  ceq1 = myCon_eq_2_1C_2AIR(x)
  ceq = c(ceq1,fe)

  return(ceq)

}

####### myTCon_ineq_1C_2AIR() ######

# Nonlinear inequalities at x

#' myTCon_ineq_1C_2AIR
#'
#' Support function, defines inequal constraint condition
#' @param x_t Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
#' @keywords internal

myTCon_ineq_1C_2AIR = function(x_t){

  return(vector())

}
