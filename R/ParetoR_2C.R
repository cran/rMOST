# Pareto-Optimization via Normal Boundary Intersection: 2 continuous objectives
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com


#' ParetoR_2C
#'
#' Command function to optimize 2 non-adverse impact objectives via NBI
#' algorithm
#' @param Rx  Matrix with intercorrelations among predictors
#' @param Rxy1  Vector with correlation between each predictor and non-adverse impact objective 1
#' @param Rxy2  Vector with correlation between each predictor and non-adverse impact objective 2
#' @param Spac  Number of Pareto points
#' @param graph If TRUE, plots will be generated for Pareto-optimal curve and
#'   predictor weights
#' @return out Pareto-Optimal solution with objective outcome values (Criterion) and
#'   predictor weights (ParetoWeights)
#' @keywords internal

ParetoR_2C = function(Rx, Rxy1, Rxy2, Spac = 10, graph = TRUE){

  # Rx_2C <<- Rx
  # Rxy1_2C <<- Rxy1
  # Rxy2_2C <<- Rxy2
  assign("Rx_2C", Rx, rMOSTenv)
  assign("Rxy1_2C", Rxy1, rMOSTenv)
  assign("Rxy2_2C", Rxy2, rMOSTenv)

  # Tolerance Level for Algorithm
  TolCon 	= 1.0000e-6 # tolerance of constraint
  TolF 	= 1.0000e-6 # tolerance of objective function
  TolX 	= 1.0000e-7 # tolerance of predictor variable

  # Do not change
  Fnum 	= 2
  Xnum = dim(Rx)[1]
  VLB = c(rep(0,Xnum))
  VUB = c(rep(1,Xnum))
  X0 = c(rep(1/Xnum,Xnum))

  ###### Find Pareto-Optimal Solution ######

  out = suppressWarnings(NBI_2C(X0, Spac, Fnum, VLB, VUB, TolX, TolF, TolCon, graph=graph))
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
#' @param graph If TRUE, plots will be generated for Pareto-optimal curve and predictor weights
#' @import nloptr
#' @return Pareto-Optimal solutions
#' @keywords internal

NBI_2C = function(X0,Spac,Fnum,VLB=vector(),VUB=vector(),TolX=1e-4,TolF=1e-4,TolCon=1e-7,graph=TRUE){

  # cat('\n Estimating Pareto-Optimal Solution ... \n')

  #------------------------------Initialize Options-------------------------------#

  Spac = Spac * 2

  X0 = assert_col_vec(X0)
  VLB = assert_col_vec(VLB)
  VUB = assert_col_vec(VUB)

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
  # xstart  = X0
  out = WeightsFun(Fnum,Spac)
  Weight = out$Weights
  Near = out$Formers
  rm(out)
  Weight = Weight/Spac

  for(i in 1:Fnum){

    temp = c(1,dim(Weight)[2])
    j = temp[i]
    # g_Weight <<- Weight[,j]
    assign("g_Weight", Weight[,j], rMOSTenv)
    fmin = 9999

    out = suppressMessages(nloptr::slsqp(x0 = X0, fn = myLinCom_2C
                                         ,lower = VLB, upper = VUB
                                         ,hin = myCon_ineq_2C
                                         ,heq = myCon_eq_2C
    ))
    x = out$par

    som = 0

    for(k in 1:nvars){
      som = som + x[k]
    }

    for(k in 1:nvars){
      x[k] = x[k]/som
    }
    # to make sum of x = 1

    ShadowX[,i] = x
    ShadowX = round(ShadowX,4)

    tempF = -myFM_2C(x)
    ShadowF[i] = round(tempF[i],4)

  }

  #------------------------------Matrix PHI-------------------------------#

  # cat('\n ----Step 2: find PHI---- \n')

  for(i in 1:Fnum){

    PHI[,i] = myFM_2C(ShadowX[,i]) + ShadowF
    PHI[i,i] = 0

  }

  # ShadowF <<- ShadowF
  # PHI <<- PHI
  assign("ShadowF", ShadowF, rMOSTenv)
  assign("PHI", PHI, rMOSTenv)

  # print(round(PHI,3))

  #Check to make sure that QPP is n-1 dimensional
  if(rcond(rMOSTenv$PHI) < 1e-8){stop(' Phi matrix singular, aborting.')}

  #------------------------------Quasi-Normal Direction-------------------------------#

  # cat('\n ----Step 3: find Quasi-Normal---- \n')

  # g_Normal <<- -PHI%*%matrix(1,Fnum,1)
  assign("g_Normal", -rMOSTenv$PHI%*%matrix(1,Fnum,1), rMOSTenv)

  #------------------------------weights-------------------------------#

  # cat('\n ----Step 4: create weights---- \n')

  out = WeightsFun(Fnum,Spac)
  Weight = out$Weight
  Near = out$Formers
  Weight = Weight/Spac
  num_w = dimFun(Weight)[2]


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
      Pareto_Fmat = cbind(Pareto_Fmat, (-rMOSTenv$PHI[,indiv_fn_index] + rMOSTenv$ShadowF))
      Pareto_Xmat = cbind(Pareto_Xmat, ShadowX[,indiv_fn_index])
      X_Near = cbind(X_Near, c(ShadowX[,indiv_fn_index],0))

    }else{

      w = rev(w)

      if(Near[k] > 0){

        xstart = X_Near[,Near[k]]
        #start X is the previous weight-order's X

      }

      #start point in F-space
      # g_StartF <<- PHI%*%w + ShadowF
      assign("g_StartF", rMOSTenv$PHI%*%w + rMOSTenv$ShadowF, rMOSTenv)

      # SOLVE NBI SUBPROBLEM

      out = suppressMessages(nloptr::slsqp(x0 = xstart, fn = myT
                                           ,lower = c(VLB,-Inf)
                                           ,upper = c(VUB,Inf)
                                           ,hin = myCon_ineq_2C
                                           ,heq = myTCon_eq_2C))

      x_trial = out$par
      f = out$value
      rm(out)

      Pareto_Fmat = cbind(Pareto_Fmat, -myFM_2C(x_trial[1:nvars]))  # Pareto optima in F-space
      som = 0

      for(k in 1:nvars){som = som + x_trial[k]}

      for(k in 1:nvars){x_trial[k] = x_trial[k]/som}

      Pareto_Xmat = cbind(Pareto_Xmat, x_trial[1:nvars])        # Pareto optima in X-space
      X_Near = cbind(X_Near,x_trial)

    }

  }

  #------------------------------Plot Solutions-------------------------------#

  #   cat('\n ----Step 6: plot---- \n')

  if(graph==TRUE){plotPareto_2C(Pareto_Fmat, Pareto_Xmat)}

  #------------------------------Output Solutions-------------------------------#

  #   Output Solution

  Pareto_Fmat = t(Pareto_Fmat)
  Pareto_Xmat = t(Pareto_Xmat)
  colnames(Pareto_Fmat) = c("Y1","Y2")
  colnames(Pareto_Xmat) = c(paste0(rep("P",nvars),1:nvars))

  # if(display_solution == TRUE){
  #
  #   solution = round(cbind(Pareto_Fmat,Pareto_Xmat),3)
  #   colnames(solution) = c("Y1","Y2", paste0(rep("P",nvars),1:nvars))
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

# User-Defined Input for NBI.r - Pareto-Optimization via Normal Boundary Intersection

# Related functions:
# myFM_2C
# myCon

###### myFM_2C() ######

#' myFM_2C
#'
#' Supporting function, defines criterion space
#' @param x Input predictor weight vector
#' @return f Criterion vector
#' @keywords internal

myFM_2C = function(x){

  # Rx <- Rx_2C
  # Rxy1 <- Rxy1_2C
  # Rxy2 <- Rxy2_2C
  Rx <- rMOSTenv$Rx_2C
  Rxy1 <- rMOSTenv$Rxy1_2C
  Rxy2 <- rMOSTenv$Rxy2_2C

  b = x

  combR = function(Rx, Ry){
    cbind(rbind(Rx,c(Ry,1)),c(Ry,1))
  }

  # Composite Validity R_cy1, Rcy2
  R_cy1 = t(c(t(b),0)%*%combR(Rx,Rxy1)%*%c(t(matrix(0,dimFun(Rx)[1],1)),1))/sqrt(t(b)%*%Rx%*%b) # DeCorte et al., 2007
  R_cy2 = t(c(t(b),0)%*%combR(Rx,Rxy2)%*%c(t(matrix(0,dimFun(Rx)[1],1)),1))/sqrt(t(b)%*%Rx%*%b) # DeCorte et al., 2007

  f = matrix(1,2,1)
  f[1,] = -R_cy1
  f[2,] = -R_cy2

  return(f)

}


####### myCon_ineq_2C() ######

# Nonlinear inequalities at x

#' myCon_ineq_2C
#'
#' Support function, defines inequal constraint condition
#' @param x Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
#' @keywords internal

myCon_ineq_2C = function(x){return(vector())}

####### myCon_eq_2C() ######

# Nonlinear equalities at x

#' myCon_eq_2C
#'
#' Support function, defines equal constraint condition
#' @param x Input predictor weight vector
#' @return Equal constraint condition for use in NBI()
#' @keywords internal

myCon_eq_2C = function(x){

  Rx <- rMOSTenv$Rx_2C

  b = x

  # Nonlinear equalities at x

  ceq = (t(b)%*%Rx%*%b) - 1

  return(ceq)

}

########################### Supporting Functions (B) ########################

# Supplementary Functions for NBI.r - Pareto-Optimization via Normal Boundary Intersection

# Function List
## assert_col_vec
## dimFun
## WeightsFun
## Weight_Generate
## myLinCom_2C
## myT
## myTCon_eq_2C
## plotPareto_2C

###### assert_col_vec() ######

#' assert_col_vec
#'
#' Support function, refines intermediate variable for use in NBI()
#' @param v Intermediate variable v
#' @return Refined variable v
#' @keywords internal

assert_col_vec = function(v){
  if(is.null(dimFun(v))){
    v=v
  }else if(dimFun(v)[1] < dimFun(v)[2]){v = t(t)}
  return(v)}

###### dimFun() ######

#' dimFun
#'
#' Support function, checks input predictor weight vector x
#' @param x Input predictor weight vector
#' @return x Checked and refined input predictor weight vector
#' @keywords internal

dimFun = function(x){
  if(is.null(dim(x))){
    return(c(0,0))
  }else(return(dim(x)))
}

###### WeightsFun() ######

#' WeightsFun
#'
#' Support function, generates all possible weights for NBI subproblems
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weights All possible weights for NBI subproblem
#' @keywords internal

WeightsFun = function(n, k){

  # package env variables
  # weight, Weights, Formers, Layer, lastone, currentone
  #
  # Generates all possible weights for NBI subproblems given:
  # n, the number of objectives
  # 1/k, the uniform spacing between two w_i (k integral)
  # This is essentially all the possible integral partitions
  # of integer k into n parts.

  assign("WeightSub", matrix(0,1,n), rMOSTenv)
  # WeightSub <<- matrix(0,1,n)
  assign("Weights", vector(), rMOSTenv)
  # Weights <<- vector()
  assign("Formers", vector(), rMOSTenv)
  # Formers <<- vector()
  assign("Layer", n, rMOSTenv)
  # Layer <<- n
  assign("lastone", -1, rMOSTenv)
  # lastone <<- -1
  assign("currentone", -1, rMOSTenv)
  # currentone <<- -1

  Weight_Generate(1, k)

  return(list(Weights = rMOSTenv$Weights, Formers = rMOSTenv$Formers))

}

###### Weight_Generate() ######

#' Weight_Generate
#'
#' Function intended to test the weight generation scheme for NBI for > 2 objectives
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weight_Generate
#' @keywords internal

Weight_Generate = function(n, k){

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
  #     num = dimFun(Weights)[2]
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
  #       num = dimFun(Weights)[2]
  #       currentone <<- num+1
  #     }
  #
  #     WeightSub[(Layer - n + 1)] <<- i
  #     Weight_Generate(n+1, k-i)
  #   }
  #
  # }

  if(n == rMOSTenv$Layer){

    if(rMOSTenv$currentone >= 0){
      rMOSTenv$Formers <- c(rMOSTenv$Formers, rMOSTenv$lastone)
      rMOSTenv$lastone <- rMOSTenv$currentone
      rMOSTenv$currentone <- -1
    }else{
      num = dimFun(rMOSTenv$Weights)[2]
      rMOSTenv$Formers <- c(rMOSTenv$Formers,num)
    }

    rMOSTenv$WeightSub[(rMOSTenv$Layer - n + 1)] <- k
    rMOSTenv$Weights <- cbind(rMOSTenv$Weights,t(rMOSTenv$WeightSub))

  }else{

    for(i in 0:k){
      if(n == (rMOSTenv$Layer - 2)){
        num = dimFun(rMOSTenv$Weights)[2]
        rMOSTenv$currentone <- num+1
      }

      rMOSTenv$WeightSub[(rMOSTenv$Layer - n + 1)] <- i
      Weight_Generate(n+1, k-i)
    }

  }

}

###### myLinCom_2C() ######

#' myLinCom_2C
#'
#' Support function
#' @param x Input predictor weight vector
#' @return f Criterion vector
#' @keywords internal

myLinCom_2C = function(x){

  # package env variable: g_Weight
  F   = myFM_2C(x)
  f = t(rMOSTenv$g_Weight)%*%F
  return(f)

}

###### myT() ######

#' myT
#'
#' Support function, define criterion space for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return f Temporary criterion space
#' @keywords internal

myT = function(x_t){

  f = x_t[length(x_t)]
  return(f)

}

###### myTCon_eq_2C() ######

#' myTCon_eq_2C
#'
#' Support function, define constraint condition for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return ceq Temporary constraint condition
#' @keywords internal

myTCon_eq_2C = function(x_t){

  # package env variables:
  # g_Normal g_StartF

  t = x_t[length(x_t)]
  x = x_t[1:(length(x_t)-1)]

  fe  = -myFM_2C(x) - rMOSTenv$g_StartF - t * rMOSTenv$g_Normal

  ceq1 = myCon_eq_2C(x)
  ceq = c(ceq1,fe)

  return(ceq)

}

###### plotPareto_2C() ######

#' plotPareto_2C
#'
#' Function for plotting Pareto-optimal curve and predictor weights
#' @param CriterionOutput Pareto-Optimal criterion solution
#' @param ParetoWeights Pareto-Optimal predictor weight solution
#' @importFrom grDevices rainbow
#' @importFrom graphics legend lines par points
#' @return Plot of Pareto-optimal curve and plot of predictor weights
#' @keywords internal

plotPareto_2C = function(CriterionOutput, ParetoWeights){

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow=c(1,2))

  Y1 = t(CriterionOutput[1,])
  Y2 = t(CriterionOutput[2,])
  X = t(ParetoWeights)

  # AI ratio - Composite Validity trade-off

  plot(Y1, Y2,
       xlim = c(min(Y1),max(Y1)),
       main = "Pareto Trade-Off Curve",
       xlab = "Rcy1",
       ylab = "Rcy2",
       type='c',col='blue')

  points(Y1, Y2,
         pch=8,col='red')

  # Predictor weights

  plot(Y1,X[,1],
       xlim=c(min(Y1),max(Y1)),ylim=c(0,1),
       main = "Predictor Weights",
       xlab = "Rcy1",
       ylab = "Predictor weight",
       type='c',col='red')
  points(Y1,X[,1],pch=8,
         col=rainbow(1))

  for(i in 2:ncol(X)){

    lines(Y1,X[,i],type='c',
          col=rainbow(1, start=((1/ncol(X))*(i-1)), alpha=1))
    points(Y1,X[,i],pch=8,
           col=rainbow(1, start=((1/ncol(X))*(i-1)), alpha=1))

  }

  legend('topleft',
         legend=c(paste0('Predictor ',1:ncol(X))),
         lty=c(rep(2,ncol(X))),lwd=c(rep(2,ncol(X))),
         col=rainbow(ncol(X)))

}
