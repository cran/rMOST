% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MOST.R
\name{MOST}
\alias{MOST}
\title{MOST}
\usage{
MOST(optProb, Rx, Rxy1, Rxy2, Rxy3, sr, prop1, prop2, d1, d2, Spac = 10)
}
\arguments{
\item{optProb}{Optimization problem. "3C" = no adverse impact objectives and three non-adverse impact objectives; "2C_1AI" = one adverse impact objective and two non-adverse impact objectives; "1C_2AI" = two adverse impact objectives and one non-adverse impact objective.}

\item{Rx}{Predictor intercorrelation matrix}

\item{Rxy1}{Needs to specify for all three types of optimization problems (optProb). Predictor criterion-related validity for non-adverse impact objective 1 (i.e., correlation between each predictor and non-adverse impact objective 1)}

\item{Rxy2}{Only specify if optimization problem is "3C" or "2C_1AI". Predictor criterion-related validity for non-adverse impact objective 2 (i.e., correlation between each predictor and non-adverse impact objective 2)}

\item{Rxy3}{Only specify if optimization problem is "3C". Predictor criterion-related validity for non-adverse impact objective 3 (i.e., correlation between each predictor and non-adverse impact objective 3)}

\item{sr}{Only specify if optimization problem is "2C_1AI" or "1C_2AI". Overall selection ratio.}

\item{prop1}{Only specify if optimization problem is "2C_1AI" or "1C_2AI". Proportion of minority1 in the applicant pool; prop1 = (# of minority1 applicants)/(total # of applicants)}

\item{prop2}{Only specify if optimization problem is "1C_2AI". Proportion of minority2 in the applicant pool; prop2 = (# of minority2 applicants)/(total # of applicants)}

\item{d1}{Only specify if optimization problem is "2C_1AI" or "1C_2AI". Vector of standardized group-mean differences between majority and minority 1 for each predictor; d1 = avg_majority - avg_minority1}

\item{d2}{Only specify if optimization problem is "1C_2AI". Vector of standardized group-mean differences between majority and minority 2 for each predictor; d2 = avg_majority - avg_minority2}

\item{Spac}{Determines the number of solutions.}
}
\value{
Pareto-Optimal solutions with objective values (e.g., C1, AI1) and the corresponding predictor weights (e.g., P1, P2)
}
\description{
Optimizes 3 objectives with normal boundary intersection algorithm
}
\details{
# Inputs required by optimization problems
Different types of optimization problems require different input parameters:
* optProb = "3C": MOST(optProb, Rx, Rxy1, Rxy2, Rxy3)
* optProb = "2C_1AI": MOST(optProb, Rx, Rxy1, Rxy2, sr, prop1, d1)
* optProb = "1C_2AI": MOST(optProb, Rx, Rxy1, sr, prop1, d1, prop2, d2)

# Notes regarding the inputs
* For personnel selection applications, all predictor-intercorrelations and criterion-related validity inputs should be corrected for range restriction and criterion unreliability to reflect the relations in the applicant sample.
* For optimization problems with 2 adverse impact objectives (i.e., optProb = "1C_2AI"), d1 and d2 should be the standardized mean difference between a minority group and the same reference group (e.g., Black-White and Hispanic-White, not Black-White and female-male)

# Optimization
* Optimization may take several minutes to run.
* Optimization may fail in some applications due to non-convergence.

For more details, please consult the vignette.
}
\examples{
# A sample optimization problem with 3 non-adverse impact objectives and 3 predictors
# For more examples, please consult the vignette.

# Specify inputs
# Predictor inter-correlation matrix (Rx)
Rx <- matrix(c(1,  .50, .50,
               .50,  1, .50,
               .50, .50,  1), 3, 3)

# Predictor-objective relation vectors (Rxy1, Rxy2, Rxy3)
# Criterion-related validities
## Criterion 1
Rxy1 <- c(-.30, 0, .30)
## Criterion 2
Rxy2 <- c(0, .30, -.30)
## Criterion 3
Rxy3 <- c(.30, -.30, 0)

# Get Pareto-optimal solutions
\donttest{
out <- MOST(optProb = "3C", Rx = Rx, Rxy1 = Rxy1, Rxy2 = Rxy2, Rxy3 = Rxy3, Spac = 10)
out
}
}
