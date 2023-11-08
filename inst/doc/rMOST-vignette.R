## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(rMOST)

## -----------------------------------------------------------------------------
## Input ##

# Predictor intercorrelation matrix
Rx <- matrix(c(  1,  .37, .51, .16, .25,
               .37,   1, .03, .31, .02,
               .51, .03,   1, .13, .34,
               .16, .31, .13,   1,-.02,
               .25, .02, .34,-.02,   1), 5, 5)

# Criterion validity of the predictors
Rxy1 <- c(.32, .52, .22, .48, .20) 
Rxy2 <- c(.30, .35, .15, .25, .10)
Rxy3 <- c(.15, .25, .30, .35, .10)

# Overall selection ratio
sr <- 0.15

# Proportion of minority applicants
prop_b <- 1/8 # Proportion of Black applicants (i.e., (# of Black applicants)/(# of all applicants))
prop_h <- 1/6 # Proportion of Hispanic applicants

# Predictor subgroup d
d_wb <- c(.39, .72, -.09, .39, .04) # White-Black subgroup difference
d_wh <- c(.17, .79, .08, .04, -.14) # White-Hispanic subgroup difference

## -----------------------------------------------------------------------------
# Example: 3 non-adverse impact objectives
out_3C = MOST(optProb = "3C", 
              # predictor intercorrelations
              Rx = Rx, 
              # predictor - objective relations
              Rxy1 = Rxy1, # non-AI objective 1
              Rxy2 = Rxy2, # non-AI objective 2
              Rxy3 = Rxy3, # non-AI objective 3
              Spac = 10)

# The first few solutions
head(out_3C)

## -----------------------------------------------------------------------------
# Example: 2 non-adverse impact objectives & 1 adverse impact objective 
out_2C_1AI = MOST(optProb = "2C_1AI", 
                  # predictor intercorrelations
                  Rx = Rx, 
                  # predictor - objective relations
                  Rxy1 = Rxy1, # non-AI objective 1
                  Rxy2 = Rxy2, # non-AI objective 2
                  d1 = d_wb, # subgroup difference for minority 1
                  # selection ratio
                  sr = sr, 
                  # proportion of minority
                  prop1 = prop_b, # minority 1
                  Spac = 10)

# The first few solutions
head(out_2C_1AI)

## -----------------------------------------------------------------------------
# Example: 1 non-adverse impact objective & 2 adverse impact objectives 
out_1C_2AI = MOST(optProb = "1C_2AI",
                  # predictor intercorrelations
                  Rx = Rx, 
                  # predictor - objective relations
                  Rxy1 = Rxy1, # non-AI objective 1
                  d1 = d_wb, # subgroup difference for minority 1
                  d2 = d_wh, # subgroup difference for minority 2
                  # selection ratio
                  sr = sr, 
                  # proportion of minority 
                  prop1 = prop_b, # minority 1
                  prop2 = prop_h, # minority 2
                  Spac = 10)

# The first few solutions
head(out_1C_2AI)

## -----------------------------------------------------------------------------
out_3C[1, 6:ncol(out_3C)]

## -----------------------------------------------------------------------------
out_3C[1, 3:5]

