## Input ##

# Bobko-Roth correlation matrix
Rx <- matrix(c(  1,  .37, .51, .16, .25,
                 .37,   1, .03, .31, .02,
                 .51, .03,   1, .13, .34,
                 .16, .31, .13,   1, -.02,
                 .25, .02, .34, -.02,   1), 5, 5)

Rxy1 <- c(.32, .52, .22, .48, .20) # Criterion validity of the predictors
Rxy2 <- c(.30, .35, .15, .25, .10)
Rxy3 <- c(.15, .25, .05, .45, .10)

# Overall selection ratio
sr <- 0.15

# Proportion of minority applicants
prop_b <- 1 / 8 # Proportion of Black applicants (i.e., (# of Black applicants)/(# of all applicants))
prop_h <- 1 / 6 # Proportion of Hispanic applicants

# Predictor subgroup d
d_wb <-
  c(.39, .72, -.09, .39, .04) # White-Black subgroup difference
d_wh <-
  c(.17, .79, .08, .04, -.14) # White-Hispanic subgroup difference

out_3C <- MOST(
  optProb = "3C",
  Rx = Rx,
  Rxy1 = Rxy1,
  Rxy2 = Rxy2,
  Rxy3 = Rxy3,
  Spac = 10
)

out_2C_1AI <-
  MOST(
    optProb = "2C_1AI",
    Rx = Rx,
    Rxy1 = Rxy1,
    Rxy2 = Rxy2,
    sr = sr,
    prop1 = prop_b,
    d1 = d_wb,
    Spac = 10
  )

out_1C_2AI <- MOST(
  optProb = "1C_2AI",
  Rx = Rx,
  Rxy1 = Rxy1,
  sr = sr,
  prop1 = prop_b,
  prop2 = prop_h,
  d1 = d_wb,
  d2 = d_wh,
  Spac = 10
)

test_that("Correct number of solutions in the example 3C problem", {
  expect_equal(dim(out_3C), c(129, 10))
})

test_that("Correct number of solutions in the example 2C_1AI problem", {
  expect_equal(dim(out_2C_1AI), c(129, 10))
})

test_that("Correct number of solutions in the example 1C_2AI problem", {
  expect_equal(dim(out_1C_2AI), c(108, 10))
})





