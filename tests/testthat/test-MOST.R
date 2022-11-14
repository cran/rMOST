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

exsol_3C <-
  structure(
    list(
      Solution = c(1, 2, 3, 4, 5, 6),
      Optimized = c(
        "C1, C2, C3",
        "C1, C2, C3",
        "C1, C2, C3",
        "C1, C2, C3",
        "C1, C2, C3",
        "C1, C2, C3"
      ),
      C1 = c(0.658, 0.658, 0.657, 0.656, 0.654, 0.653),
      C2 = c(0.4,
             0.405, 0.408, 0.41, 0.413, 0.414),
      C3 = c(0.424, 0.423, 0.421,
             0.419, 0.416, 0.413),
      P1 = c(0.029, 0.059, 0.084, 0.106, 0.126,
             0.144),
      P2 = c(0.399, 0.406, 0.413, 0.422, 0.43, 0.44),
      P3 = c(0.094,
             0.085, 0.078, 0.072, 0.068, 0.065),
      P4 = c(0.346, 0.353, 0.356,
             0.358, 0.357, 0.355),
      P5 = c(0.161, 0.157, 0.153, 0.149, 0.145,
             0.141)
    ),
    row.names = c(NA, 6L),
    class = "data.frame"
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

exsol_2C_1AI <-
  structure(
    list(
      Solution = c(1, 2, 3, 4, 5, 6),
      Optimized = c(
        "C1, C2, AI1",
        "C1, C2, AI1",
        "C1, C2, AI1",
        "C1, C2, AI1",
        "C1, C2, AI1",
        "C1, C2, AI1"
      ),
      C1 = c(0.22, 0.22, 0.22, 0.22, 0.22, 0.638),
      C2 = c(0.15,
             0.15, 0.15, 0.15, 0.15, 0.417),
      AI1 = c(1.147, 1.147, 1.147,
              1.147, 1.147, 0.269),
      P1 = c(0, 0, 0, 0, 0, 0.176),
      P2 = c(0,
             0, 0, 0, 0, 0.456),
      P3 = c(1, 1, 1, 1, 1, 0),
      P4 = c(0, 0, 0,
             0, 0, 0.27),
      P5 = c(0, 0, 0, 0, 0, 0.098)
    ),
    row.names = c(NA,
                  6L),
    class = "data.frame"
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

exsol_1C_2AI <-
  structure(
    list(
      Solution = c(1, 2, 3, 4, 5, 6),
      Optimized = c(
        "C1, AI1, AI2",
        "C1, AI1, AI2",
        "C1, AI1, AI2",
        "C1, AI1, AI2",
        "C1, AI1, AI2",
        "C1, AI1, AI2"
      ),
      C1 = c(0.652, 0.52, 0.489, 0.458, 0.426, 0.395),
      AI1 = c(0.355, 0.629, 0.684, 0.738, 0.794, 0.849),
      AI2 = c(0.483,
              0.687, 0.714, 0.74, 0.768, 0.794),
      P1 = c(0.025, 0, 0, 0, 0,
             0),
      P2 = c(0.328, 0.151, 0.131, 0.112, 0.092, 0.072),
      P3 = c(0.164,
             0.45, 0.501, 0.551, 0.601, 0.651),
      P4 = c(0.333, 0.28, 0.262,
             0.244, 0.226, 0.208),
      P5 = c(0.151, 0.119, 0.106, 0.093, 0.081,
             0.068)
    ),
    row.names = c(NA, 6L),
    class = "data.frame"
  )

# Inputs for non-convergence scenario

error_Rx <- matrix(.1, nrow = 3, ncol = 3, dimnames = list(paste0("P", 1:3), paste0("P", 1:3)))
diag(error_Rx) <- 1

error_Rxy <- matrix(c(-.2, 0, .2,
                      0, .2, -.2,
                      .2, -.2, 0),
                    nrow = 3, ncol = 3, byrow = T,
                    dimnames = list(paste0("C", 1:3), paste0("P", 1:3)))

test_that("Correct number of solutions in the example 3C problem", {
  expect_equal(dim(out_3C), c(129, 10))
})

test_that("Correct number of solutions in the example 2C_1AI problem", {
  expect_equal(dim(out_2C_1AI), c(129, 10))
})

test_that("Correct number of solutions in the example 1C_2AI problem", {
  expect_equal(dim(out_1C_2AI), c(108, 10))
})

test_that("Should trigger an equality constraint error", {
  expect_error(MOST(
    optProb = "3C",
    Rx = error_Rx,
    Rxy1 = error_Rxy["C1", ],
    Rxy2 = error_Rxy["C2", ],
    Rxy3 = error_Rxy["C3", ],
    Spac = 10
  ), "equality constraints in x0 returns NA")
})





