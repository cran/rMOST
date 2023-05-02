# By convention, .onLoad() and friends are usually saved in a file called R/zzz.R

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to rMOST.\nrMOST estimates Pareto-optimal solutions for personnel selection with 3 objectives using Normal Boundary Intersection algorithm.")
}


.onLoad <- function(libname, pkgname) {

  # Variables that need to be accessed by multiple functions simultaneously need to be stored in an environment accessible by all functions (e.g., package environment)

  # Do not use Global environment to respect R landscape of users
  assign("rMOSTenv", new.env(parent = emptyenv()), parent.env(environment()))
  assign("rMOSTenv_3C", new.env(parent = emptyenv()), parent.env(environment()))
  assign("rMOSTenv_2C", new.env(parent = emptyenv()), parent.env(environment()))
  assign("rMOSTenv_1C_1AIR", new.env(parent = emptyenv()), parent.env(environment()))
  assign("rMOSTenv_1C_2AIR", new.env(parent = emptyenv()), parent.env(environment()))
  assign("rMOSTenv_2C_1AIR", new.env(parent = emptyenv()), parent.env(environment()))
  # rMOSTenv <- new.env(parent = emptyenv())
  # Variables and functions can be retrieved by rMOSTenv_**$[VARNAME]

}
