.onLoad <- function(libname, pkgname) {
  ## separator between type (sd/var/cor/cov) and following elements in
  ## ran_pars terms
  options(broom.mixed.sep1 = "__")
  ## separator between elements in cor/cov terms
  options(broom.mixed.sep2 = ".")
}
