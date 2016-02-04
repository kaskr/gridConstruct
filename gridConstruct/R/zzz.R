.First.lib <- function(lib, pkg) {
  library.dynam("gridConstruct", pkg, lib)
}

##.onAttach <- function(libname, pkgname)
##.onLoad <- function(libname, pkgname)
