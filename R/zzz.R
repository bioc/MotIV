.onAttach <- function(libname, pkgname){
  data(jaspar2010)
  data(jaspar2010_scores)
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.13")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}
