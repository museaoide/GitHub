print.eqsys <- function(x) {
  for (iq in 1:length(x)) {
   cat(paste(names(x)[iq], ifelse(attr(x,"forward")[iq], "*", ""), sep=""), "\n", as.character(x[iq]),"\n\n")
  }
  cat("##------------\n\n",sep="")
  cat("vlist\n", attr(x, "vlist"),"\n\n")
  cat("param\n", attr(x,"param"),"\n\n")
  cat("shock\n", attr(x,"shock"), "\n\n")
  invisible(x)
}
