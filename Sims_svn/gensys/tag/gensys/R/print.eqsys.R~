print.eqsys <- function(x) {
  for (iq in 1:length(x)) {
   cat(paste(names(x)[iq], ifelse(attr(x,"forward")[iq], "*", ""), sep=""), "\n", as.character(x[iq]),"\n\n")
  }
  cat("##------------\n\n",sep="")
  vlist <- attr(x, "vlist")
  if (length(vlist) ==0) vlist <- "NONE"
  cat("vlist\n", vlist,"\n\n")
  param <- attr(x,"param")
  if (length(param) == 0) param <- "NONE"
  cat("param\n", param,"\n\n")
  shock <- attr(x,"shock")
  if (length(shock) == 0) shock <- "NONE"
  cat("shock\n", shock, "\n\n")
  invisible(x)
}
