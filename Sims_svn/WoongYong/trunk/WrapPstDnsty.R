#WrapPstDnsty <- function(param,eqs,z,trf,dfactor,prior=TRUE){
WrapPstDnsty <- function(param,...,prior.param){
	pstout <- try(PstDnsty(param,...,prior.param))
	if (inherits(pstout,"try-error")) {
		fv <- Inf
	} else {
		if (is.null(prior.param)) {fv <- pstout$lkh} else {fv <- (pstout$lkh + sum(pstout$prd))}
	}
  return(fv)
}