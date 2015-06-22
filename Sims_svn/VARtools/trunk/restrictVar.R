#' restrict VAR
#'
#' Posterior odds on linear restrictions on a VAR.
#'
#' The calculations are done as though the prior on the restricted space is
#' proporitonal to the prior density on the unrestricted space, confined to
#' the restricted space.  This makes sense if the restriction is approximating
#' a restriction to a small neighborhood of the restricted space.
#' For the case of \code{SVARhtskdmdd} output as input to this function, the
#' results are conditional on the modal \code{A0} and \code{lmd}, which must
#' be included in the \code{fitdata} argument. For the
#' case of \code{mgnldensity()} output, the results condition on the modal
#' residual covariance matrix, as if it were known.  
#' Restrictions can be specified as rows of rmat, with coefficients applied to
#' elements of \code{By} and \code{Bx} stacked as they are in \code{xxi} (and
#' then repeated across the equation index), or
#' they can be specified in \code{yzrone}, \code{xzrone}.  Each zero element
#' of \code{yzrone} or \code{xzrone} generates a restriction that sets the
#' corresponding coefficient in \code{By} or \code{Bx} to zero (or to a
#' constant, if \code{!is.null(const)}).  Both kinds of restrictions can be
#' non-trivial in the same call.
#'
#' @param fitdata Output from a var estimate with a proper prior from
#'                \code{mgnldnsty()} or \code{SVARhtskmdd()})
#' @param type Either \code{mdd} or \code{SVmdd}, depending on which
#'             program produced \code{fitdata}.
#'
#' @param const    The right hand side of \code{rmat %*% coeff = const}, not
#'                 the constant in the var.
#' 
#' @param cyzr,cxzr  If using \code{yzrone}, \code{xzrone} with non-trivial
#'                   constants, leave \code{const=NULL} and specify constants
#'                   with \code{cyzr} and \code{cxzr}
#'
#' @param sc The posterior odds ratio.  It will be similar, at least in sign,
#'           to the Schwarz criterion on small models, but if large numbers
#'           of coefficients are being tested at once, the prior, which the
#'           Schwarz criterion ignores, may be very important.
#'------------------------------------------------
#' @return  A list with elements
#'   \describe{
#'    \item{chiSquared}{The usual chi-squared statistic}
#'    \item{df}{Degrees of freedom of the chi-squared statistic (number of
#'             restrictions)}
#'    \item{sc}{Posterior odds ratio against the restriction.)
#'    \item{pval}{p-value of the chi-squared}
#'    \item{ singxxi,singsig}{If either of these are TRUE, there might be
#'     numerical error}
#'   }
#' @import tensor
#' @export
#' 
restrictVAR <- function(fitdata, type=c("mdd", "SVmdd"), rmat=NULL,
                        yzrone=NULL, xzrone=NULL,const=NULL, cyzr=NULL,
                        cxzr=NULL) {
    if (length(type) > 1) type <- type[1]
    if (type == "SVmdd") {
        ##  require("tensor")  #not needed if packaged, because of @import
        bvw <- fitdata
        vout2 <- list(post=bvw$vout$var, prior=bvw$vout$varp)
    } else {                            ##mdd
        vout2 <- list(post=fitdata$var, prior=fitdata$varp)
    }
    neq <- dim(vout2[[1]]$By)[1]
    ny <- dim(vout2[[1]]$By)[2]
    lags <- dim(vout2[[1]]$By)[3]
    nx <- dim(vout2[[1]]$Bx)[2]
    ncf <- ny * lags + nx
    if (is.null(rmat)) {
        rmat <- matrix(0, 0, ncf *neq)
    }
    if (!is.null(yzrone)) {
        byz <- which(yzrone == 0, arr.ind=TRUE)
        nrstr <- dim(byz)[1]
        if (is.null( cyzr)) cyzr <- array(0, dim(yzrone))
        for (ir in 1:nrstr ) {
            newrow <- rep(0, neq * ncf)
            newrow[(byz[ir,1] - 1) * ncf + (byz[ir, 3] -1) * ny +
                       byz[ir, 2]] <- 1
            rmat <- rbind(rmat,newrow)
        }
        const <- c(const, cyzr[byz])
    }
    if (!is.null(xzrone)) {
        bxz <- which(xzrone == 0, arr.ind=TRUE )
        nrstr <- dim(bxz)[1]
        if (is.null(cxzr)) cxzr <- matrix(0, neq, nx)
        for (ir in 1:nrstr)  {
            newrow <- rep(0,ncf * neq)
            newrow[(bxz[ir,1] - 1) * ncf + ny * lags + bxz[ir, 2]] <- 1
            rmat <- rbind(rmat, newrow)
        }
        const <- c(const, cxzr[bxz])
    }
    svdr <- svd(rmat)
    if (max(abs(svdr$d)) > 1e10 * min(abs(svdr$d))) {
        error("restrictions not full rank")
    }
    result <- list(post=NULL, prior=NULL)
    for (ivo in 1:2) {
        vout <- vout2[[ivo]]
        T <- dim(vout$u)[1]
        if (type == "mdd") {
            sig <- cov(vout$u)
            svdsig <- svd(sig)
            singsig <- (max(svdsig$d) > 1e10 * min(svdsig$d))
            if(singsig) warning("Near singular sig matrix in restrictVAR")
            svdxxi <- svd(vout$xxi)
            singxxi <- (max(svdxxi$d) > 1e10 * min(svdxxi$d))
            if(singxxi) warning("Near singular xxi matrix in restrictVAR")
            ## schwarz <- rmat %*% kronecker(svdsig$u %*% diag(1/sqrt(svdsig$d)), svdxxi$u %*% diag(1/sqrt(svdxxi$d)))
            ##schwarz <- kronecker((1/sqrt(svdsig$d)) * t(svdsig$u), (1/sqrt(svdxxi$d)) * t(svdxxi$u)) %*% rv  #2013.5.9
            ## sqrtVb <- kronecker(sqrt(svdsig$d) * t(svdsig$u), 1/sqrt(svdxxi$d)) * t(svdxxi$u)
            ## line above seems to be a mistake, since xxi is already x'x-inverse
            sqrtVb <- kronecker(sqrt(svdsig$d) * t(svdsig$u), sqrt(svdxxi$d) *
                                    t(svdxxi$u))
            dgVb <- apply(sqrtVb^2, 2, sum)
            rmatC <- rmat %*% diag(sqrt(T * dgVb))
            sqrtVbC <- sqrtVb %*% diag(1/sqrt(T * dgVb))
            lndetVb <- sum(log(svdsig$d)) * dim(vout$xxi)[1] + sum(log(svdxxi$d)) *
                dim(sig)[1]
            lndetVbC <- lndetVb - sum(log(dgVb * T))
            ## } else if (type == "KF") {          #type=="KF"
            ## svdVb <- svd(vout$Vb)
            ## sqrtVb <- sqrt(diag(svdVb$d)) %*% t(svdVb$u)
            ## dgVb <- diag(vout$Vb)
            ## rmatC <- rmat %*% diag(sqrt(T * dgVb))
            ## sqrtVbC <- sqrtVb %*% diag(1/sqrt(T * dgVb))
            ## lndetVb <- sum(log(svdVb$d))
            ## lndetVbC <- lndetVb - sum(log(dgVb * T))
            ## ## schwarz <- rmat %*% svdVb$u %*% diag(1/sqrt(svdVb$d)) #below is more efficient version for large Vb
            ## ## schwarz <- (1/sqrt(svdVb$d)) * (t(svdVb$u) %*% rv)
        } else {            # (type == "SVmdd")
            neq <- dim(vout$By)[1]
            ##-----------------------
            ## temp fix to use Karthik's output
            if(with(vout, exists("xx"))) {
                vout$xxi <- array(0, dim(vout$xx))
                for (iq in 1:neq) vout$xxi[ , , iq] <- solve(as.matrix(vout$xx[ , , iq]))
            }
            ##-------------------
            nX <- dim(vout$xxi)[1]
            Vb <- matrix(0, nX * neq, nX * neq)
            for (iq in 1:neq) {
                Vb[nX * (iq-1) + 1:nX, nX * (iq-1) + 1:nX] <- vout$xxi[ , , iq]
            }
            A0i <- solve(bvw$A)
            Vb <- kronecker(A0i, diag(nX)) %*% Vb %*% kronecker(t(A0i), diag(nX))
            svdVb <- svd(Vb)
            sqrtVb <- sqrt(diag(svdVb$d)) %*% t(svdVb$u)
            ## dgVb <- diag(Vb)
            ## rmatC <- rmat %*% diag(sqrt(T * dgVb))
            ## sqrtVbC <- sqrtVb %*% diag(1/sqrt(T *dgVb))
            ## lndetVb <- sum(log(svdVb$d))
            ## lndetVbC <- lndetVb - sum(log(dgVb * T))
        }
        svdvr <- svd(sqrtVb %*% t(rmat))
        ## svdvrC <- svd(sqrtVbC %*% t(rmatC)) #result == line above?
        ## vdim1 <- dim(svdvr$u)[1]
        ## svdvrp <- svd(diag(vdim1) - svdvr$u %*% t(svdvr$u), nu=vdim1 - dim(rmat)[1])
        ## svdvrpC <- svd(diag(vdim1) - svdvrC$u %*% t(svdvrC$u), nu=vdim1 -
        ##                   dim(rmat)[1])
        ## svdvrpuv <- svd(crossprod(svdvrp$u, t(sqrtVb))) 
        ## svdvrpuvC <- svd(crossprod(svdvrpC$u, t(sqrtVbC))) 
        ## lndetUR <- (sum(log(svdvrpuv$d))
        ## lndetURC <- sum(log(svdvrpuvC$d))
        df <- dim(rmat)[1]
        ## schwarz <- -2 * sum(log(diag(chol(crossprod(schwarz)))))   +  df * log(2 * pi)
        ## schwarz <- lndetVb - 2 * lndetUR + df * log(2 * pi)
        ## schwarzC <- lndetVbC - 2 * lndetURC + df * log(2 * pi)
        schwarX <- diag(1/sqrt(svdVb$d)) %*% t(svdVb$v) %*% svdr$v
        schwarX <- -2 * sum(log(abs(svd(schwarX, nu=0)$d)))
        if(is.null(const)) const <- rep(0, dim(rmat)[1])
        if(type == "SVmdd") {
            if (exists("tensor")) {
                vout$By <- tensor(A0i, vout$By, 2, 1)
            } else {                    #workaround
                vby <- A0i %*% matrix(vout$By, nrow=dim(vout$By)[1])
                vout$By <- array(vby, dim(vout$By))
            }
            vout$Bx <- A0i %*% vout$Bx
        }
        stackedcf <- c(t(cbind(matrix(vout$By, nrow=neq), vout$Bx)))
        gap <- rmat %*% stackedcf - const
        ##svdv <- svd(rmat %*% vout$Vb %*% t(rmat))
        chstat <- (1/svdvr$d) * t(svdvr$v) %*%  gap
        chstat <- crossprod(chstat)
        result[[ivo]] <-
            list(chiSquared=chstat, df=df,   #sc=schwarz, 
                          pval=pchisq(chstat,df),
                          ## sc2 = schwarz - (ncf*neq-df)*log(1 - df/(neq*ncf)),
                          ## scC=schwarzC,
                 schwarX=schwarX )
    }
  result[[3]] <- list(or = result[[1]]$chiSquared - result[[2]]$chiSquared + result[[1]]$schwarX -
                          result[[2]]$schwarX)
  return(result)       
}
