library("R.matlab")

### Makes a list of all relevant parameters that characterize the
### model we are iteration on, and also some environment variables
### relevant for calling dynare.
makedynframework <- function(tempdir="/tmp",     # ramdisk recommended
                             dynpath,            # path of dynare
                             var,                # vector of variables (char)
                             varexo,             # exogenous variables (char)
                             parameters,         # parameter names (char)
                             model,              # list of strings or filename
                             initval,            # make up something
                             order=2,            # order of approximation
                             dynparams="--sim 0 --per 0",
                                                 # parameters for dynare
                             statevars           # state variables
                             ) {
  ## checks
  stopifnot(order >= 1)
  ## read model if given as a filename, otherwise append "model;",
  ## "end;", and semicolon at the end of each line
  if (length(model)==1)                 # this is a filename
    model <- readLines(model)
  else {                                # not a filename
    model <- sapply(model,function(a) sprintf("%s;",a)) # semicolons
    model <- c("model;",model,"end;")
  }
  ## construct framework
  dynfwork <- list(tempdir=tempdir,
                   dynpath=dynpath,
                   var=var,
                   varexo=varexo,
                   parameters=parameters,
                   model=model,
                   initval=initval,
                   order=as.integer(order),
                   dynparams=dynparams,
                   statevars=statevars
                   )
}

### Output lists of the form "keyword a, b, c, d;"
catenum <- function(con,keyword,var,append=TRUE) {
    cat(file=con,keyword,sapply(var[-length(var)], # a, b, c, ...
          function(a) sprintf("%s,",a)),
        sprintf("%s;\n",var[length(var)]),append=append) # d; (last variable)
}

### Output lists of the form "a=1; b=2, c=3;", with linebreaks
catvalues <- function(con, varnames, values) {
  stopifnot(length(varnames)==length(values))
  writeLines(con=con,mapply(function(a,b)
        sprintf("%s=%f;",a,b),varnames,values))
}

## This function generates a model file with the given parameters,
## calls dynare, analyzes the journal (eg for indeterminacy), reads
## the matrices in Matlab format, and returns a list.
calldynare <- function(dynfw,params,vcovmat) {
  with(dynfw, {
    ## test for consistency
    stopifnot(is.matrix(vcovmat) && (nrow(vcovmat)==ncol(vcovmat))
              && is.finite(vcovmat))
    stopifnot(is.vector(params) && is.finite(params))
    ## construct model file
    con <- file(sprintf("%s/tmp.mod",tempdir),"w") # open connection
    catenum(con,"var",var)
    catenum(con,"varexo",varexo)
    catenum(con,"parameters",parameters)
    catvalues(con,parameters,params)    # parameters
    cat(file=con,"\n")
    writeLines(model,con)               # model
    cat(file=con,"\n\ninitval;\n")
    catvalues(con,var,initval)          # initval
    cat(file=con,"end;\n\nvcov = [\n")
    for (i in 1:nrow(vcovmat)) {        # variance matrix
      cat("  ", vcovmat[i,], file=con)
      if (i != nrow(vcovmat))           # semicolon for all but last line
        cat(";", file=con)
    }
    cat(file=con,"\n];\n\n")
    cat(file=con,sprintf("order=%d;\n",order)) # order
    close(con)                          # close connection
    ## remove journal and mat files now (if any remained from previous
    ## calls), so if dynare fails, reading them will generate and error
    unlink(sprintf("%s/tmp.jnl",tempdir))
    unlink(sprintf("%s/tmp.mat",tempdir))
    ## call dynare
    dynout <- system(sprintf("%s %s %s/tmp.mod",dynpath,dynparams,
                             tempdir),TRUE,FALSE)
    ## analyze the journal file
    journal <- readLines(sprintf("%s/tmp.jnl",tempdir))
    ## the following is just an example
    isunstable <- length(grep("model not stable", journal)) > 0
    ## read the results, return as a list
    if (isunstable)
      dyn <- NULL
    else
      dyn <- readMat(sprintf("%s/tmp.mat",tempdir))
    list(dyn=dyn,isunstable=isunstable,dynout=dynout)
  })
}

## We need the mapping from the folded tensor to the unfolded one.
## The folding algorithm is explained in Section 230 (page 91) of the
## documentation file tl.pdf, url:
## http://www.cepremap.cnrs.fr/juillard/mambo/download/manual/dynare++/tl.pdf
## The function is actually a wrapper, to account for the fact that R
## starts indices from 1, not 0.  ii is a vector if indices (starting
## from ONE!), n is the total number of variables.  We also need to
## sort the variables.

getoffset <- function(ii,n) {
  getoffsetinternal(sort(ii)-1,n)+1
}

getoffsetinternal <- function(ii,n) {   # this is for indexing from 0
  k <- length(ii)
  if (k==0)
    0
  else
    choose(n+k-1,k)-choose(n-ii[1]+k-1,k)+
      getoffsetinternal(ii[-1]-ii[1],n-ii[1])
}

## ## I wrote some basic functions to test the code above, they are not
## ## needed for anything else.
## ## This section is commented out, as I only used it for testing.
## genoffset <- function(n,k,i) {
##   if (k==1)
##     matrix(i:n,n-i+1,1)
##   else
##     do.call(rbind,lapply(i:n,function (j) cbind(j,genoffset(n,k-1,j))))
## }
## testoffset <- function(n,k) {
##   oo <- genoffset(n,k,1)
##   trueo <- as.vector(apply(oo,1,function(o)
##                            getoffset(o[sample(1:k)],n)))
##   all.equal(trueo,1:choose(n+k-1,k))
## }
## testrandomoffsets <- function() {
##   for (i in 1:1000) {
##     if (isTRUE(testoffset(sample(6,1),sample(6,1))))
##       cat(".")                          # good
##     else
##       cat("!")                          # bad
##   }
##   cat("\n")
## }

## Converts a flat vector index into an array index.  k is the
## dimension, n is the number of elements along each dimension.  As
## opposed to a base conversion, numbers are reversed for an array.
## Also, we add and subtract 1 (R indexing).
mapindex <- function(m,n,k) {
  mm <- vector(mode="numeric",length=k)
  m <- m-1
  for (i in seq(along=mm)) {
    mm[i] <- m %% n
    m <- m %/% n
  }
  mm+1
}

## unfolds a single tensor (a row of the folded matrix)
unfoldtensor <- function(folded,        # the folded tensor (from dynare)
                         k,             # the order
                         n) {           # total number of variables
  a <- vector(mode="numeric",length=n^k)
  for (i in seq(along=a))
    a[i] <- folded[getoffset(mapindex(i,n,k),n)]
  array(a,dim=rep(n,k))
}

## match variables names
matchvarnames <- function(varnames,dynnames) {
  ## check
  stopifnot(length(varnames)==length(dynnames))
  ## strip padding spaces
  dynnames <- sapply(dynnames,function(a) gsub(" ","",a))
  ## find matches
  pp <- vector(mode="numeric",length=length(varnames))
  for (i in seq(along=pp)) {
    m <- grep(varnames[i],dynnames)     # find matches
    if (length(m)!=1)
      stop(sprintf("There are %d matches for variable %s in matchvarnames!\n",
                   length(m),varnames[i]))
    pp[i] <- m[1]                       # store match
  }
  pp
}
