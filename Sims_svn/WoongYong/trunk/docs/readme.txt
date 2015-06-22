There are two main files:

FTmodel.runMaxPstDnsty.R     -> maximize posterior density
FTmodel.runPstAnalysis.R     -> generate some plots and marginal posterior density plots at a mode

***** directory structure

Below a working directory, the main files require directories to save results as follows:

working_directory/result/1
                         2
                         3

where 1 is for full sample, 2 is for the first half and 3 is for the second half.
I stored a current estimate for the full sample in ./FTmodel/result/1.

***** results

Results from posterior density maximization (FTmodel.runMaxPstDnsty.R) will be saved as 

FTmodel.(subsample).xhat.(try no).RData - estimates and optimization outputs
FTmodel.(subsample).xhat.(try no).log   - logs

For example,

FTmodel.1.xhat.102.RData
FTmodel.1.xhat.102.log

The log files will be saved only if R is told so (output_=1). 

Plots by FTmodel.runPstAnalysis.R will be saved as

FTmodel.(subsample).xhat.(try no).forecastanderrors.pdf - one-step ahead forecasts and errors
FTmodel.(subsample).xhat.(try no).irfs.pdf              - impulse response functions
FTmodel.(subsample).xhat.(try no).residuals.pdf         - residual plots (line, scatter)
FTmodel.(subsample).xhat.(try no).sshat.pdf             - one-step ahead predicted state vector

and estimates and residual comparisons with BVAR is saved as

FTmodel.(subsample).xhat.(try no).PstAnalysis.pdf

***** running main files

1) runMaxPstDnsty

The followings need to be set:

machine_ 	# 1: Laptop on C, 2: Laptop on H, 3: Server
data_ 		# 1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1
tryno_    # trial number (specified by user)
output_   # 1: output to screen, 2: output to file (sink)

div_      # growth rate restriction for gensysct (default <= 0.001)

2) runPstAnalysis

machine_ 	# 1: Laptop on C, 2: Laptop on H, 3: Server
data_ 		# 1: full sample 1954Q4~2005Q1, 2: 1954Q4~1979Q2, 3: 1983Q1~2005Q1
tryno_    # trial number (specified by user)
dev_ 			# graphic device 1: X11(), 2: pdf()
computeCuts_ 	# 1: compute marginal posterior density cuts 0: don't

div_      # growth rate restrictions for gensysct

***** to compute the posterior density of an estimate

Run FTmodel.runTestEnv.R and then, for example,

load("./result/1/FTmodel.1.xhat.102.RData")
estparam <- TrfParamBack(pstout.optim$par)
estparam["thet"] <- estparam["rho"] - (1-estparam["sig"])*estparam["gamm"]
estparam["rbar"] <- estparam["abar"] <- estparam["gamm"] + estparam["thet"] + estparam["pibar"]
pstd <- PstDnsty(param=estparam,SysEq=SysEq,z=y,bp0=bp0,trf=FALSE,dfactor=1,div=div_,prior.param=prior.param)

Two lines definining theta, rbar and abar are not needed here. But they are necessary if you want
to use gensysct() or SolveGensys().

***** files

* main files
FTmodel.runMaxPstDnsty.R     -> maximize posterior density
FTmodel.runPstAnalysis.R     -> generate some plots and marginal posterior density plots at a mode

* data and model
datawork.R                   -> compile data
FTmodel.data.RData           -> data (y1: without bond, y2: with bond)
FTmodel.SysEq.RData          -> equations, variables and parameters

* files to keep the main files from getting messy
FTmodel.runLibraries.R       -> load sources
FTmodel.prior.10n.R          -> prior definition for subsamples (n=1,2,3)

* functions
PstDnsty6.R       -> compute posterior densities (function - PstDnsty())
PstDnsty7.R       -> same as PstDnsty6.R, but returns estimated state vectors and et cetra
WrapPstDnsty.R    -> a wrapper for PstDnsty()
SolveGensys.R     -> call by PstDnsty() or other  routines to solve the model, also returns steady state values
Pnlty.R           -> penalty function, called by SolveGensys()
PriorDnsty.R      -> compute prior densities (need prior definition)

TrfParam.R        -> reparameterization for bounded optimization
TrfParamBack.R    -> back to original parameterization

numHess2.R        -> direct computation of numerical Hessian
numgrad2.R        -> two-sided numerical derivative

numIntCnst.R      -> numerically compute integration of matrix exponential for constant
numIntCov.R       -> numerically compute integration of covariance matrix (dtInnovCov)

logitbeta.R       -> logit-transformed-like beta distribution (for phi0)

* miscellaneous

FTmodel.runSysEq.R           -> generate the system of equations
FTmodel.runTestEnv.R         -> set working directory, load sources and load data to use PstDnsty() or SolveGensys()
FTmodel.runDescriptiveStat.R -> calculate some descriptive statistics