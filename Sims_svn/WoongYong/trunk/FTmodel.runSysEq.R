# lam <- log(lam)
# variables : c, cc(cdot), b, tau, r, a, p, w(pdot), lam
# shocks    : epsm, epsr, epspc, epstau
# parameters: phi0 0.5, phi1 0.4, phi2 0.75, gamm 0.02, rho 0.01, sig 2, bet 0.1, delt 0.2, omega 1, psi 2, pibar 0.03, cbar 0.01
# thet = rho - (1-sig)*gamm
# rbar = rho + sig * gamm + pibar

### AR shocks

eqchab <- expression(mpolicy=rdot - (-phi0 * ((r - rbar) - phi1 * (pdot - pibar) - phi2 * cdot) + epsm),
IS=-lamdot - (a - adot/a - pdot - gamm - thet + xir),
termstruc=r - (a - adot/a),
phcurve=wdot - (bet * ((pdot - pibar) - delt * (c - cbar)) - xipc),
GBC=bdot - ((a - adot/a - pdot - gamm)*b - tau),
fpolicy=taudot - (omega*exp(c)*cdot + epstau),
lamdef=exp(lam) - exp(-sig * c-.5*psi*(1 - sig)*cdot^2)*(1 + psi*((1-sig)*cdot^2 - (thet + psi*(1-sig)*cdot*ccdot)*cdot + ccdot)),
#lamdef=lam - (-sig * c -.5*psi*(1 - sig)*cdot^2 +log(1 + psi*((1-sig)*cdot^2 - (thet + psi*(1-sig)*cdot*ccdot)*cdot + ccdot))),
ISshock=xirdot - (-rhor*xir + epsr),
PCshock=xipcdot - (-rhopc*xipc + epspc),
wdef=pdot - w,
ccdef=cc - cdot)

class(eqchab) <- c("eqsys","expression")

vars <- c("r", "a", "p", "w", "c", "cc", "b", "tau", "lam","xir","xipc")
vards <- paste(vars,"dot",sep="")
shocks <- c("epsm","epstau","epsr","epspc")
expeqs <- c(2,3,4,7)
params <- c("phi0","phi1","phi2","gamm","rho","sig","bet","delt","omega","psi","pibar","cbar","taubar","rhor","rhopc","sig2m","sig2tau","sig2r","sig2pc","rhob","sig2b")

SysEq <- list(eqs=eqchab,vars=vars,vards=vards,shocks=shocks,expeqs=expeqs,params=params)