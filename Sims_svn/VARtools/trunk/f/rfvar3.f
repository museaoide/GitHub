subroutine rfvar3( By,Bx,u,xxi, ydata,lags,xdata,breaks,lambda,mu)
implicit none
real*8, intent=out, allocatable :: By, Bx, u, xxi
real*8, intent=in :: ydata(:,:)
real*8, intent=in, optional :: xdata(:,:), lambda, mu
integer, intent=in :: lags 
integer, intent=in, optional :: breaks(:)
logical ::  isx
integer :: nbreaks,n,smpl(0),Tsmpl
real*8 :: X, ybar, xbar
!---------------------------------------------
! This algorithm goes for accuracy without worrying about memory requirements.
! ydata:   dependent variable data matrix
! xdata:   exogenous variable data matrix
! lags:    number of lags
! breaks:  rows in ydata and xdata after which there is a break.  This allows for
!          discontinuities in the data (e.g. war years) and for the possibility of
!          adding dummy observations to implement a prior.  This must be a column vector.
!          Note that a single dummy observation becomes lags+1 rows of the data matrix,
!          with a break separating it from the rest of the data.  The function treats the 
!          first lags observations at the top and after each "break" in ydata and xdata as
!          initial conditions. 
! lambda:  weight on "co-persistence" prior dummy observations.  This expresses
!          belief that when data on *all* y's are stable at their initial levels, they will
!          tend to persist at that level.  lambda=5 is a reasonable first try.  With lambda<0,
!          constant term is not included in the dummy observation, so that stationary models
!          with means equal to initial ybar do not fit the prior mean.  With lambda>0, the prior
!          implies that large constants are unlikely if unit roots are present.
! mu:      weight on "own persistence" prior dummy observation.  Expresses belief
!          that when y_i has been stable at its initial level, it will tend to persist
!          at that level, regardless of the values of other variables.  There is
!          one of these for each variable.  A reasonable first guess is mu=2.
!      The program assumes that the first lags rows of ydata and xdata are real data, not dummies.
!      Dummy observations should go at the end, if any.  If pre-sample x's are not available,
!      repeating the initial xdata(lags+1,:) row or copying xdata(lags+1:2*lags,:) into 
!      xdata(1:lags,:) are reasonable subsititutes.  These values are used in forming the
!      persistence priors.
! Code written by Christopher Sims.  This version 6/15/03.
![T,nvar]=size(ydata);
T=size(ydata,1)
nvar=size(ydata,2)
!nox=isempty(xdata);
isx=present(xdata)
if(isx) then
   T2=size(xdata,1)
   nx=size(xdata,2)
else
   T2=T
   nx=0
endif
dim xbar(1,nx), ybar(1,nvar)
% note that x must be same length as y, even though first part of x will not be used.
% This is so that the lags parameter can be changed without reshaping the xdata matrix.
if(T2 ~= T) then
   print *,'Mismatch of x and y data lengths')
   return
endif
if (.not.present(breaks)) then
   nbreaks=0
else
   nbreaks=size(breaks);
end
breaks=[0,breaks,T];
do nb=1, nbreaks+1
   smpl=[smpl,[breaks(nb)+lags+1:breaks(nb+1)]]
end do
Tsmpl=size(smpl)
dim X(Tsmpl,nvar,lags)
do is=1,Tsmpl
    X(is,:,:)=transpose(ydata(smpl(is)-1:smpl(is)-lags),:))
end do
X=reshape(source=[X xdata(smpl,:)], shape=[Tsmpl,nvar*lags+size(xdata,2)])
y=ydata(smpl,:)
% Everything now set up with input data for y=Xb+e 
% ------------------Form persistence dummies-------------------
if (lambda /= 0 .or. mu>0) then
   ybar=sum(ydata(1:lags,:),1)/lags
   if (isx )
      xbar=sum(xdata(1:lags,:),1)/lags
   else
      xbar = []
   endif
   if (lambda /= 0) then
      if (lambda<0) then
         lambda=-lambda
         xbar = 0.0  
      endif
      xdum=lambda*[reshape(source=[],shape=[1,nvar*lags],pad=ybar), xbar]
      ydum=lambda*ybar;
      y=[y;ydum];
      X=[X(:,:);xdum];
   end
   if mu>0
      xdum=[repmat(diag(ybar),1,lags) zeros(nvar,nx)]*mu;
      ydum=mu*diag(ybar);
      X=[X;xdum];
      y=[y;ydum];
   end
end
[vl,d,vr]=svd(X(:,:),0);
di=1../diag(d);
B=vl'*y;
B=(vr.*repmat(di',nvar*lags+nx,1))*B;
u=y-X(:,:)*B;
xxi=vr.*repmat(di',nvar*lags+nx,1);
xxi=xxi*xxi';
B=reshape(B,[nvar*lags+nx,nvar]); % rhs variables, equations
By=B(1:nvar*lags,:);
By=reshape(By,nvar,lags,nvar);% variables, lags, equations
By=permute(By,[3,1,2]); %equations, variables, lags to match impulsdt.m
if nox
   Bx=[];
else
   Bx=B(nvar*lags+(1:nx),:)';
end
%logintlh=matrictint(u'*u,xxi,size(X,1)-nvar-1)-.5*nvar*(nvar+1)*log(2*pi);
var.By=By;var.Bx=Bx;var.u=u;var.xxi=xxi;%var.logintlh=logintlh;  
% Desired features: 1) automatic dummies for vcv prior
%                   2) automatic calculation of integrated pdf, accounting
%                      for the dummy variables as a prior
%                   3) automatic dummies for "Minnesota prior"
