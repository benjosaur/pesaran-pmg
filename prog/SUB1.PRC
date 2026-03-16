/*  **********************************************************
*   Procedures used in the ML estimation option 1 for the main
*   program "MAIN.PRG", which estimates the pooled ML dynamic
*   heterogenous panels.
*   ********************************************************** */

/*  ********************************************************  */
    Proc (6) = OLS1(Y,X,Z,PLAG,QLAG,k,N,T,yname,xname,zname);
/*  ********************************************************
*   Estimating the individual unrestricted ECM model
*   to obtain the OLS estimators of the parameters,
*   then obtain the estimators of the LR parameters, theta
*
*   dy(i) = a(i) + phi(i)*y(i,-1) + r_1(i)*dy(i,-1)
*         + ... + r_p-1(i)*dy(i,-p+1)
*         + b(i)*x(i) + c_0(i)*dx(i) + c_1(i)*cx(i,-1)
*         + ... + c_q-1(i)*dx(i,-q+1) + e(i),
*
*   where  i = 1, ..., N,
*
*  Also compute the sum of unrestricted log-likelihood functions
*
*   Y:  SUMT-vector of dependent variable
*   X: SUMTxk  matrix of stochastic independent variables
*   Z: SUMTxnz matrix of fixed regressors
*
*   PLAG: lag order of the dependent variable:  Nx1
*   QLAG: lag orders of the independent variables, X: Nxk
*   ************************************************************* */

Local theta,thta,sumlike,ini,fini,mxi,Tis,mxj,nreg,
      ytemp,xtemp,numsi,y0,y1,dy,x0,ww,inpt,sig2,reg,
      tmp,tmp1,tmp2,tmp3,tmp4,tmp5,SRb,faihat,like,
      varb,stdeb,tb,olsb,olsseb,olstb,olsres,beta,
      dtheta,h,vtheta,setheta,ttheta,
      dep,yhat,res,olsdres,ecm,epara,dtheta1,dtheta2,
      coname,cfname,ztemp;

Theta = zeros(k,n);     @ Theta @
OLSb = 0;                 @ initialization of the OLS estimators, s.e and t @
OLSseb = 0;
OLStb = 0;
sumlike = 0;              @ restricted sum of log-likelihood @
numsi = zeros(n,1);       @ number of other parameters for each group @
OLSDRES = zeros(11,n);    @ Diagnostic Results including R2, LL and so on @
CFname = 0;               @ Name of all the LR and SR parameters @

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;

   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ******************************************
*  Generation of the data used in estimation
*
*  y0: current value of dep vble
*  dy: first difference of dep vble
*  y1: lagged value of dep vble
*  x0: current value of indep vbles, x
*  ww: dy's, dx1's, dx2's and fixed regressor including intercept
*  ************************************************************** */

   ytemp = y[ini:fini];     @ y for ith group @
   xtemp = x[ini:fini,.];    @ x for ith group @
   ztemp = z[ini:fini,.];    @ x for ith group @

   {y0, dy, y1, x0, ww, coname}
   = NDGP1(ytemp,xtemp,ztemp,plag,qlag,k,mxi,tis,yname,xname,zname);

   numsi[mxi] = cols(ww);  @ number of other parameters @
   CFname = CFname|Coname;

/* **************************************
*  The OLS Estimation of the UR ECM Model
*  ************************************** */

@  Estimation of the OLS parameters @

   if plag[mxi] eq 0;  @ no lagged dep vble @
      reg = x0~ww;     @ regressors x, (differences), fixed regressors @
      dep = y0;
   else;
      reg = y1~x0~ww; @ regressors, y(-1), x1, (differences), fixed regressors @
      dep = dy;
   endif;

   tmp1 = invpd(reg'reg);
   SRb  = tmp1*reg'dep;          @ OLS parameters @

   tmp2 = (dep - reg*SRb);      @ residual @
   sig2 = tmp2'tmp2/(rows(dep)-cols(reg));  @ e'e/(T-k) @

   varb  = sig2.*tmp1;
   stdeb = sqrt(diag(varb));
   tb    = srb./stdeb;

   sig2 = tmp2'tmp2/rows(dep);  @ e'e/T @
   like = -(rows(dy)/2)*ln(2*pi*sig2) - rows(dy)/2;  @ log-likelihood @
   sumlike = sumlike + like;  @ sum of log-likelihood @

@  computation of the LR parameters @

   if plag[mxi] eq 0;     @ no lagged dep vble @

      beta = SRb[1:k];        @ beta(i) @
      thta = beta;

      setheta = stdeb[1:k];      @ std. err. of theta @
      ttheta  = thta./setheta;   @ t-ratios of theta @

   else;

      faihat = SRb[1];          @ fai(i) @
      beta = SRb[2:k+1];        @ beta(i) @

      thta = -beta./faihat;     @ theta(i) @

      @  Computation of the standard errors and t-values
        for the LR parameters using the delta method @

      dtheta1 = zeros(k,1);  @ first differential  w.r.t. phi @
      h = 1;
      do until h > k;
         dtheta1[h,1] = beta[h]/faihat^2;
         h = h + 1;
      endo;

      dtheta2 = eye(k);     @ first differential w.r.t beta @
      dtheta2 = -(1/faihat).*dtheta2;
      dtheta  = dtheta1~dtheta2;

      vtheta  = dtheta*varb[1:k+1,1:k+1]*dtheta'; @ variance of theta @
      setheta = sqrt(diag(vtheta));                 @ std. err. of theta @
      ttheta  = thta./setheta;                     @ t-ratios of theta @

   endif;

   Theta[.,mxi] = thta;

@  computation of R2 and diagnostic statistics using UR ECM @

   dep = dy;

   if plag[mxi] eq 0; epara = -1|srb;
                      nreg = k + numsi[mxi]; @ # of regressors @
   else;              epara = srb;
                      nreg = k + 1 + numsi[mxi]; @ # of regressors @
   endif;

   reg = y1~x0~ww;
   yhat = reg*epara;
   res  = tmp2;

   { OLSDRES[.,mxi] } = DIAGNOS(dep,reg,res,yhat,nreg);

@  Saving the results for the estimates, std. errors and t-ratios
*  in the order of theta, fai, beta and si:
*  the dimension is k+1+k+numsi(i) @

   if plag[mxi] eq 0;

      olsb   = olsb|thta|-1|srb;
      olsseb = olsseb|setheta|"NA"|stdeb;
      olstb  = olstb|ttheta|"NA"|tb;
   else;
      olsb   = olsb|thta|srb;
      olsseb = olsseb|setheta|stdeb;
      olstb  = olstb|ttheta|tb;
   endif;

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

OLSRES = olsb~olsseb~olstb;       @ Results for estimates, std. err and t @

@ Result in the order of theta1, theta2, phi, beta1, beta2, other para
* for group 1, 2, ..., N:
* # of parameters in each group is k+1+k+numsi(i) @

OLSRES = olsres[2:rows(olsres),.];
CFname = CFname[2:rows(CFname),.];

RETP(OLSRES,OLSDRES,theta,sumlike,numsi,cfname);ENDP;

/*  ***************************************************  */
    Proc (2) = PMLEBS1(Y,X,Z,PLAG,QLAG,k,N,T,OLDTHETA);
/*  ***************************************************
*   Computing the Pooled ML estimatiors of the LR parameters as well as
*   the SR parameters given the LR parameters by `BACK-SUBSTITUTION'
*   ************************************************************* */

Local sumden,sumnum,sumlike,ini,fini,mxi,Tis,mxj,
      ytemp,xtemp,y0,y1,dy,x0,ww,inpt,sig2,reg,
      tmp,tmp1,tmp2,tmp3,tmp4,tmp5,SRb,faihat,beta2hat,T1den,T1num,like,
      theta,tmp0,ztemp;

sumden = zeros(k,k);    @ denominator of Theta1 @
sumnum = zeros(k,1);     @ numerator of Theta1 @
sumlike = 0;              @ restricted sum of log-likelihood @

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;

   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ******************************************
*  Generation of the data used in estimation
*
*  y0: current value of dep vble
*  dy: first difference of dep vble
*  y1: lagged value of dep vble
*  x0: current value of indep vbles, x
*  ww: dy's, dx1's, dx2's and intercept
*  ****************************************** */

   ytemp = y[ini:fini];     @ y for ith group @
   xtemp = x[ini:fini,.];   @ x for ith group @
   ztemp = z[ini:fini,.];   @ z for ith group @

   {y0, dy, y1, x0, ww}
   = DGP1(ytemp,xtemp,ztemp,plag,qlag,k,mxi,tis);

/* **********************************************
*  ML Estimation of the SR and the LR parameters
*  *********************************************** */

@  Estimation of phi and other SR parameters including fixed regressors @

   if plag[mxi] eq 0;  faihat = -1;
                       tmp0 = y0 - x0*oldtheta;
                       SRb = invpd(ww'ww)*ww'tmp0;
                       tmp2 = tmp0 - ww*SRb;      @ residual @
   else;
                       tmp = y1 - x0*oldtheta;
                       reg = tmp~ww;
                       tmp1 = invpd(reg'reg);
                       SRb = tmp1*reg'dy;
                       faihat = SRb[1];            @ phi @
                       tmp2 = dy - reg*SRb;      @ residual @
   endif;

   sig2 = tmp2'tmp2/rows(dy);  @ e'e/T @

   like = -(rows(dy)/2)*ln(2*pi*sig2) - rows(dy)/2;  @ log-likelihood @
   sumlike = sumlike + like;  @ sum of log-likelihood @

@  Estimation of theta @

   tmp3  = x0'x0 - x0'ww*invpd(ww'ww)*ww'x0;
   T1den = (faihat^2/sig2).*tmp3;
   tmp4  = dy - faihat*y1;
   tmp5  = x0'tmp4 - x0'ww*invpd(ww'ww)*ww'tmp4;
   T1num = -(faihat/sig2).*tmp5;

   sumden = sumden + T1den;   @ sum of denominator @
   sumnum = sumnum + T1num;   @ sum of numerator @

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

theta = invpd(sumden)*sumnum;

RETP(theta,sumlike);ENDP;

/*  **********************************************************  */
    Proc (3) = PMLENR1(Y,X,Z,PLAG,QLAG,k,N,T,OLDTHETA,OLDPHI);
/*  **********************************************************
*   Computing the Pooled ML estimatiors of the LR parameters
*   and the ECM parameters by `NEWTON-RAPHSON'
*   ****************************************** */

Local sumden,sumnum,sumlike,ini,fini,mxi,Tis,mxj,ytemp,xtemp,
      y0,y1,dy,x0,ww,inpt,sig2,reg,tmp,tmp1,tmp2,tmp3,tmp4,tmp5,
      jnk1,xi,phifd,phisd,thphisd,phi,oldpara,ndep,res,para,faihat,
      Tden1,Tden2,Tden,Tnum,like,theta,tmp0,ztemp;

sumden = zeros(k,k);     @ k by k matrix of the second der for Theta block @
sumnum = zeros(k,1);     @ k by 1 vectot with respect to Theta @

phifd = zeros(n,1);      @ n by 1 vector with respect to phi_1,..,phi_N @
phisd = zeros(n,n);      @ n by n matrix of the second der for phi_1,..,phi_N @

thphisd = zeros(k,n);    @ k by n matrix with respect to theta and then phi @

sumlike = 0;             @ restricted sum of log-likelihood @

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;

   faihat = oldphi[mxi];

   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ******************************************
*  Generation of the data used in estimation
*
*  y0: current value of dep vble
*  dy: first difference of dep vble
*  y1: lagged value of dep vble
*  x0: current value of indep vbles, x
*  ww: dy's, dx1's, dx2's and intercept
*  ****************************************** */

   ytemp = y[ini:fini];     @ y for ith group @
   xtemp = x[ini:fini,.];   @ x for ith group @
   ztemp = z[ini:fini,.];   @ z for ith group @

   {y0, dy, y1, x0, ww}
   = DGP1(ytemp,xtemp,ztemp,plag,qlag,k,mxi,tis);

/* **********************************************
*  ML Estimation of the SR and the LR parameters
*  *********************************************** */

@  Construction of the residual, sigma, error correction term
   given theta(j-1) and phi(j-1) @

   xi = y1 - x0*oldtheta;                  @ ec term @
   ndep = dy - xi*faihat;             @ dy - phi*ecm @
   res = ndep - ww*invpd(ww'ww)*ww'ndep;   @ residual @
   sig2 = res'res/rows(dy);  @ sigma^2 @

   like = -(rows(dy)/2)*ln(2*pi*sig2) - rows(dy)/2;  @ log-likelihood @
   sumlike = sumlike + like;  @ sum of restricted log-likelihood @

@  Estimation of theta and phi by NR @
@  Construction of the gradient vector @

   jnk1 = x0'res - x0'ww*invpd(ww'ww)*ww'res;   @ X(i)'Hi*e(i) @
   jnk1 = -(faihat/sig2).*jnk1;
   sumnum = sumnum + jnk1;     @ gradient for theta @

   jnk1 = xi'res - xi'ww*invpd(ww'ww)*ww'res;   @ xi(i)'Hi*e(i) @
   phifd[mxi] = (1/sig2)*jnk1; @ gradient for phi @

@  Construction of the Hesssian Matrix @

   jnk1  = x0'x0 - x0'ww*invpd(ww'ww)*ww'x0;   @ X(i)'Hi*X(i) @
   jnk1 = (faihat^2/sig2).*jnk1;

   sumden = sumden + jnk1;   @ for theta @

   jnk1 = xi'xi - xi'ww*invpd(ww'ww)*ww'xi;  @ xi(i)'Hi*xi(i) @
   phisd[mxi,mxi] = (1/sig2)*jnk1;

   jnk1 = x0'xi - x0'ww*invpd(ww'ww)*ww'xi;  @ X(i)'Hi*xi(i) @
   thphisd[.,mxi] = -(faihat/sig2).*jnk1;

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

TNUM  = sumnum|phifd;      @ k+n by 1 @
TDEN1 = sumden~thphisd;    @ k by (k+n) @
TDEN2 = thphisd'~phisd;    @ n by (k+n) @
TDEN = TDEN1|Tden2;        @ (k+n) by (k+n) @

oldpara = oldtheta|oldphi;
para = oldpara + invpd(TDEN)*Tnum;
theta = para[1:k];
phi = para[k+1:rows(para)];

RETP(theta,phi,sumlike);ENDP;

/*  ***************************************************  */
    Proc (1) = INIPHI1(Y,X,Z,PLAG,QLAG,k,N,T,OLDTHETA);
/*  **********************************************************
*   Computing the initial estimates of the ECM parameter
*   which will be used in `NEWTON-RAPHSON'
*   ****************************************** */

Local phiini,ini,fini,mxi,Tis,mxj,ytemp,xtemp,ztemp,xi,den,num,
      y0,y1,dy,x0,ww,inpt,sig2,reg,tmp,tmp1,tmp2,tmp3,tmp4,tmp5;

phiini = zeros(n,1);     @ n by 1 vector of initial estimates of phi @

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;

   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ******************************************
*  Generation of the data used in estimation
*
*  y0: current value of dep vble
*  dy: first difference of dep vble
*  y1: lagged value of dep vble
*  x0: current value of indep vbles, x
*  ww: dy's, dx1's, dx2's and intercept
*  ****************************************** */

   ytemp = y[ini:fini];     @ y for ith group @
   xtemp = x[ini:fini,.];   @ x for ith group @
   ztemp = z[ini:fini,.];   @ z for ith group @

   {y0, dy, y1, x0, ww}
   = DGP1(ytemp,xtemp,ztemp,plag,qlag,k,mxi,tis);

@  The error correction term given theta(j-1) @

   xi = y1 - x0*oldtheta;                  @ ec term @

@  Estimates of phi @

   den = xi'xi - xi'ww*invpd(ww'ww)*ww'xi;   @ denominator @
   num = xi'dy - xi'ww*invpd(ww'ww)*ww'dy;   @ numerator @

   if plag[mxi] eq 0; phiini[mxi] = -1;
   else;              phiini[mxi] =  num/den;
   endif;

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

RETP(phiini);ENDP;

/*  ******************************************************  */
    Proc (4) = PMLESR1(Y,X,Z,PLAG,QLAG,k,N,T,NUMSI,THETA);
/*  ******************************************************
*   Final estimation results using the converged theta1
*   *************************************************** */

Local sumlike,ini,fini,mxi,Tis,p_cnt,q1_cnt,q2_cnt,mxj,sumsi,
      g12ini,g12fini,g13ini,g13fini,g14ini,g14fini,gg1,gg2,
      g2ini,g2fini,rg23ini,rg23fini,cg23ini,cg23fini,
      rg24ini,rg24fini,cg24ini,cg24fini,g4ini,g4fini,
      g3ini,g3fini,rg34fini,rg34ini,cg34ini,cg34fini,
      g11,g12,g13,g14,g22,g23,g24,g33,g34,g44,a2tmp,a3tmp,a4tmp,tmp0,
      g1tmp,g12tmp,g13tmp,g14tmp,g2tmp,g23tmp,g24tmp,g3tmp,g34tmp,g4tmp,
      ytemp,xtemp,y0,y1,dy,x0,ww,inpt,sig2,reg,
      tmp,tmp1,tmp2,tmp3,tmp4,tmp5,SRb,faihat,beta2hat,T1den,T1num,like,
      allpara,se_all,t_all,gg,si,
      beta,btmp,sebtmp,tbtmp,ini12,ini13,ini22,ini23,ini33,
      fini12,fini13,fini22,fini23,fini33,kk,dbeta,dbeta1,dbeta2,h,ggtmp,
      vtmp11,vtmp12,vtmp13,vtmp22,vtmp23,vtmp33,vtmp1,vtmp2,vtmp3,vtmp,
      vbeta,sebeta,tbeta,sebta,tbta,seall,tall,pmleres,PMLECOV,
      chk2,dtmp,pmledres,bta,nreg,dep,yhat,res,betmp,ztemp,
      epara,ecm,na,nfai0,nn,a3ntmp,hh,atmp,htmp,hflag,a3se,a3t,a3temp,hini;

PMLEDRES = zeros(11,n);

beta = zeros(k,n);  @ SR parameters @

@ initialization for phi and si @
a3tmp = 0;
a4tmp = 0;

nfai0 = 0;  @ number of phi(i) = -1 @

sumsi = sumc(numsi); @ sum of other parameters for all groups @
                     @ note that numsi could be different!! @

@ Information matrix G @

@ Theta1 block @

G11 = zeros(k,k);
G13 = zeros(k,n);
G14 = zeros(k,sumsi);

@ Fai block @

G33 = zeros(n,n);
G34 = zeros(n,sumsi);

@ Si block @

G44 = zeros(sumsi,sumsi);

@ initial and end positions @

g13fini = 0;     g13ini = 1;
g14fini = 0;     g14ini = 1;

g3fini   = 0;      g3ini = 1;
rg34fini = 0;    rg34ini = 1;
cg34fini = 0;    cg34ini = 1;

g4fini = 0;      g4ini = 1;

sumlike = 0;              @ restricted sum of log-likelihood @

ini  = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;

   Tis  = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

   ytemp = y[ini:fini];     @ y for ith group @
   xtemp = x[ini:fini,.];   @ x for ith group @
   ztemp = z[ini:fini,.];   @ z for ith group @

   {y0, dy, y1, x0, ww}
   = DGP1(ytemp,xtemp,ztemp,plag,qlag,k,mxi,tis);

/* *********************************************
*  ML Estimation of the SR and the LR parameters
*  ********************************************* */

@  Estimation of phi @

   if plag[mxi] eq 0;  faihat = -1;
                       nfai0 = nfai0 + 1;
                       tmp0 = y0 - x0*theta;
                       si = invpd(ww'ww)*ww'tmp0;
                       tmp2 = tmp0 - ww*si;      @ residual @
                       nreg = k + numsi[mxi];    @ # of regressors @

   else;               tmp = y1 - x0*theta;
                       reg = tmp~ww;
                       tmp1 = invpd(reg'reg);
                       SRb = tmp1*reg'dy;
                       faihat = SRb[1];            @ phi @
                       si = SRb[2:rows(srb)];
                       tmp2 = dy - reg*SRb;      @ residual @
                       nreg = k + 1 + numsi[mxi];    @ # of regressors @
   endif;

   sig2 = tmp2'tmp2/rows(dy);  @ e'e/T @

   like = -(rows(dy)/2)*ln(2*pi*sig2) - rows(dy)/2;  @ log-likelihood @
   sumlike = sumlike + like;  @ sum of log-likelihood @

@  computation of the SR parameters: beta @

   betmp = -faihat.*theta;
   beta[.,mxi] = betmp;        @ beta(i) @

@  Generation of the final ecm term: y(-1) - x*theta @

   ecm = y1 - x0*theta;

@  computation of R2 and diagnostic statistics using UR ECM @

   if plag[mxi] eq 0;  faihat = -1;
                       nreg = k + numsi[mxi];    @ # of regressors @
   else;               faihat = SRb[1];
                       nreg = k + 1+ numsi[mxi];    @ # of regressors @
   endif;
   dep = dy;
   reg = y1~x0~ww;
   epara = faihat|betmp|si;
   yhat = reg*epara;
   res = dep - yhat;

   { PMLEDRES[.,mxi] } = DIAGNOS(dep,reg,res,yhat,nreg);

@  Estimation of the Covariance matrix for all parameters @

   a3tmp = a3tmp|faihat;    @ stacking fai(i) @
   a4tmp = a4tmp|si;        @ stacking si(i) @

@ Theta1 block of information matrix @

   G1tmp = (faihat^2/sig2).*(x0'x0);
   G11 = G11 + G1tmp;

   G13tmp  = -(faihat/sig2).*(x0'ecm);
   g13fini = g13fini + 1;
   G13[.,g13ini:g13fini] = G13tmp;
   g13ini  = g13ini + 1;

   G14tmp  = -(faihat/sig2).*(x0'ww);
   g14fini = g14fini + numsi[mxi];
   G14[.,g14ini:g14fini] = G14tmp;
   g14ini  = g14ini + numsi[mxi];

@ Fai block of information matrix @

   G3tmp  = (1/sig2).*(ecm'ecm);
   g3fini = g3fini + 1;
   if plag[mxi] eq 0; G33[g3ini:g3fini,g3ini:g3fini] = 0; @ no varaince @
   else;              G33[g3ini:g3fini,g3ini:g3fini] = G3tmp;
   endif;
   g3ini  = g3ini + 1;

   G34tmp   = (1/sig2).*(ecm'ww);
   rg34fini = rg34fini + 1;
   cg34fini = cg34fini + numsi[mxi];
   G34[rg34ini:rg34fini,cg34ini:cg34fini] = G34tmp;
   rg34ini  = rg34ini + 1;
   cg34ini  = cg34ini + numsi[mxi];

@ Si block of information matrix @

   G4tmp  = (1/sig2).*(ww'ww);
   g4fini = g4fini + numsi[mxi];
   G44[g4ini:g4fini,g4ini:g4fini] = G4tmp;
   g4ini  = g4ini + numsi[mxi];

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

@ delete unnecessary initial values @

a3tmp = a3tmp[2:rows(a3tmp)];  @ n by 1    @
a4tmp = a4tmp[2:rows(a4tmp)];  @ sumsi by 1 @

nn = n - nfai0;        @ number of cases where phi =\ -1 @
if nn eq 0;
else;   a3ntmp = zeros(1,2);
        @ initialization for nnx2 matrix of phi(i) and i @
        hh = 1;
        do until hh > n;

           if a3tmp[hh] eq -1;
           else; atmp = a3tmp[hh]~hh;
                 a3ntmp = a3ntmp|atmp;
           endif;
           hh = hh + 1;
        endo;
        a3ntmp = a3ntmp[2:rows(a3ntmp),.];
endif;

@ whole information matrix @

   g1tmp = G11~G13~G14;
   g3tmp = G13'~G33~G34;
   g4tmp = G14'~G34'~G44;

   GG = g1tmp|g3tmp|g4tmp;

@  Deal with the case where some of phi(i) = -1 @

   if nfai0 eq 0;
   else;

      hflag = 0;
GREPEAT:
      hh = k + 1;
      do until hh > k + n;

         htmp = GG[hh,hh];
         if htmp eq 0;   @ delete the block of information for phi(i) = -1 @
            GG1 = GG[1:hh-1,1:hh-1]~GG[1:hh-1,hh+1:cols(GG)];
            GG2 = GG[hh+1:rows(GG),1:hh-1]~GG[hh+1:rows(GG),hh+1:rows(GG)];
            GG = GG1|GG2;
            hflag = hflag + 1;
            goto DASI;
         endif;
         hh = hh + 1;

      endo;
DASI:
      if hflag eq nfai0;
      else; goto GREPEAT;
      endif;
   endif;

GGtmp = invpd(GG);
PMLECOV = GGtmp[1:k,1:k];

@ standard errors and t-values @

if nn eq 0; allpara = theta|a4tmp;
else;       allpara = theta|a3ntmp[.,1]|a4tmp;
endif;

se_all = sqrt(diag(GGtmp));
t_all = allpara./se_all;

@ including the case of no standard errors for phi(i) = -1 @

a3se = zeros(n,1);
a3t  = zeros(n,1);
if nn eq 0;
else;
   hini = k + 1;
   hh = 1;
   do until hh > nn;

      a3temp = a3ntmp[hh,2];
      a3se[a3temp] = se_all[hini];
      a3t[a3temp] = t_all[hini];

      hini = hini + 1;
      hh = hh + 1;
   endo;
endif;

@  Computation of the standard errors and t-values
   for the SR parameters (using the delta method) @

   sebeta = zeros(k,n);
   tbeta  = zeros(k,n);

   ini13 = k+1;    fini13 = k;
   ini33 = k+1;    fini33 = k;

   ini = 1;
   fini = 0;
   mxi = 1;
   do until mxi > n;

      faihat = a3tmp[mxi];

      if plag[mxi] eq 0;   @ no lagged dep vble @

         vbeta = GGtmp[1:k,1:k];
         @ Variance for theta as well as beta @

      else;

         vtmp11 = GGtmp[1:k,1:k];
         fini13 = fini13 + 1;
         vtmp13 = GGtmp[1:k,ini13:fini13];
         ini13  = ini13 + 1;
         fini33 = fini33 + 1;
         vtmp33 = GGtmp[ini33:fini33,ini33:fini33];
         ini33  = ini33 + 1;

         vtmp1 = vtmp11~vtmp13;
         vtmp3 = vtmp13'~vtmp33;
         vtmp  = vtmp1|vtmp3;    @ variance for theta and phi @

         dbeta1 = eye(k);     @ first differential w.r.t beta @
         dbeta1 = -faihat.*dbeta1;

         dbeta2 = zeros(k,1);  @ first differential  w.r.t. phi @
         h = 1;
         do until h > k;
            dbeta2[h,1] = -theta[h];
            h = h + 1;
         endo;
         dbeta = dbeta1~dbeta2;

         vbeta = dbeta*vtmp*dbeta';   @ variance of beta @

      endif;

      bta = beta[.,mxi];
      sebta = sqrt(diag(vbeta));  @ std. err. of beta @
      tbta = bta./sebta;        @ t-ratios of beta @

@  Saving the results for the estimates, std. errors and t-ratios @

      sebeta[.,mxi] = sebta;
      tbeta[.,mxi] = tbta;

   mxi = mxi + 1;
endo;

beta = vec(beta);
sebeta = vec(sebeta);
tbeta = vec(tbeta);

@ All parameters in the order of theta1, theta2(1),...theta2(N),
  phi(1), ..., phi(N), and si(1), ..., si(N):
  The dimension is (k1+k2*n)+n+(k1+k2*n)+sumnumsi-vector @

allpara = theta|a3tmp|beta|a4tmp;
fini = k;
seall = se_all[1:fini]|a3se|sebeta|se_all[fini+nn+1:rows(se_all)];
tall = t_all[1:fini]|a3t|tbeta|t_all[fini+nn+1:rows(t_all)];

PMLERES = allpara~seall~tall;

RETP(PMLERES,PMLEDRES,sumlike,PMLECOV);ENDP;

/* ******************************************************* */
   PROC (5) = DGP1(ytemp,xtemp,ztemp,plag,qlag,k,mxi,tis);
/* *******************************************************
*  Generation of the data used in estimation
*
*  y0: current value of dep vble
*  dy: first difference of dep vble
*  y1: lagged value of dep vble
*  x0: current value of indep vbles, x
*  ww: dy's, dx1's and intercept
*  ************************************* */

Local y1temp,dytemp,x1temp,dxtemp,dy1temp,dy2temp,dy3temp,dystemp,wflag,
      dx0temp,dx1temp,dx2temp,dx3temp,dxtmp,dxflag,
      ltemp,mxlag,sini,sig2,mxj,y0,y1,dy,x0,ww;

 @ Generation y, dy, y(-1), x and ww @

   y1temp  = lag(ytemp);       @ y(-1) @
   dytemp  = ytemp - y1temp;   @ dy @

   x1temp  = lag(xtemp);       @ x(-1) @
   dxtemp  = xtemp - x1temp;   @ dx @

 @ Generation of W such as lagged differences and intercept @

 @ Generation of lagged dy's @

   if     plag[mxi] eq 2; dy1temp = lag(dytemp);
                          dystemp = dy1temp;
   elseif plag[mxi] eq 3; dy1temp = lag(dytemp); dy2temp = lag(dy1temp);
                          dystemp = dy1temp~dy2temp;
   elseif plag[mxi] eq 4; dy1temp = lag(dytemp); dy2temp = lag(dy1temp);
                          dy3temp = lag(dy2temp);
                          dystemp = dy1temp~dy2temp~dy3temp;
   endif;

   if plag[mxi] > 1;  WW = dystemp; wflag = 1;
   else;              wflag = 0;    @ W's empty yet @
   endif;

 @ Generation of current and lagged dx's @

   dxtmp = zeros(Tis,1);
   mxj = 1;
   do until mxj > k;

      if     QLAG[mxi,mxj] eq 1;  dx0temp = dxtemp[.,mxj];
                                  dxtmp = dxtmp~dx0temp;
      elseif QLAG[mxi,mxj] eq 2;  dx0temp = dxtemp[.,mxj];
                                  dx1temp = lag(dx0temp);
                                  dxtmp = dxtmp~dx0temp~dx1temp;
      elseif QLAG[mxi,mxj] eq 3;  dx0temp = dxtemp[.,mxj];
                                  dx1temp = lag(dx0temp);
                                  dx2temp = lag(dx1temp);
                          dxtmp = dxtmp~dx0temp~dx1temp~dx2temp;
      elseif QLAG[mxi,mxj] eq 4;  dx0temp = dxtemp[.,mxj];
                                  dx1temp = lag(dx0temp);
                                  dx2temp = lag(dx1temp);
                                  dx3temp = lag(dx2temp);
                  dxtmp = dxtmp~dx0temp~dx1temp~dx2temp~dx3temp;
      endif;

      mxj = mxj + 1;
   endo;

   if cols(dxtmp) eq 1;  dxflag = 0;       @ No dx's in the regression @
   else; dxtmp = dxtmp[.,2:cols(dxtmp)];
         @ column dimension is q_1(i) + q_2(i) + ... + q_k(i) @
         dxflag = 1;
   endif;

   if dxflag eq 0;
   else;   if wflag eq 1;  WW = WW~dxtmp;
           else;           WW = dxtmp;
           endif;
           wflag = 1;  @ W's not empty anymore @
   endif;

   if wflag eq 0; WW = ztemp;
   else;          WW = WW~ztemp;
   endif;

 @ sample adjustment @

   ltemp = plag[mxi]|qlag[mxi,.]';
   mxlag = maxc(ltemp);      @ maximum lag order for each group @

   if mxlag eq 0; sini = 2;
   else;          sini = mxlag + 1;
   endif;

   y0 = ytemp[sini:Tis];        @ y @
   y1 = y1temp[sini:Tis];       @ y(-1) @
   dy = dytemp[sini:Tis];       @ dy @

   x0 = xtemp[sini:Tis,.];      @ x @
   ww = ww[sini:Tis,.];         @ other variables including intercept @

RETP(y0,dy,y1,x0,ww); ENDP;

/* ************************************************************************** */
   PROC (6) = NDGP1(ytemp,xtemp,ztemp,plag,qlag,k,mxi,tis,yname,xname,zname);
/* **************************************************************************
*  See DGP1
*  ******** */

Local y1temp,dytemp,x1temp,dxtemp,dy1temp,dy2temp,dy3temp,dystemp,wflag,
      dx0temp,dx1temp,dx2temp,dx3temp,dxtmp,dxflag,
      ltemp,mxlag,sini,sig2,mxj,y0,y1,dy,x0,ww,s,
      dyname,dxname,dxtmp0,dxtmp1,srname,srname1,lrname1,coname;

 @ Generation y, dy, y(-1), x and ww @

   y1temp  = lag(ytemp);       @ y(-1) @
   dytemp  = ytemp - y1temp;   @ dy @

   x1temp  = lag(xtemp);       @ x(-1) @
   dxtemp  = xtemp - x1temp;   @ dx @

 @ Generation of W such as lagged differences and intercept @

 @ Generation of lagged dy's @

   SRname = "ini";
   if     plag[mxi] eq 2; dy1temp = lag(dytemp);
                          dystemp = dy1temp;
                          s = seqa(1,1,1);
                dyname = 0 $+ "d" $+ yname $+ "(-" $+ ftocv(s,1,0) $+ ")";
                          SRname = SRname|dyname;
   elseif plag[mxi] eq 3; dy1temp = lag(dytemp); dy2temp = lag(dy1temp);
                          dystemp = dy1temp~dy2temp;
                          s = seqa(1,1,2);
                dyname = 0 $+ "d" $+ yname $+ "(-" $+ ftocv(s,1,0) $+ ")";
                          SRname = SRname|dyname;
   elseif plag[mxi] eq 4; dy1temp = lag(dytemp); dy2temp = lag(dy1temp);
                          dy3temp = lag(dy2temp);
                          dystemp = dy1temp~dy2temp~dy3temp;
                          s = seqa(1,1,3);
                dyname = 0 $+ "d" $+ yname $+ "(-" $+ ftocv(s,1,0) $+ ")";
                          SRname = SRname|dyname;
   endif;

   if plag[mxi] > 1;
      WW = dystemp; wflag = 1;
   else; wflag = 0;    @ W's empty yet @
   endif;

 @ Generation of current and lagged dx1's @

   dxtmp = zeros(Tis,1);
   mxj = 1;
   do until mxj > k;

      if     QLAG[mxi,mxj] eq 1;
             dx0temp = dxtemp[.,mxj];
             dxtmp = dxtmp~dx0temp;

             dxname = 0 $+ "d" $+ xname[mxj];
             SRname = SRname|dxname;

      elseif QLAG[mxi,mxj] eq 2;
             dx0temp = dxtemp[.,mxj];
             dx1temp = lag(dx0temp);
             dxtmp = dxtmp~dx0temp~dx1temp;

             dxtmp0 = 0 $+ "d" $+ xname[mxj];
             s = seqa(1,1,1);
             dxtmp1 = 0 $+ "d" $+ xname[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      elseif QLAG[mxi,mxj] eq 3;
             dx0temp = dxtemp[.,mxj];
             dx1temp = lag(dx0temp);
             dx2temp = lag(dx1temp);
             dxtmp = dxtmp~dx0temp~dx1temp~dx2temp;

             dxtmp0 = 0 $+ "d" $+ xname[mxj];
             s = seqa(1,1,2);
             dxtmp1 = 0 $+ "d" $+ xname[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      elseif QLAG[mxi,mxj] eq 4;
             dx0temp = dx1temp[.,mxj];
             dx1temp = lag(dx0temp);
             dx2temp = lag(dx1temp);
             dx3temp = lag(dx2temp);
             dxtmp = dxtmp~dx0temp~dx1temp~dx2temp~dx3temp;

             dxtmp0 = 0 $+ "d" $+ xname[mxj];
             s = seqa(1,1,3);
             dxtmp1 = 0 $+ "d" $+ xname[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      endif;

      mxj = mxj + 1;
   endo;

   if cols(dxtmp) eq 1;  dxflag = 0;       @ No dx's in the regression @
   else; dxtmp = dxtmp[.,2:cols(dxtmp)];
         @ column dimension is q_1(i) + q_2(i) + ... + q_k(i) @
         dxflag = 1;
   endif;

   if dxflag eq 0;
   else;   if wflag eq 1;  WW = WW~dxtmp;
           else;           WW = dxtmp;
           endif;
           wflag = 1;  @ W's not empty anymore @
   endif;

   if wflag eq 0; WW = ztemp;
   else;          WW = WW~ztemp;
   endif;
   SRname = SRname|Zname;
   SRname = SRname[2:rows(SRname)];

 @ sample adjustment @

   ltemp = plag[mxi]|qlag[mxi,.]';
   mxlag = maxc(ltemp);      @ maximum lag order for each group @

   if mxlag eq 0; sini = 2;
   else;          sini = mxlag + 1;
   endif;

   y0 = ytemp[sini:Tis];        @ y @
   y1 = y1temp[sini:Tis];       @ y(-1) @
   dy = dytemp[sini:Tis];       @ dy @

   x0 = xtemp[sini:Tis,.];      @ x @
   ww = ww[sini:Tis,.];         @ other variables including intercept @

COname = Xname|"Phi"|Xname|SRname;

RETP(y0,dy,y1,x0,ww,Coname); ENDP;

