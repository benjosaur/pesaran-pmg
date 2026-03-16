/*  **********************************************************************  */
    Proc (7)
    = OLS2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,k1,k2,N,T,Yname,X1name,X2name,Zname);
/*  **********************************************************************
*   Estimating the individual unrestricted ECM model
*   to obtain the OLS estimators of the parameters,
*   then obtain the estimators of the LR parameters,
*   theta1(i) and theta2(i).
*
*   dy(i) = a(i) + phi(i)*y(i,-1) + r_1(i)*dy(i,-1)
*         + ... + r_p-1(i)*dy(i,-p+1)
*         + b(i)*x(i) + c_0(i)*dx(i) + c_1(i)*cx(i,-1)
*         + ... + c_q-1(i)*dx(i,-q+1) + e(i),
*
*   Y:  SUMT-vector of dependent variable
*   X1: SUMTxk1 matrix of independent variables
*      where the LR parameters are restricted to be the same
*   X2: SUMTxk2 matrix of independent variables
*      where the LR parameters are allowed to be different
*
*   PLAG: lag order of the dependent variable:  Nx1
*   Q1LAG: lag orders of the independent variables, X1: Nxk1
*   Q2LAG: lag orders of the independent variables, X2: Nxk2
*   ************************************************************* */

Local theta1,theta2,sumlike,ini,fini,mxi,Tis,mxj,ztemp,
      ytemp,x1temp,x2temp,numsi,y0,y1,dy,x10,x20,ww,inpt,sig2,reg,
      tmp,tmp1,tmp2,tmp3,tmp4,tmp5,SRb,faihat,beta2hat,T1den,T1num,like,
      varb,stdeb,tb,olsb,olsseb,olstb,olsres,t1hat,t2hat,beta,
      b1hat,b2hat,theta,kk,dtheta,h,vtheta,setheta,ttheta,
      dep,yhat,res,olsdres,ecm,epara,dtheta1,dtheta2,nreg,coname,cfname,
      dtmp,RSQ,RBARSQ,LL,AIC,SBC,HQ,chiSC,chiFF,chiNO,chiHE,sige;

Theta1 = zeros(k1,n);     @ Theta1 @
Theta2 = zeros(k2,n);     @ Theta2 @
OLSb = 0;                 @ initialization of the OLS estimators, s.e and t @
OLSseb = 0;
OLStb = 0;
sumlike = 0;              @ restricted sum of log-likelihood @
numsi = zeros(n,1);       @ number of other parameters for each group @
OLSDRES = zeros(11,n);    @ Diagnostic Results including R2, LL and so on @
CFname = "ini";

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;

   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ****************************************************
*  Generation of the data used in estimation: See NDGP2
*  **************************************************** */

   ytemp  = y[ini:fini];       @ y for ith group @
   x1temp = x1[ini:fini,.];    @ x1 for ith group @
   x2temp = x2[ini:fini,.];    @ x2 for ith group @
   ztemp  = z[ini:fini,.];     @ z for ith group @

   {y0, dy, y1, x10, x20, ww, COname} = NDGP2(ytemp,x1temp,x2temp,ztemp,
               plag,q1lag,q2lag,k1,k2,mxi,tis,yname,x1name,x2name,zname);

   CFname = CFname|COname;

   numsi[mxi] = cols(ww);  @ number of other parameters @

/* **************************************
*  The OLS Estimation of the UR ECM Model
*  ************************************** */

@  Estimation of the OLS parameters @

   if plag[mxi] eq 0;       @ no lagged dep vble @
      reg = x10~x20~ww;     @ regressors x1, x2, (differences), inpt @
      dep = y0;
   else;
      reg = y1~x10~x20~ww;  @ regressors, y(-1), x1, x2, (differences), inpt @
      dep = dy;
   endif;

   tmp1 = invpd(reg'reg);
   SRb = tmp1*reg'dep;          @ OLS parameters @

   tmp2 = (dep - reg*SRb);      @ residual @
   sig2 = tmp2'tmp2/(rows(dep)-cols(reg));  @ e'e/(T-m) @

   varb = sig2.*tmp1;
   stdeb = sqrt(diag(varb));
   tb = srb./stdeb;

   sig2 = tmp2'tmp2/rows(dep);  @ e'e/T @
   like = -(rows(dy)/2)*ln(2*pi*sig2) - rows(dy)/2;  @ log-likelihood @
   sumlike = sumlike + like;  @ sum of log-likelihood @

@  computation of the LR parameters @

   if plag[mxi] eq 0;     @ no lagged dep vble @

      b1hat = SRb[1:k1];        @ beta1(i) @
      b2hat = SRb[k1+1:k1+k2];  @ beta2(i) @
      beta = b1hat|b2hat;

      t1hat = b1hat;     @ theta1(i) @
      t2hat = b2hat;     @ theta2(i) @

      theta   = t1hat|t2hat;
      setheta = stdeb[1:k1+k2];       @ std. err. of theta @
      ttheta  = theta./setheta;       @ t-ratios of theta @

   else;

      faihat = SRb[1];            @ fai(i) @
      b1hat = SRb[2:k1+1];        @ beta1(i) @
      b2hat = SRb[k1+2:k1+k2+1];  @ beta2(i) @
      beta = b1hat|b2hat;

      t1hat = -b1hat./faihat;     @ theta1(i) @
      t2hat = -b2hat./faihat;     @ theta2(i) @

      @  Computation of the standard errors and t-values
        for the LR parameters using the delta method @

      theta = t1hat|t2hat;
      kk = k1 + k2;

      dtheta1 = zeros(kk,1);  @ first differential  w.r.t. phi @
      h = 1;
      do until h > kk;
         dtheta1[h,1] = beta[h]/faihat^2;
         h = h + 1;
      endo;

      dtheta2 = eye(kk);     @ first differential w.r.t beta @
      dtheta2 = -(1/faihat).*dtheta2;
      dtheta  = dtheta1~dtheta2;

      vtheta  = dtheta*varb[1:kk+1,1:kk+1]*dtheta'; @ variance of theta @
      setheta = sqrt(diag(vtheta));                 @ std. err. of theta @
      ttheta  = theta./setheta;                     @ t-ratios of theta @

   endif;

   Theta1[.,mxi] = t1hat;
   Theta2[.,mxi] = t2hat;

@  computation of R2 and diagnostic statistics @

   dep = dy;
   reg = y1~x10~x20~ww;
   if plag[mxi] eq 0; epara = -1|srb;
                      nreg  = k1 + k2 + numsi[mxi];  @ # of regressors @
   else;              epara = srb;
                      nreg  = k1 + k2 + 1 + numsi[mxi];  @ # of regressors @
   endif;
   yhat = reg*epara;
   res = tmp2;

   { OLSDRES[.,mxi] } = DIAGNOS(dep,reg,res,yhat,nreg);

@  Saving the results for the estimates, std. errors and t-ratios
*  in the order of theta, fai, beta and si:
*  the dimension is (k1+k2)+1+(k1+k2)+numsi(i) @

   if plag[mxi] eq 0;

      olsb = olsb|theta|-1|srb;
      olsseb = olsseb|setheta|"NA"|stdeb;
      olstb = olstb|ttheta|"NA"|tb;

   else;

      olsb = olsb|theta|srb;
      olsseb = olsseb|setheta|stdeb;
      olstb = olstb|ttheta|tb;
   endif;

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

OLSRES = olsb~olsseb~olstb;       @ Results for estimates, std. err and t @

@ Result in the order of theta1, theta2, phi, beta1, beta2, other para
* for group 1, 2, ..., N:
* # of parameters in each group is k1+k2+1+k1+k2+numsi(i) @

OLSRES = olsres[2:rows(olsres),.];
CFname = CFname[2:rows(CFname)];

RETP(OLSRES,OLSDRES,theta1,theta2,sumlike,numsi,CFname);ENDP;

/*  ************************************************************************  */
    Proc (3) =
    PMLEBS2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,k1,k2,N,T,OLDTHTA1,OLDTHTA2,ini_op,iter);
/*  ***************************************************************************
*   Computing the Pooled ML estimatiors of the LR parameters as well as
*   the SR parameters given the LR parameters by `BACK-SUBSTITUTION'
*   See OLS2
*   ************************************************************* */

Local sumden,sumnum,sumlike,ini,fini,mxi,Tis,mxj,ztemp,
      ytemp,x1temp,x2temp,numsi,y0,y1,dy,x10,x20,ww,inpt,sig2,reg,
      tmp,tmp1,tmp2,tmp3,tmp4,tmp5,SRb,faihat,beta2hat,T1den,T1num,like,
      theta1,theta2,t2hat,oldt2,t2den,t2num,tmp0;

sumden = zeros(k1,k1);    @ denominator of Theta1 @
sumnum = zeros(k1,1);     @ numerator of Theta1 @
Theta2 = zeros(k2,n);
sumlike = 0;              @ restricted sum of log-likelihood @

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;

   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ***************************************************
*  Generation of the data used in estimation: See DGP2
*  *************************************************** */

   ytemp  = y[ini:fini];       @ y for ith group @
   x1temp = x1[ini:fini,.];    @ x1 for ith group @
   x2temp = x2[ini:fini,.];    @ x2 for ith group @
   ztemp  = z[ini:fini,.];     @ z for ith group @

   {y0, dy, y1, x10, x20, ww}
   = DGP2(ytemp,x1temp,x2temp,ztemp,plag,q1lag,q2lag,k1,k2,mxi,tis);

/* **********************************************
*  ML Estimation of the SR and the LR parameters
*  *********************************************** */

@  Estimation of phi and other SR dynamic parameters including intercept @

   if iter eq 1;
      if ini_op eq 1; oldt2 = oldthta2[.,mxi];
      else;           oldt2 = oldthta2;
      endif;
   else;              oldt2 = oldthta2[.,mxi];
   endif;

   if plag[mxi] eq 0;  faihat = -1;
                       tmp0 = y0 - x10*oldthta1 - x20*oldt2;
                       SRb = invpd(ww'ww)*ww'tmp0;
                       tmp2 = tmp0 - ww*SRb;      @ residual @
   else;
                       tmp = y1 - x10*oldthta1 - x20*oldt2;
                       reg = tmp~ww;
                       tmp1 = invpd(reg'reg);
                       SRb = tmp1*reg'dy;
                       @ SR parameters given Theta1 and Theta2 @
                       faihat = SRb[1];            @ phi @
                       tmp2 = dy - reg*SRb;      @ residual @
   endif;

   sig2 = tmp2'tmp2/rows(dy);  @ e'e/T @

   like = -(rows(dy)/2)*ln(2*pi*sig2) - rows(dy)/2;  @ log-likelihood @
   sumlike = sumlike + like;  @ sum of log-likelihood @

@  Estimation of theta2 @

   tmp3 = x20'x20 - x20'ww*invpd(ww'ww)*ww'x20;
   T2den = (faihat^2/sig2).*tmp3;
   tmp4 = dy - faihat*(y1 - x10*oldthta1);
   tmp5 = x20'tmp4 - x20'ww*invpd(ww'ww)*ww'tmp4;
   T2num = -(faihat/sig2).*tmp5;

   T2hat = invpd(T2den)*T2num;
   Theta2[.,mxi] = T2hat;

@  Estimation of theta1 @

   tmp3 = x10'x10 - x10'ww*invpd(ww'ww)*ww'x10;
   T1den = (faihat^2/sig2).*tmp3;
   tmp4 = dy - faihat*(y1 - x20*T2hat);
   tmp5 = x10'tmp4 - x10'ww*invpd(ww'ww)*ww'tmp4;
   T1num = -(faihat/sig2).*tmp5;

   sumden = sumden + T1den;   @ sum of denominator @
   sumnum = sumnum + T1num;   @ sum of numerator @

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

theta1 = invpd(sumden)*sumnum;

RETP(theta1,theta2,sumlike);ENDP;

/*  **********************************************************  */
    Proc (4) = PMLENR2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,k1,k2,N,T,
                       OLDTHTA1,OLDTHTA2,OLDPHI,ini_op,iter);
/*  ***************************************************************
*   Computing the Pooled ML estimatiors of the LR parameters and
*   the ECM parameters by `Newton-Rapson'
*   See OLS2
*   ************************************************************ */

Local sumden,sumnum,sumlike,ini,fini,mxi,Tis,mxj,ztemp,
      ytemp,x1temp,x2temp,numsi,y0,y1,dy,x10,x20,ww,inpt,sig2,reg,
      tmp,tmp1,tmp2,tmp3,tmp4,tmp5,SRb,faihat,beta2hat,T1den,T1num,like,
      theta1,theta2,t2hat,oldt2,t2den,t2num,tmp0,t1ini,t1fini,ecm,
      t2ini,t2fini,theta2fd,phifd,theta2sd,phisd,th1th2sd,th1phisd,th2phisd,
      jnk1,jnk2,jnk3,jnk4,para,oldpara,otheta2,xi,ndep,res,
      tnum,tden,tden1,tden2,tden3,phi,thta2,thi;

sumden = zeros(k1,k1);  @ k1 by k1 matrix of the second der for Theta1 block @
sumnum = zeros(k1,1);  @ k1 by 1 vectot with respect to Theta1 @

/*
theta2fd = zeros(k2*n,1);
@ k2*n by 1 vector with respect to theta2_1,...,theta_2(N) @
*/

phifd = zeros(n,1);      @ n by 1 vector with respect to phi_1,..,phi_N @

theta2sd = zeros(k2*n,k2*n);
@ k2*n by k2*n matrix of Hessian for theta2_1,...,theta_2(N) @

phisd = zeros(n,n);      @ n by n matrix of the second der for phi_1,..,phi_N @

th1th2sd = zeros(k1,k2*n);
@ k1 by k2*n matrix with respect to theta1 and then theta2 @

th1phisd = zeros(k1,n);  @ k1 by n matrix with respect to theta1 and then phi @

th2phisd = zeros(k2*n,n);
@ k2*n by n matrix with respect to theta2 and then phi @

Theta2 = zeros(k2,n);
sumlike = 0;              @ restricted sum of log-likelihood @

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
t1ini = 1; t1fini = k1;
t2ini = 1; t2fini = k2;
mxi = 1;
do until mxi > N;

   faihat = oldphi[mxi];

   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ***************************************************
*  Generation of the data used in estimation: See DGP2
*  *************************************************** */

   ytemp  = y[ini:fini];       @ y for ith group @
   x1temp = x1[ini:fini,.];    @ x1 for ith group @
   x2temp = x2[ini:fini,.];    @ x2 for ith group @
   ztemp  = z[ini:fini,.];     @ z for ith group @

   {y0, dy, y1, x10, x20, ww}
   = DGP2(ytemp,x1temp,x2temp,ztemp,plag,q1lag,q2lag,k1,k2,mxi,tis);

/* **********************************************
*  ML Estimation of the SR and the LR parameters
*  *********************************************** */

@  Construction of the residual, sigma, error correction term
   given theta1(j-1), theta2(j-1) and phi(j-1) @

   if iter eq 1;
      if ini_op eq 1; oldt2 = oldthta2[.,mxi];
      else;           oldt2 = oldthta2;
      endif;
   else;              oldt2 = oldthta2[.,mxi];
   endif;

   if mxi eq 1;      OTHETA2 = oldt2;
   else;             OTHETA2 = Otheta2|oldt2;
   endif;

   xi = y1 - x10*oldthta1 - x20*oldt2;   @ EC term @
   ndep = dy - xi*faihat;                @ dy - phi*ecm @
   res = ndep - ww*invpd(ww'ww)*ww'ndep; @ residual @
   sig2 = res'res/rows(dy);              @ sigma^2 @

   like = -(rows(dy)/2)*ln(2*pi*sig2) - rows(dy)/2;  @ log-likelihood @
   sumlike = sumlike + like;  @ sum of log-likelihood @

@  Estimation of theta1 theta2 and phi @
@  Construction of the gradient vector @

   jnk1 = x10'res - x10'ww*invpd(ww'ww)*ww'res;   @ X1(i)'Hi*e(i) @
   jnk1 = -(faihat/sig2).*jnk1;
   sumnum = sumnum + jnk1;   @ for theta1 @

   jnk1 = x20'res - x20'ww*invpd(ww'ww)*ww'res;   @ X2(i)'Hi*e(i) @
   jnk2 = -(faihat/sig2).*jnk1;
   if mxi eq 1; theta2fd = jnk2;
   else;        theta2fd = theta2fd|jnk2;
   endif;

   jnk1 = xi'res - xi'ww*invpd(ww'ww)*ww'res;   @ xi(i)'Hi*e(i) @
   phifd[mxi] = (1/sig2)*jnk1;

@  Construction of the Hesssian Matrix @

   jnk3  = x10'x10 - x10'ww*invpd(ww'ww)*ww'x10;   @ X1(i)'Hi*X1(i) @
   jnk3 = (faihat^2/sig2).*jnk3;
   sumden = sumden + jnk3;   @ for theta1: k1 by k1 @

   jnk4  = x20'x20 - x20'ww*invpd(ww'ww)*ww'x20;   @ X2(i)'Hi*X2(i) @
   jnk4 = (faihat^2/sig2).*jnk4;
   theta2sd[t2ini:t2fini,t2ini:t2fini] = jnk4;

   jnk4 = xi'xi - xi'ww*invpd(ww'ww)*ww'xi;  @ xi(i)'Hi*xi(i) @
   phisd[mxi,mxi] = (1/sig2)*jnk4;

   jnk4 = x10'x20 - x10'ww*invpd(ww'ww)*ww'x20;  @ X1(i)'Hi*X2(i) @
   th1th2sd[.,t2ini:t2fini] = (faihat^2/sig2).*jnk4;

   jnk4 = x10'xi - x10'ww*invpd(ww'ww)*ww'xi;  @ X1(i)'Hi*xi(i): k1 by 1 @
   th1phisd[.,mxi] = -(faihat/sig2).*jnk4;

   jnk4 = x20'xi - x20'ww*invpd(ww'ww)*ww'xi;  @ X2(i)'Hi*xi(i): k2 by 1 @
   th2phisd[t2ini:t2fini,mxi] = -(faihat/sig2).*jnk4;

   ini = ini + Tis;
   mxi = mxi + 1;

   t2ini = t2ini + k2; t2fini = t2fini + k2;

endo;

TNUM  = sumnum|theta2fd|phifd;          @ k1+k2*n+n by 1 @
TDEN1 = sumden~th1th2sd~th1phisd;       @ k1 by (k1+k2*n+n) @
TDEN2 = th1th2sd'~theta2sd~th2phisd;    @ k2 by (k1+k2*n+n) @
TDEN3 = th1phisd'~th2phisd'~phisd;      @ k2 by (k1+k2*n+n) @
TDEN = TDEN1|Tden2|Tden3;        @ (k1+k2*n+n) by (k1+k2*n+n) @

/* output on; format 5,1; "CHK" Tden;""; */

oldpara = oldthta1|OTHETA2|oldphi;
para = oldpara + invpd(TDEN)*Tnum;
theta1 = para[1:k1];
thta2 = para[k1+1:k1+k2*n];
phi = para[k1+k2*n+1:rows(para)];

@ convert to matrix: k2 by n @

theta2 = zeros(k2,n);
ini = 1;
fini = k2;
thi = 1;
do until thi > n;

   theta2[.,thi] = thta2[ini:fini];
   thi = thi + 1;
   ini = ini + k2;
   fini = fini + k2;
endo;

RETP(theta1,theta2,phi,sumlike);ENDP;

/*  ********************************************************  */
    Proc (1) = INIPHI2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,k1,k2,N,T,
                       OLDTHTA1,OLDTHTA2,ini_op);
/*  ************************************************************
*   Computing the initial estimates of the ECM parameter
*   which will be used in `NEWTON-RAPHSON'
*   ****************************************** */

Local phiini,ini,fini,mxi,Tis,mxj,ztemp,xi,den,num,
      ytemp,x1temp,x2temp,numsi,y0,y1,dy,x10,x20,ww,inpt,sig2,reg,
      theta1,theta2,t2hat,oldt2,t2den,t2num,tmp0,t1ini,t1fini,ecm,
      t2ini,t2fini;

phiini = zeros(n,1);      @ n by 1 vector of initial estimates of phi @

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;

   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ***************************************************
*  Generation of the data used in estimation: See DGP2
*  *************************************************** */

   ytemp  = y[ini:fini];       @ y for ith group @
   x1temp = x1[ini:fini,.];    @ x1 for ith group @
   x2temp = x2[ini:fini,.];    @ x2 for ith group @
   ztemp  = z[ini:fini,.];     @ z for ith group @

   {y0, dy, y1, x10, x20, ww}
   = DGP2(ytemp,x1temp,x2temp,ztemp,plag,q1lag,q2lag,k1,k2,mxi,tis);

@  The error correction term given theta1(j-1) and theta2(j-1) @

   if ini_op eq 1; oldt2 = oldthta2[.,mxi];
   else;           oldt2 = oldthta2;
   endif;

   xi = y1 - x10*oldthta1 - x20*oldt2;     @ ec term @

@  Estimates of phi @

   den = xi'xi - xi'ww*invpd(ww'ww)*ww'xi;   @ denominator @
   num = xi'dy - xi'ww*invpd(ww'ww)*ww'dy;   @ numerator @

   if plag[mxi] eq 0; phiini[mxi] = -1;
   else;              phiini[mxi] = num/den;
   endif;

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

RETP(phiini);ENDP;

/*  *******************************************************************  */
    Proc (4)
    = PMLESR2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,k1,k2,N,T,NUMSI,THETA1,THETA2);
/*  *******************************************************************
*   See comments in the procedure PMLE2
*   Final estimation results using the converged theta1
*   *************************************************** */

Local sumlike,ini,fini,mxi,Tis,p_cnt,q1_cnt,q2_cnt,mxj,sumsi,
      g12ini,g12fini,g13ini,g13fini,g14ini,g14fini,gg1,gg2,
      g2ini,g2fini,rg23ini,rg23fini,cg23ini,cg23fini,
      rg24ini,rg24fini,cg24ini,cg24fini,g4ini,g4fini,dbeta1,dbeta2,
      g3ini,g3fini,rg34fini,rg34ini,cg34ini,cg34fini,
      g11,g12,g13,g14,g22,g23,g24,g33,g34,g44,a2tmp,a3tmp,a4tmp,tmp0,
      g1tmp,g12tmp,g13tmp,g14tmp,g2tmp,g23tmp,g24tmp,g3tmp,g34tmp,g4tmp,
      ytemp,x1temp,x2temp,y0,y1,dy,x10,x20,ww,inpt,sig2,reg,ztemp,
      tmp,tmp1,tmp2,tmp3,tmp4,tmp5,SRb,faihat,beta2hat,T1den,T1num,like,
      allpara,se_all,t_all,gg,si,ecm1,oldt2,t2hat,t2den,t2num,
      beta1,beta2,beta,btmp,sebtmp,tbtmp,ini12,ini13,ini22,ini23,ini33,
      fini12,fini13,fini22,fini23,fini33,kk,dbeta,theta,h,ggtmp,
      vtmp11,vtmp12,vtmp13,vtmp22,vtmp23,vtmp33,vtmp1,vtmp2,vtmp3,vtmp,
      vbeta,sebeta,tbeta,seall,tall,pmleres,PMLECOV,
      sebeta1,sebeta2,tbeta2,tbeta1,chk2,dtmp,pmledres,nreg,
      RSQ,RBARSQ,LL,AIC,SBC,HQ,chiSC,chiFF,chiNO,chiHE,sige,dep,yhat,res,
      epara,ecm,na,nfai0,nn,a3ntmp,hh,atmp,htmp,hflag,a3se,a3t,a3temp,hini;


PMLEDRES = zeros(11,n);

beta1 = zeros(k1,n);  @ SR parameters @
beta2 = zeros(k2,n);

a2tmp = 0;  @ initialization for theta2, phi and si @
a3tmp = 0;
a4tmp = 0;

nfai0 = 0;  @ number of phi(i) = -1 @

sumsi = sumc(numsi); @ sum of other parameters for all groups @
                     @ note that numsi could be different!! @

@ Information matrix G @

@ Theta1 block @

G11 = zeros(k1,k1);
G12 = zeros(k1,k2*n);
G13 = zeros(k1,n);
G14 = zeros(k1,sumsi);

@ Theta2 block @

G22 = zeros(k2*n,k2*n);
G23 = zeros(k2*n,n);
G24 = zeros(k2*n,sumsi);

@ Fai block @

G33 = zeros(n,n);
G34 = zeros(n,sumsi);

@ Si block @

G44 = zeros(sumsi,sumsi);

@ initial and end positions @

g12fini = 0;     g12ini = 1;
g13fini = 0;     g13ini = 1;
g14fini = 0;     g14ini = 1;

g2fini = 0;      g2ini = 1;
rg23fini = 0;    rg23ini = 1;
cg23fini = 0;    cg23ini = 1;
rg24fini = 0;    rg24ini = 1;
cg24fini = 0;    cg24ini = 1;

g3fini = 0;      g3ini = 1;
rg34fini = 0;    rg34ini = 1;
cg34fini = 0;    cg34ini = 1;

g4fini = 0;      g4ini = 1;

sumlike = 0;              @ restricted sum of log-likelihood @

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;

   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ***************************************************
*  Generation of the data used in estimation: See DGP2
*  *************************************************** */

   ytemp  = y[ini:fini];      @ y for ith group @
   x1temp = x1[ini:fini,.];   @ x1 for ith group @
   x2temp = x2[ini:fini,.];   @ x2 for ith group @
   ztemp  = z[ini:fini,.];    @ z for ith group @

   {y0, dy, y1, x10, x20, ww}
   = DGP2(ytemp,x1temp,x2temp,ztemp,plag,q1lag,q2lag,k1,k2,mxi,tis);

/* *********************************************
*  ML Estimation of the SR and the LR parameters
*  ********************************************* */

@  Estimation of phi and fai2 @

   T2hat = theta2[.,mxi];

   if plag[mxi] eq 0;  faihat = -1;
                       nfai0 = nfai0 + 1;
                       tmp0 = y0 - x10*theta1 - x20*t2hat;
                       si = invpd(ww'ww)*ww'tmp0;
                       tmp2 = tmp0 - ww*si;      @ residual @
                       nreg = k1 + k2 + numsi[mxi];  @ # of regressors @

   else;               tmp = y1 - x10*theta1 - x20*t2hat;
                       reg = tmp~ww;
                       tmp1 = invpd(reg'reg);
                       SRb = tmp1*reg'dy;
                       @ SR parameters given Theta1 and Theta2 @
                       faihat = SRb[1];            @ phi @
                       si = SRb[2:rows(srb)];
                       tmp2 = dy - reg*SRb;      @ residual @
                       nreg = k1 + k2 + 1 + numsi[mxi];  @ # of regressors @
   endif;

   sig2 = tmp2'tmp2/rows(dy);  @ e'e/T @

   like = -(rows(dy)/2)*ln(2*pi*sig2) - rows(dy)/2;  @ log-likelihood @
   sumlike = sumlike + like;  @ sum of log-likelihood @

@  computation of the SR parameters: beta @

   beta1[.,mxi] = -faihat.*theta1;        @ beta1(i) @
   beta2[.,mxi] = -faihat.*T2hat;         @ beta2(i) @

@  Generation of the final ecm term: y(-1) - x1*theta1 - x2*theta2 @

   ecm = y1 - x10*theta1 - x20*t2hat;
   reg = ecm~ww;

@  computation of R2 and diagnostic statistics @

   dep = dy;
   reg = y1~x10~x20~ww;
   epara = faihat|beta1[.,mxi]|beta2[.,mxi]|si;
   yhat = reg*epara;
   res = dep - yhat;

   { PMLEDRES[.,mxi] } = DIAGNOS(dep,reg,res,yhat,nreg);

@  Estimation of the Covariance matrix for all parameters @

   a2tmp = a2tmp|t2hat;     @ stacking theta2(i) @
   a3tmp = a3tmp|faihat;    @ stacking fai(i) @
   a4tmp = a4tmp|si;        @ stacking si(i) @

@ Theta1 block of information matrix @

   G1tmp = (faihat^2/sig2).*(x10'x10);
   G11 = G11 + G1tmp;

   G12tmp = (faihat^2/sig2).*(x10'x20);
   g12fini = g12fini + k2;
   G12[.,g12ini:g12fini] = G12tmp;
   g12ini = g12ini + k2;

   G13tmp = -(faihat/sig2).*(x10'ecm);
   g13fini = g13fini + 1;
   G13[.,g13ini:g13fini] = G13tmp;
   g13ini = g13ini + 1;

   G14tmp = -(faihat/sig2).*(x10'ww);
   g14fini = g14fini + numsi[mxi];
   output on; format /rd 6,3;
   G14[.,g14ini:g14fini] = G14tmp;
   g14ini = g14ini + numsi[mxi];

@ Theta2 block of information matrix @

   G2tmp = (faihat^2/sig2).*(x20'x20);
   g2fini = g2fini + k2;
   G22[g2ini:g2fini,g2ini:g2fini] = G2tmp;
   g2ini = g2ini + k2;

   G23tmp = -(faihat/sig2).*(x20'ecm);
   rg23fini = rg23fini + k2;
   cg23fini = cg23fini + 1;
   G23[rg23ini:rg23fini,cg23ini:cg23fini] = G23tmp;
   rg23ini = rg23ini + k2;
   cg23ini = cg23ini + 1;

   G24tmp = -(faihat/sig2).*(x20'ww);
   rg24fini = rg24fini + k2;
   cg24fini = cg24fini + numsi[mxi];
   G24[rg24ini:rg24fini,cg24ini:cg24fini] = G24tmp;
   rg24ini = rg24ini + k2;
   cg24ini = cg24ini + numsi[mxi];

@ Fai block of information matrix @

   G3tmp = (1/sig2).*(ecm'ecm);
   g3fini = g3fini + 1;
   if plag[mxi] eq 0; G33[g3ini:g3fini,g3ini:g3fini] = 0; @ no varaince @
   else;              G33[g3ini:g3fini,g3ini:g3fini] = G3tmp;
   endif;
   g3ini = g3ini + 1;

   G34tmp = (1/sig2).*(ecm'ww);
   rg34fini = rg34fini + 1;
   cg34fini = cg34fini + numsi[mxi];
   G34[rg34ini:rg34fini,cg34ini:cg34fini] = G34tmp;
   rg34ini = rg34ini + 1;
   cg34ini = cg34ini + numsi[mxi];

@ Si block of information matrix @

   G4tmp = (1/sig2).*(ww'ww);
   g4fini = g4fini + numsi[mxi];
   G44[g4ini:g4fini,g4ini:g4fini] = G4tmp;
   g4ini = g4ini + numsi[mxi];

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

@ delete unnecessary initial values @

a2tmp = a2tmp[2:rows(a2tmp)];  @ k2*n by 1 @
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

   g1tmp = G11~G12~G13~G14;
   g2tmp = G12'~G22~G23~G24;
   g3tmp = G13'~G23'~G33~G34;
   g4tmp = G14'~G24'~G34'~G44;

   GG = g1tmp|g2tmp|g3tmp|g4tmp;

@  Deal with the case where some of phi(i) = -1 @

   if nfai0 eq 0;
   else;

      hflag = 0;
GREPEAT:
      hh = k1 + k2*n + 1;
      do until hh > k1+ k2*n + n;

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
PMLECOV = GGtmp[1:k1,1:k1];

@ standard errors and t-values @

if nn eq 0; allpara = theta1|a2tmp|a4tmp;
else;       allpara = theta1|a2tmp|a3ntmp[.,1]|a4tmp;
endif;

se_all = sqrt(diag(GGtmp));
t_all = allpara./se_all;

@ including the case of no standard errors for phi(i) = -1 @

a3se = zeros(n,1);
a3t  = zeros(n,1);
if nn eq 0;
else;
   hini = k1 + k2*n + 1;
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

   sebeta1 = zeros(k1,n);
   tbeta1 = zeros(k1,n);
   sebeta2 = zeros(k2,n);
   tbeta2 = zeros(k2,n);

   ini12 = k1+1;         fini12 = k1;
   ini13 = k1+k2*n+1;    fini13 = k1+k2*n;
   ini22 = k1+1;         fini22 = k1;
   ini23 = k1+k2*n+1;    fini23 = k1+k2*n;
   ini33 = k1+k2*n+1;    fini33 = k1+k2*n;

   ini = 1;
   fini = 0;
   kk = k1 + k2;
   mxi = 1;
   do until mxi > n;

      beta = beta1[.,mxi]|beta2[.,mxi];
      fini = fini + k2;
      theta = theta1|a2tmp[ini:fini];
      ini = ini + k2;
      faihat = a3tmp[mxi];

      if plag[mxi] eq 0;   @ no lagged dep vble @

         vtmp11 = GGtmp[1:k1,1:k1];
         fini12 = fini12 + k2;
         vtmp12 = GGtmp[1:k1,ini12:fini12];
         ini12 = ini12 + k2;
         fini22 = fini22 + k2;
         vtmp22 = GGtmp[ini22:fini22,ini22:fini22];
         ini22 = ini22 + k2;

         vtmp1 = vtmp11~vtmp12;
         vtmp2 = vtmp12'~vtmp22;

         vbeta = vtmp1|vtmp2;    @ Variance for theta1 and theta2 as well
                                as beta1 and beta2 @

      else;

         vtmp11 = GGtmp[1:k1,1:k1];             @ k1 by k1 @
         fini12 = fini12 + k2;
         vtmp12 = GGtmp[1:k1,ini12:fini12];     @ k1 by k2 @
         ini12 = ini12 + k2;
         fini13 = fini13 + 1;
         vtmp13 = GGtmp[1:k1,ini13:fini13];     @ k1 by 1 @
         ini13 = ini13 + 1;

         fini22 = fini22 + k2;
         vtmp22 = GGtmp[ini22:fini22,ini22:fini22];   @ k2 by k2 @
         fini23 = fini23 + 1;
         vtmp23 = GGtmp[ini22:fini22,ini23:fini23];   @ k2 by 1 @
         ini23  = ini23 + 1;
         ini22  = ini22 + k2;
         fini33 = fini33 + 1;
         vtmp33 = GGtmp[ini33:fini33,ini33:fini33];   @ 1 by 1 @
         ini33  = ini33 + 1;

         vtmp1 = vtmp11~vtmp12~vtmp13;
         vtmp2 = vtmp12'~vtmp22~vtmp23;
         vtmp3 = vtmp13'~vtmp23'~vtmp33;
         vtmp = vtmp1|vtmp2|vtmp3;    @ variance for theta and phi @

         dbeta1 = eye(kk);     @ first differential w.r.t beta @
         dbeta1 = -faihat.*dbeta1;

         dbeta2 = zeros(kk,1);  @ first differential  w.r.t. phi @
         h = 1;
         do until h > kk;

            dbeta2[h,1] = -theta[h];
            h = h + 1;
         endo;
         dbeta = dbeta1~dbeta2;

         vbeta = dbeta*vtmp*dbeta';   @ variance of beta @

      endif;

      sebeta = sqrt(diag(vbeta));  @ std. err. of beta @
      tbeta = beta./sebeta;        @ t-ratios of beta @

@  Saving the results for the estimates, std. errors and t-ratios @

      sebeta1[.,mxi] = sebeta[1:k1];
      sebeta2[.,mxi] = sebeta[k1+1:k1+k2];
      tbeta1[.,mxi] = tbeta[1:k1];
      tbeta2[.,mxi] = tbeta[k1+1:k1+k2];

   mxi = mxi + 1;
endo;

beta1 = vec(beta1);
sebeta1 = vec(sebeta1);
tbeta1 = vec(tbeta1);

beta2 = vec(beta2);
sebeta2 = vec(sebeta2);
tbeta2 = vec(tbeta2);

@ All parameters in the order of theta1, theta2(1),...theta2(N),
* phi(1), ..., phi(N), and si(1), ..., si(N):
* The dimension is (k1+k2*n)+n+(k1+k2*n)+sumnumsi-vector @

allpara = theta1|a2tmp|a3tmp|beta1|beta2|a4tmp;
fini = k1 + k2*n;
seall = se_all[1:fini]|a3se|sebeta1|sebeta2|se_all[fini+nn+1:rows(se_all)];
tall = t_all[1:fini]|a3t|tbeta1|tbeta2|t_all[fini+nn+1:rows(t_all)];

PMLERES = allpara~seall~tall;

RETP(PMLERES,PMLEDRES,sumlike,PMLECOV);ENDP;

/* ************************************************************************** */
   PROC (6) = DGP2(ytemp,x1temp,x2temp,ztemp,plag,q1lag,q2lag,k1,k2,mxi,tis);
/* **************************************************************************
*  Generation of the data used in estimation
*
*  y0: current value of dep vble
*  dy: first difference of dep vble
*  y1: lagged value of dep vble
*  x10: current value of indep vbles, x1
*  x20: current value of indep vbles, x2
*  ww: dy's, dx1's, dx2's and intercept
*  ************************************* */

Local y1temp,dytemp,x11temp,dx1temp,x21temp,dx2temp,
      dy1temp,dy2temp,dy3temp,dystemp,wflag,
      dx10temp,dx11temp,dx12temp,dx13temp,dx1tmp,dx1flag,
      dx20temp,dx21temp,dx22temp,dx23temp,dx2tmp,dx2flag,
      WW,ltemp,mxlag,sini,sig2,mxj,y0,y1,dy,x10,x20;

 @ Generation y, dy, y(-1), x1, x2 and ww @

   y1temp  = lag(ytemp);        @ y(-1) @
   dytemp  = ytemp - y1temp;    @ dy @

   x11temp  = lag(x1temp);      @ x1(-1) @
   dx1temp  = x1temp - x11temp; @ dx1 @

   x21temp  = lag(x2temp);      @ x2(-1) @
   dx2temp  = x2temp - x21temp; @ dx2 @

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

 @ Generation of current and lagged dx1's @

   dx1tmp = zeros(Tis,1);
   mxj = 1;
   do until mxj > k1;

      if     Q1LAG[mxi,mxj] eq 1;  dx10temp = dx1temp[.,mxj];
                                   dx1tmp = dx1tmp~dx10temp;
      elseif Q1LAG[mxi,mxj] eq 2;  dx10temp = dx1temp[.,mxj];
                                   dx11temp = lag(dx10temp);
                                   dx1tmp = dx1tmp~dx10temp~dx11temp;
      elseif Q1LAG[mxi,mxj] eq 3;  dx10temp = dx1temp[.,mxj];
                                   dx11temp = lag(dx10temp);
                                   dx12temp = lag(dx11temp);
                          dx1tmp = dx1tmp~dx10temp~dx11temp~dx12temp;
      elseif Q1LAG[mxi,mxj] eq 4;  dx10temp = dx1temp[.,mxj];
                                   dx11temp = lag(dx10temp);
                                   dx12temp = lag(dx11temp);
                                   dx13temp = lag(dx12temp);
                  dx1tmp = dx1tmp~dx10temp~dx11temp~dx12temp~dx13temp;
      endif;

      mxj = mxj + 1;
   endo;

   if cols(dx1tmp) eq 1;  dx1flag = 0;       @ No dx1's in the regression @
   else; dx1tmp = dx1tmp[.,2:cols(dx1tmp)];
         @ column dimension is q_1(i) + q_2(i) + ... + q_k(i) @
         dx1flag = 1;
   endif;

   if dx1flag eq 0;
   else;   if wflag eq 1;  WW = WW~dx1tmp;
           else;           WW = dx1tmp;
           endif;
           wflag = 1;  @ W's not empty anymore @
   endif;

 @ Generation of lagged dx2's @

   dx2tmp = zeros(Tis,1);
   mxj = 1;
   do until mxj > k2;

      if     Q2LAG[mxi,mxj] eq 1;  dx20temp = dx2temp[.,mxj];
                                   dx2tmp = dx2tmp~dx20temp;
      elseif Q2LAG[mxi,mxj] eq 2;  dx20temp = dx2temp[.,mxj];
                                   dx21temp = lag(dx20temp);
                                   dx2tmp = dx2tmp~dx20temp~dx21temp;
      elseif Q2LAG[mxi,mxj] eq 3;  dx20temp = dx2temp[.,mxj];
                                   dx21temp = lag(dx20temp);
                                   dx22temp = lag(dx21temp);
                          dx2tmp = dx2tmp~dx20temp~dx21temp~dx22temp;
      elseif Q2LAG[mxi,mxj] eq 4;  dx20temp = dx2temp[.,mxj];
                                   dx21temp = lag(dx20temp);
                                   dx22temp = lag(dx21temp);
                                   dx23temp = lag(dx22temp);
                 dx2tmp = dx2tmp~dx20temp~dx21temp~dx22temp~dx23temp;
      endif;

      mxj = mxj + 1;
   endo;

   if cols(dx2tmp) eq 1; dx2flag = 0;       @ No dx2's in the regression @
   else; dx2tmp = dx2tmp[.,2:cols(dx2tmp)];
         @ column dimension is q_1(i) + q_2(i) + ... + q_k(i) @
         dx2flag = 1;
   endif;

   if wflag eq 0;  if dx2flag eq 0;
                   else;  WW = dx2tmp; wflag = 1;
                   endif;
   else;           if dx2flag eq 0;
                   else;  WW = WW~dx2tmp;
                   endif;
   endif;

   if wflag eq 0; WW = ztemp;
   else;          WW = WW~ztemp;
   endif;

 @ sample adjustment @

   ltemp = plag[mxi]|q1lag[mxi,.]'|q2lag[mxi,.]';
   mxlag = maxc(ltemp);      @ maximum lag order for each group @

   if mxlag eq 0; sini = 2;
   else;          sini = mxlag + 1;
   endif;

   y0 = ytemp[sini:Tis];        @ y @
   y1 = y1temp[sini:Tis];       @ y(-1) @
   dy = dytemp[sini:Tis];       @ dy @

   x10 = x1temp[sini:Tis,.];    @ x1 @
   x20 = x2temp[sini:Tis,.];    @ x2 @
   ww = ww[sini:Tis,.];         @ other variables including intercept @

RETP(y0,dy,y1,x10,x20,ww); ENDP;

/* ************************************************************************* */
   PROC (7) = NDGP2(ytemp,x1temp,x2temp,ztemp,plag,q1lag,q2lag,k1,k2,mxi,tis,
                   yname,x1name,x2name,Zname);
/* *******************************************************************
*  Generation of the data used in estimation: See DGP2
*  *************************************************** */

Local y1temp,dytemp,x11temp,dx1temp,x21temp,dx2temp,
      dy1temp,dy2temp,dy3temp,dystemp,wflag,
      dx10temp,dx11temp,dx12temp,dx13temp,dx1tmp,dx1flag,
      dx20temp,dx21temp,dx22temp,dx23temp,dx2tmp,dx2flag,
      inpt,ltemp,mxlag,sini,sig2,mxj,dx1name,dx2name,dyname,
      y0,y1,dy,x10,x20,ww,coname,srname1,srname2,srname,lrname1,lrname2,
      dxtmp0,dxtmp1,dztmp0,dztmp1,s;

 @ Generation y, dy, y(-1), x1, x2 and ww @

   y1temp = lag(ytemp);        @ y(-1) @
   dytemp = ytemp - y1temp;    @ dy @

   x11temp = lag(x1temp);      @ x1(-1) @
   dx1temp = x1temp - x11temp; @ dx1 @

   x21temp = lag(x2temp);      @ x2(-1) @
   dx2temp = x2temp - x21temp; @ dx2 @

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

   dx1tmp = zeros(Tis,1);
   mxj = 1;
   do until mxj > k1;

      if     Q1LAG[mxi,mxj] eq 1;
             dx10temp = dx1temp[.,mxj];
             dx1tmp = dx1tmp~dx10temp;

             dx1name = 0 $+ "d" $+ x1name[mxj];
             SRname = SRname|dx1name;

      elseif Q1LAG[mxi,mxj] eq 2;
             dx10temp = dx1temp[.,mxj];
             dx11temp = lag(dx10temp);
             dx1tmp = dx1tmp~dx10temp~dx11temp;

             dxtmp0 = 0 $+ "d" $+ x1name[mxj];
             s = seqa(1,1,1);
             dxtmp1 = 0 $+ "d" $+ x1name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dx1name = dxtmp0|dxtmp1;
             SRname = SRname|dx1name;

      elseif Q1LAG[mxi,mxj] eq 3;
             dx10temp = dx1temp[.,mxj];
             dx11temp = lag(dx10temp);
             dx12temp = lag(dx11temp);
             dx1tmp = dx1tmp~dx10temp~dx11temp~dx12temp;

             dxtmp0 = 0 $+ "d" $+ x1name[mxj];
             s = seqa(1,1,2);
             dxtmp1 = 0 $+ "d" $+ x1name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dx1name = dxtmp0|dxtmp1;
             SRname = SRname|dx1name;

      elseif Q1LAG[mxi,mxj] eq 4;
             dx10temp = dx1temp[.,mxj];
             dx11temp = lag(dx10temp);
             dx12temp = lag(dx11temp);
             dx13temp = lag(dx12temp);
             dx1tmp = dx1tmp~dx10temp~dx11temp~dx12temp~dx13temp;

             dxtmp0 = 0 $+ "d" $+ x1name[mxj];
             s = seqa(1,1,3);
             dxtmp1 = 0 $+ "d" $+ x1name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dx1name = dxtmp0|dxtmp1;
             SRname = SRname|dx1name;

      endif;

      mxj = mxj + 1;
   endo;

   if cols(dx1tmp) eq 1;  dx1flag = 0;       @ No dx1's in the regression @
   else; dx1tmp = dx1tmp[.,2:cols(dx1tmp)];
         @ column dimension is q_1(i) + q_2(i) + ... + q_k(i) @
         dx1flag = 1;
   endif;

   if dx1flag eq 0;
   else;   if wflag eq 1;  WW = WW~dx1tmp;
           else;           WW = dx1tmp;
           endif;
           wflag = 1;  @ W's not empty anymore @
   endif;

 @ Generation of current and lagged dx2's @

   dx2tmp = zeros(Tis,1);
   mxj = 1;
   do until mxj > k2;

      if     Q2LAG[mxi,mxj] eq 1;
             dx20temp = dx2temp[.,mxj];
             dx2tmp = dx2tmp~dx20temp;

             dx2name = 0 $+ "d" $+ x2name[mxj];
             SRname = SRname|dx2name;

      elseif Q2LAG[mxi,mxj] eq 2;
             dx20temp = dx2temp[.,mxj];
             dx21temp = lag(dx20temp);
             dx2tmp = dx2tmp~dx20temp~dx21temp;

             dztmp0 = 0 $+ "d" $+ x2name[mxj];
             s = seqa(1,1,1);
             dztmp1 = 0 $+ "d" $+ x2name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dx2name = dztmp0|dztmp1;
             SRname = SRname|dx2name;

      elseif Q2LAG[mxi,mxj] eq 3;
             dx20temp = dx2temp[.,mxj];
             dx21temp = lag(dx20temp);
             dx22temp = lag(dx21temp);
             dx2tmp = dx2tmp~dx20temp~dx21temp~dx22temp;

             dztmp0 = 0 $+ "d" $+ x2name[mxj];
             s = seqa(1,1,2);
             dztmp1 = 0 $+ "d" $+ x2name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dx2name = dztmp0|dztmp1;
             SRname = SRname|dx2name;

      elseif Q2LAG[mxi,mxj] eq 4;
             dx20temp = dx2temp[.,mxj];
             dx21temp = lag(dx20temp);
             dx22temp = lag(dx21temp);
             dx23temp = lag(dx22temp);
             dx2tmp = dx2tmp~dx20temp~dx21temp~dx22temp~dx23temp;

             dztmp0 = 0 $+ "d" $+ x2name[mxj];
             s = seqa(1,1,3);
             dztmp1 = 0 $+ "d" $+ x2name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dx2name = dztmp0|dztmp1;
             SRname = SRname|dx2name;

      endif;

      mxj = mxj + 1;
   endo;

   if cols(dx2tmp) eq 1; dx2flag = 0;       @ No dx2's in the regression @
   else; dx2tmp = dx2tmp[.,2:cols(dx2tmp)];
         @ column dimension is q_1(i) + q_2(i) + ... + q_k(i) @
         dx2flag = 1;
   endif;

   if wflag eq 0;  if dx2flag eq 0;
                   else;  WW = dx2tmp; wflag = 1;
                   endif;
   else;           if dx2flag eq 0;
                   else;  WW = WW~dx2tmp;
                   endif;
   endif;

   if wflag eq 0; WW = ztemp;
   else;          WW = WW~ztemp;
   endif;
   SRname = SRname|Zname;
   SRname = SRname[2:rows(SRname)];

 @ sample adjustment @

   ltemp = plag[mxi]|q1lag[mxi,.]'|q2lag[mxi,.]';
   mxlag = maxc(ltemp);      @ maximum lag order for each group @

   if mxlag eq 0; sini = 2;
   else;          sini = mxlag + 1;
   endif;

   y0 = ytemp[sini:Tis];        @ y @
   y1 = y1temp[sini:Tis];       @ y(-1) @
   dy = dytemp[sini:Tis];       @ dy @

   x10 = x1temp[sini:Tis,.];    @ x1 @
   x20 = x2temp[sini:Tis,.];    @ x2 @
   ww = ww[sini:Tis,.];         @ other variables including intercept @

COname = X1name|X2name|"Phi"|X1name|X2name|SRname;

RETP(y0,dy,y1,x10,x20,ww,Coname); ENDP;

