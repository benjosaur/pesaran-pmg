/*  ***************************************************
*   Sub-procedures computing the static and the dynamic
*   fixed effects estimates and the standard errors
*   *************************************************** */

/*  *************************** */
    PROC(6) = SFE(Y,X,Z,n,T,k);
/*  ****************************************************
*   Static fixed effects estimation of the LR parameters
*
*   Model:
*
*   dy(i) = -(y1(i) - theta*x(i)) + a(i) + e(i)  (1)
*
*   estimate (1) by FE, and obtain the diagnostics by
*
*   dQY = -QY1 + beta*QX + E
*
*   Note here that x may contain fixed regressors excluding intercept.
*   ****************************************************************** */

Local sfedres,mxi,tis,y0,x0,inpt,tmpy,tmpx,fden,fnum,sumden,sumnum,ini,fini,
      qy,qx,qtemp,qres,qsig2,qvar,ftheta,qse,qt,qint,meandep,meanreg,dep,reg,
      tini,tfini,tmreg,tmdep,tdep,treg,ei1,Rvar,RQvar,RQse,RQt,sum,res,nreg,
      qdep,qreg,qu,temp,pmatrix,qindse,qindt,qind,sferes,yhat,dtmp,ei,difft,
      RSQ,RBARSQ,LL,AIC,SBC,HQ,chiSC,chiFF,chiNO,chiHE,sige,sfeind,rsferes,
      dy,y1,qy1,qdy,tmpy1,tmpdy,qecm,meany,meanx,mnx,epara,z0,ntmp,qobs,
      tchk;

@ Initialization @

ntmp = k + cols(z)-1;

sumden = zeros(ntmp,ntmp);
sumnum = zeros(ntmp,1);

meanY = zeros(n,1);           @ mean of Y(1), ..., Y(n) @
meanX = zeros(n,ntmp);           @ mean of X(1), ..., X(n) @

ini  = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @

mxi = 1;
do until mxi > N;

   Tis  = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ******************************************
*  Generation of the data used in estimation
*
*  y0: current value of dep vble
*  x0: current value of indep vbles, x
*  z0: current value of fixed regressors
*  ****************************************** */

   y0 = y[ini:fini];       @ y for ith group @
   y1 = lag(y0);           @ y(-1) for ith group @
   dy = y0 - y1;           @ dy for ith group @
   x0 = x[ini:fini,.];     @ x for ith group @
   z0 = z[ini:fini,.];     @ z for ith group @

   @ Combine x0 with z0 excluding intercepts @

   if cols(z0) eq 1;  @ case of only intercept @
   else;  x0 = x0~z0[.,1:cols(z0)-1];
   endif;

   @ sample adjustment matched with other estimation procedures @

   y0 = y0[2:rows(y0)];
   y1 = y1[2:rows(y1)];
   dy = dy[2:rows(dy)];
   x0 = x0[2:rows(x0),.];
   inpt = ones(rows(y0),1);

   @ Generation of the demeaned Y and X @

   tmpy  = y0 - inpt*invpd(inpt'inpt)*inpt'y0;  @ demeaned y @
   tmpy1 = y1 - inpt*invpd(inpt'inpt)*inpt'y1;  @ demeaned y1 @
   tmpdy = dy - inpt*invpd(inpt'inpt)*inpt'dy;  @ demeaned dy @
   tmpx  = x0 - inpt*invpd(inpt'inpt)*inpt'x0;  @ demeaned x @

@ whole demeaned samples used @

   if mxi eq 1;
     QY  = tmpy; QX = tmpx; QY1 = tmpy1; Qdy = tmpdy;
   else;
     QY  = QY|tmpy;    QX  = QX|tmpx;
     QY1 = QY1|tmpy1;  Qdy = Qdy|tmpdy;
   endif;

   @ Mean of Y and X for ith group @

   meanY[mxi]  = meanc(y0);           @ mean of Y  for group i @
   mnX = zeros(1,ntmp);               @ mean of X for group i @
   ei1 = 1;
   do until ei1 > k;
      mnX[1,ei1] = meanc(X0[.,ei1]);
      ei1 = ei1 + 1;
   endo;
   meanX[mxi,.] = mnX;                @ mean of X for group i @

/* ****************************************************
   Static Fixed Effects Estimation of the LR parameters
   **************************************************** */

   Fden  = x0'x0 - x0'inpt*invpd(inpt'inpt)*inpt'x0;
   Fnum  = x0'y0 - x0'inpt*invpd(inpt'inpt)*inpt'y0;

   sumden = sumden + Fden;   @ sum of denominator @
   sumnum = sumnum + Fnum;   @ sum of numerator @

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

Ftheta = invpd(sumden)*sumnum;   @ fixed effects estimator @

Qres  = Qy - Qx*Ftheta;                 @ residual @
nreg  = rows(Ftheta) + n;               @ # of parameters @
Qobs  = rows(Qy);                       @ observation used @
Qsig2 = Qres'Qres/(Qobs-nreg);      @ variance @
Qtemp = invpd(QX'QX);
Qvar  = Qsig2.*Qtemp;
Qse   = sqrt(diag(Qvar));                 @ standard error @
Qt    = Ftheta./Qse;                      @ t-ratios @

SFERES = Ftheta~Qse~Qt;

@  computation of R2 and diagnostic statistics and so on @

   reg = qy1~qx;
   dep = QdY;
   epara = -1|Ftheta;
   yhat = reg*epara;

   { SFEDRES } = DIAGNOS(dep,reg,Qres,yhat,nreg);

@ computation of robust standard errors of White @

sum = 0;
tini = 1;
tfini = 0;
ei = 1;
do until ei > n;

   diffT = T[ei]-1;
   tfini = tfini + diffT;
   Qdep = Qy[tini:tfini];       @ Tx1 @
   Qreg = Qx[tini:tfini,.];     @ Txk @
   tini = tini + diffT;

   Qu = Qdep - Qreg*Ftheta;     @ Tx1 @
   temp = Qreg'Qu*Qu'Qreg;      @ kxk @
   sum = sum + temp;

   ei = ei + 1;
endo;

RQvar = Qtemp*sum*Qtemp;
RQse  = sqrt(diag(RQvar));
RQt   = Ftheta./RQse;

RSFERES = Ftheta~RQse~RQt;

@ estimation of fixed effects @

Qind   = zeros(n,1);
Qindse = zeros(n,1);
ei = 1;
do until ei > n;
   Qind[ei] = meanY[ei] - meanx[ei,.]*Ftheta;  @ fixed effects @
   Qindse[ei] = sqrt(Qsig2/T[ei]);
   ei = ei + 1;
endo;
Qindt = Qind./Qindse;

SFEIND = Qind~Qindse~Qindt;   @ Results for the individual effects @

Ftheta = Ftheta[1:k];  @ collecting LR para on x only @

RETP(SFERES,SFEIND,RSFERES,SFEDRES,Qobs,Ftheta); ENDP;

/*  ***********************************
*   Case for the ML estimation option 1
*   *********************************** */

/*  ******************************************** */
    PROC(6) = DFE1(Y,X,Z,n,T,k,numot,plag,qlag);
/*  ******************************************************
*   Dynamic fixed effects estimation of the LR parameters
*
*   Model:
*
*   dy(i) = phi*y1(i) + beta*x(i) + a(i) + W*si + e(i)  (1)
*
*   estimate (1) by FE, and obtain the diagnostics,
*   then obtain the LR parameters by
*
*   theta = -beta/phi
*   *********************************** */

Local qy0,qdy,qx0,qy1,qww,meandep,dimind,meanind,ini,fini,mxi,tis,
      ytemp,xtemp,y0,dy,y1,x0,ww,coname,tmpy0,tmpdy,tmpy1,tmpx0,tmpww,
      freg,tmpf,fden,fnum,sumnum,sumden,qreg,ei,mnind,fpara,qtheta,qres,
      nreg,qsig2,qvar,qse,qt,qtemp,dferes,tini,tfini,sum,deptmp,indtmp,temp,qu,
      rqse,rqt,rqvar,qind,qindse,qindt,dfeind,ei1,difft,vtheta,
      beta,setheta,ttheta,faihat,dtheta,dtheta1,dtheta2,h,qobs,
      Rsetheta,Rttheta,rvtheta,dep,qecm,reg,yhat,qpara,res,inpt,wflag,
      RDFERES,dfedres,thres,rthres,pp,mxlag,ztemp;

@ Initialization @

pp = plag[1];
if     pp eq 0;
       dimind = k + numot;
elseif pp > 0; QY1 = 0;   @ when there are lagged dep vbles @
       dimind = k + 1 + numot;  @ dimension of indep variables
                                  excluding fixed effects @
endif;
sumden = zeros(dimind,dimind);
sumnum = zeros(dimind,1);

meandep = zeros(n,1);           @ mean of dep(1), ..., dep(n) @
meanind = zeros(n,dimind);    @ mean of ind(1), ..., ind(n) @

ini  = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @

mxi = 1;
do until mxi > N;

   Tis  = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ******************************************
*  Generation of the data used in estimation
*  See FDGP1
*  ****************************************** */

   ytemp = y[ini:fini];     @ y for ith group @
   xtemp = x[ini:fini,.];   @ x for ith group @
   ztemp = z[ini:fini,.];   @ z for ith group @

   {y0, dy, y1, x0, ww, wflag, mxlag}
   = FDGP1(ytemp,xtemp,ztemp,plag,qlag,k,tis,mxi);

   @ Generation of the demeaned ones @

   inpt = ones(rows(dy),1);

   tmpy0 = y0 - inpt*invpd(inpt'inpt)*inpt'y0;  @ demeaned y @
   tmpy1 = y1 - inpt*invpd(inpt'inpt)*inpt'y1;  @ demeaned y1 @
   tmpdy = dy - inpt*invpd(inpt'inpt)*inpt'dy;  @ demeaned dy @
   tmpx0 = x0 - inpt*invpd(inpt'inpt)*inpt'x0;  @ demeaned x @

   @ whole demeaned samples used @

   if mxi eq 1;
      QY0 = tmpy0; QX0 = tmpx0; QY1 = tmpy1; Qdy = tmpdy;
   else;
      QY1 = QY1|tmpy1;   Qdy = Qdy|tmpdy;
      QY0 = QY0|tmpy0;   QX0 = QX0|tmpx0;
   endif;

   if wflag eq 1;
      tmpww = ww - inpt*invpd(inpt'inpt)*inpt'ww;  @ demeaned ww @
      if mxi eq 1; QWW = tmpww;
      else;        QWW = QWW|tmpww;
      endif;
   endif;

/* *****************************************************
*  Dynamic Fixed Effects Estimation of the SR parameters
*  ***************************************************** */

   if pp eq 0;    @ when no lagged dep vble @
      if     wflag eq 0; Freg = x0;
      elseif wflag eq 1; Freg = x0~ww;  @ regressors @
      endif;
   else;
      if     wflag eq 0; Freg = y1~x0;
      elseif wflag eq 1; Freg = y1~x0~ww;  @ regressors @
      endif;
   endif;

   tmpf = Freg - inpt*invpd(inpt'inpt)*inpt'Freg;
   if mxi eq 1;  Qreg = tmpf;
   else;         Qreg = Qreg|tmpf;  @ whole demeaned regressors @
   endif;

   Fden  = Freg'Freg - Freg'inpt*invpd(inpt'inpt)*inpt'Freg;

   if pp eq 0; Fnum  = Freg'y0 - Freg'inpt*invpd(inpt'inpt)*inpt'y0;
   else;       Fnum  = Freg'dy - Freg'inpt*invpd(inpt'inpt)*inpt'dy;
   endif;

   sumden = sumden + Fden;   @ sum of denominator @
   sumnum = sumnum + Fnum;   @ sum of numerator @

@ Mean of dep and indep vatiabels for ith group @

   meandep[mxi] = meanc(dy);         @ mean of dep vble for group i @
   mnind = zeros(1,dimind);          @ mean of indep vbles for group i @
   ei1 = 1;
   do until ei1 > dimind;
      mnind[1,ei1] = meanc(Freg[.,ei1]);
      ei1 = ei1 + 1;
   endo;
   meanind[mxi,.] = mnind;

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

@ fixed effects estimator of all parameters @

Fpara  = invpd(sumden)*sumnum;

@ Fixed effects estimator of LR parameters @

If     pp eq 0;   Qtheta = Fpara[1:k];
elseif pp eq 0;   Qtheta = Fpara[2:k+1];
endif;

if pp eq 0; Qres = Qy0 - Qreg*Fpara;            @ residual @
else;       Qres = Qdy - Qreg*Fpara;            @ residual @
endif;

Qobs  = rows(Qdy);                   @ observations @
nreg  = rows(Fpara) + n;             @ # of parameters @
Qsig2 = Qres'Qres/(Qobs-nreg);  @ variance @
Qtemp = invpd(Qreg'Qreg);
Qvar  = Qsig2.*Qtemp;

Qse   = sqrt(diag(Qvar));            @ standard error @
Qt    = Fpara./Qse;                  @ t-ratios @

DFERES = Fpara~Qse~Qt;

@ computation of robust standard errors of White @

sum = 0;
tini = 1;
tfini = 0;
ei = 1;
do until ei > n;

   if mxlag eq 0; diffT = T[ei] - 1;
   else;          diffT = T[ei] - mxlag;
   endif;
   tfini = tfini + diffT;
   if mxlag eq 0; deptmp = Qy0[tini:tfini];     @ Tx1 @
   else;          deptmp = Qdy[tini:tfini];     @ Tx1 @
   endif;
   indtmp = Qreg[tini:tfini,.];   @ Tx(# of para) @
   tini = tini + diffT;

   Qu = deptmp - indtmp*Fpara;    @ Tx1 @
   temp = indtmp'Qu*Qu'indtmp;    @ (# of para)x(# of para) @
   sum = sum + temp;

   ei = ei + 1;
endo;

RQvar = Qtemp*sum*Qtemp;
RQse  = sqrt(diag(RQvar));
RQt   = Fpara./RQse;

RDFERES = Fpara~RQse~RQt;

@ estimation of fixed effects @

Qind = zeros(n,1);
Qindse = zeros(n,1);
ei = 1;
do until ei > n;
   Qind[ei] = meandep[ei] - meanind[ei,.]*Fpara;
   Qindse[ei] = sqrt(Qsig2/T[ei]);
   ei = ei + 1;
endo;
Qindt = Qind./Qindse;

DFEIND = Qind~Qindse~Qindt;   @ Results for the individual fixed effects @

@  computation of the LR parameters @

if pp eq 0;     @ no lagged dep vble @

      beta   = Fpara[1:k];        @ beta(i) @
      Qtheta = beta;

      setheta = Qse[1:k];         @ std. err. of theta @
      ttheta  = Qtheta./setheta;  @ t-ratios of theta @
      Thres = Qtheta~setheta~ttheta;

      Rsetheta = RQse[1:k];         @ std. err. of theta @
      Rttheta  = Qtheta./Rsetheta;  @ t-ratios of theta @
      RThres = Qtheta~Rsetheta~Rttheta;

   else;

      faihat = Fpara[1];          @ fai(i) @
      beta = Fpara[2:k+1];        @ beta(i) @

      Qtheta = -beta./faihat;     @ theta(i) @

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

    @ Using OLS @

      vtheta  = dtheta*Qvar[1:k+1,1:k+1]*dtheta'; @ variance of theta @
      setheta = sqrt(diag(vtheta));               @ std. err. of theta @
      ttheta  = Qtheta./setheta;                  @ t-ratios of theta @

      Thres = Qtheta~setheta~ttheta;

    @ Using Robust @

      Rvtheta  = dtheta*RQvar[1:k+1,1:k+1]*dtheta'; @ variance of theta @
      Rsetheta = sqrt(diag(Rvtheta));               @ std. err. of theta @
      Rttheta  = Qtheta./Rsetheta;                  @ t-ratios of theta @

      RThres = Qtheta~Rsetheta~Rttheta;

   endif;

@  computation of R2 and diagnostic statistics and so on using UR ECM @

dep = QdY;

IF pp eq 0; Qpara = -1|Fpara;
else;       Qpara = Fpara;
endif;

if     wflag eq 0;  reg = QY1~QX0;
elseif wflag eq 1;  reg = QY1~QX0~Qww;
endif;

yhat = reg*Qpara;
res = dep - yhat;

{ DFEDRES } = DIAGNOS(dep,reg,res,yhat,nreg);

@  Save the estimation results @

DFERES = Thres|DFERES;
@ in the order of theta, phi, beta and other excluding intercepts @
RDFERES = RThres|RDFERES;

RETP(DFERES,DFEIND,RDFERES,DFEDRES,Qobs,QTHETA);ENDP;

/* ****************************************************** */
   PROC (7) = FDGP1(ytemp,xtemp,ztemp,pp,qlag,k,tis,mxi);
/* ******************************************************
*  Generation of the data used in estimation
*
*  y0: current value of dep vble
*  dy: first difference of dep vble
*  y1: lagged value of dep vble
*  x0: current value of indep vbles, x
*  ww: dy's, dx1's and other fixed regressors excluding intercept.
*  *************************************************************** */

Local y1temp,dytemp,x1temp,dxtemp,dy1temp,dy2temp,dy3temp,dystemp,wflag,
      dx0temp,dx1temp,dx2temp,dx3temp,dxtmp,dxflag,
      inpt,ltemp,mxlag,sini,sig2,mxj,y0,y1,dy,x0,ww;

 @ Generation y, dy, y(-1), x and ww @

   y1temp  = lag(ytemp);       @ y(-1) @
   dytemp  = ytemp - y1temp;   @ dy @

   x1temp  = lag(xtemp);       @ x(-1) @
   dxtemp  = xtemp - x1temp;   @ dx @

 @ Generation of W such as lagged differences and intercept @

 @ Generation of lagged dy's @

   if     pp eq 2; dy1temp = lag(dytemp);
                          dystemp = dy1temp;
   elseif pp eq 3; dy1temp = lag(dytemp); dy2temp = lag(dy1temp);
                          dystemp = dy1temp~dy2temp;
   elseif pp eq 4; dy1temp = lag(dytemp); dy2temp = lag(dy1temp);
                          dy3temp = lag(dy2temp);
                          dystemp = dy1temp~dy2temp~dy3temp;
   endif;

   if pp > 1;  WW = dystemp; wflag = 1;
   else;       wflag = 0;    @ W's empty yet @
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

   if wflag eq 0;  if dxflag eq 0;
                   else;  WW = dxtmp; wflag = 1;
                   endif;
   else;           if dxflag eq 0;
                   else;  WW = WW~dxtmp;
                   endif;
   endif;

@ Including other fixed regressors @

   if wflag eq 0;  if cols(ztemp) eq 1;  @ case of intercept only @
                   else;  WW = ztemp[.,1:cols(ztemp)-1]; wflag = 1;
                   endif;
   else;           if cols(ztemp) eq 1;
                   else;  WW = WW~ztemp[.,1:cols(ztemp)-1];
                   endif;
   endif;

 @ sample adjustment @

   ltemp = pp|qlag[mxi,.]';
   mxlag = maxc(ltemp);      @ maximum lag order @

   if mxlag eq 0; sini = 2;
   else;          sini = mxlag + 1;
   endif;

   y0 = ytemp[sini:Tis];        @ y @
   y1 = y1temp[sini:Tis];       @ y(-1) @
   dy = dytemp[sini:Tis];       @ dy @

   x0 = xtemp[sini:Tis,.];      @ x @
   if wflag eq 0; ww = 0;
   else; ww = ww[sini:Tis,.];   @ other variables excluding intercept @
   endif;

RETP(y0,dy,y1,x0,ww,wflag,mxlag); ENDP;

/*  ***********************************
*   Case for the ML estimation option 2
*   *********************************** */

/*  *********************************************************** */
    PROC(7) = DFE2(Y,X1,X2,Z,n,T,k1,k2,numot,plag,q1lag,q2lag);
/*  ***********************************************************
*   Dynamic fixed effects estimation of the LR parameters
*
*   Model:
*
*   dy(i) = phi*y1(i) + beta*x(i) + a(i) + W*si + e(i)  (1)
*
*   estimate (1) by FE, and obtain the diagnostics
*   then obtain the LR parameters by
*
*   theta1 = -beta1/phi
*   theta2 = -beta2/phi
*   *********************************** */

Local qy0,qdy,qx0,qy1,qww,meandep,dimind,meanind,ini,fini,mxi,tis,
      ytemp,xtemp,y0,dy,y1,x0,ww,coname,tmpy0,tmpdy,tmpy1,tmpx0,tmpww,
      freg,tmpf,fden,fnum,sumnum,sumden,qreg,ei,mnind,fpara,qtheta,qres,
      nreg,qsig2,qvar,qse,qt,qtemp,dferes,tini,tfini,sum,deptmp,indtmp,temp,qu,
      rqse,rqt,rqvar,qind,qindse,qindt,dfeind,ei1,difft,vtheta,
      beta,setheta,ttheta,faihat,dtheta,dtheta1,dtheta2,h,qobs,
      Rsetheta,Rttheta,rvtheta,dep,qecm,reg,yhat,qpara,res,inpt,wflag,
      thres,rthres,pp,x1temp,x2temp,x10,x20,tmpx10,tmpx20,qx10,qx20,beta1,beta2,
      RDFERES,dfedres,qtheta1,qtheta2,mxlag,ztemp;

@ Initialization @

pp = plag[1];
if     pp eq 0;
       dimind = k1 + k2 + numot;
elseif pp > 0; QY1 = 0;   @ when there are lagged dep vbles @
       dimind = k1 + k2 + 1 + numot;  @ dimension of indep variables
                                        excluding fixed effects @
endif;
sumden = zeros(dimind,dimind);
sumnum = zeros(dimind,1);

meandep = zeros(n,1);           @ mean of dep(1), ..., dep(n) @
meanind = zeros(n,dimind);    @ mean of ind(1), ..., ind(n) @

ini  = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @

mxi = 1;
do until mxi > N;

   Tis  = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

/* ******************************************
*  Generation of the data used in estimation
*  See FDGP2
*  ****************************************** */

   ytemp  = y[ini:fini];        @ y for ith group @
   x1temp = x1[ini:fini,.];     @ x1 for ith group @
   x2temp = x2[ini:fini,.];     @ x2 for ith group @
   ztemp  = z[ini:fini,.];      @ z for ith group @

   {y0, dy, y1, x10, x20, ww, wflag, mxlag }
   = FDGP2(ytemp,x1temp,x2temp,ztemp,plag,q1lag,q2lag,k1,k2,mxi,tis);

   @ Generation of the demeaned ones @

   inpt = ones(rows(dy),1);

   tmpy0 = y0 - inpt*invpd(inpt'inpt)*inpt'y0;  @ demeaned y @
   tmpy1 = y1 - inpt*invpd(inpt'inpt)*inpt'y1;  @ demeaned y1 @
   tmpdy = dy - inpt*invpd(inpt'inpt)*inpt'dy;  @ demeaned dy @

   tmpx10 = x10 - inpt*invpd(inpt'inpt)*inpt'x10;  @ demeaned x1 @
   tmpx20 = x20 - inpt*invpd(inpt'inpt)*inpt'x20;  @ demeaned x2 @

   @ whole demeaned samples used @

   if mxi eq 1;
      QY0 = tmpy0; QX10 = tmpx10; QY1 = tmpy1; Qdy = tmpdy; QX20 = tmpx20;
   else;
      QY1 = QY1|tmpy1;   Qdy = Qdy|tmpdy;
      QY0 = QY0|tmpy0;   QX10 = QX10|tmpx10;   QX20 = QX20|tmpx20;
   endif;

   if wflag eq 1;
      tmpww = ww - inpt*invpd(inpt'inpt)*inpt'ww;  @ demeaned ww @
      if mxi eq 1; QWW = tmpww;
      else;        QWW = QWW|tmpww;
      endif;
   endif;

/* *****************************************************
*  Dynamic Fixed Effects Estimation of the SR parameters
*  ***************************************************** */

   if pp eq 0;    @ when no lagged dep vble @
      if     wflag eq 0; Freg = x10~x20;
      elseif wflag eq 1; Freg = x10~x20~ww;  @ regressors @
      endif;
   else;
      if     wflag eq 0; Freg = y1~x10~x20;
      elseif wflag eq 1; Freg = y1~x10~x20~ww;  @ regressors @
      endif;
   endif;

   tmpf = Freg - inpt*invpd(inpt'inpt)*inpt'Freg;
   if mxi eq 1;  Qreg = tmpf;
   else;         Qreg = Qreg|tmpf;  @ whole demeaned regressors @
   endif;

   Fden  = Freg'Freg - Freg'inpt*invpd(inpt'inpt)*inpt'Freg;

   if pp eq 0; Fnum  = Freg'y0 - Freg'inpt*invpd(inpt'inpt)*inpt'y0;
   else;       Fnum  = Freg'dy - Freg'inpt*invpd(inpt'inpt)*inpt'dy;
   endif;

   sumden = sumden + Fden;   @ sum of denominator @
   sumnum = sumnum + Fnum;   @ sum of numerator @

@ Mean of dep and indep vatiabels for ith group @

   meandep[mxi] = meanc(dy);         @ mean of dep vble for group i @
   mnind = zeros(1,dimind);          @ mean of indep vbles for group i @
   ei1 = 1;
   do until ei1 > dimind;
      mnind[1,ei1] = meanc(Freg[.,ei1]);
      ei1 = ei1 + 1;
   endo;
   meanind[mxi,.] = mnind;

   ini = ini + Tis;
   mxi = mxi + 1;
endo;

@ fixed effects estimator of all parameters @

Fpara  = invpd(sumden)*sumnum;

if pp eq 0; Qres = Qy0 - Qreg*Fpara;            @ residual @
else;       Qres = Qdy - Qreg*Fpara;            @ residual @
endif;

Qobs  = rows(Qdy);                   @ observations @
nreg  = rows(Fpara) + n;             @ # of parameters @
Qsig2 = Qres'Qres/(Qobs-nreg);       @ variance @
Qtemp = invpd(Qreg'Qreg);
Qvar  = Qsig2.*Qtemp;

Qse   = sqrt(diag(Qvar));            @ standard error @
Qt    = Fpara./Qse;                  @ t-ratios @

DFERES = Fpara~Qse~Qt;

@ computation of robust standard errors of White @

sum = 0;
tini = 1;
tfini = 0;
ei = 1;
do until ei > n;

   if mxlag eq 0; diffT = T[ei] - 1;
   else;          diffT = T[ei] - mxlag;
   endif;
   tfini = tfini + diffT;
   if mxlag eq 0; deptmp = Qy0[tini:tfini];       @ Tx1 @
   else;          deptmp = Qdy[tini:tfini];      @ Tx1 @
   endif;
   indtmp = Qreg[tini:tfini,.];   @ Tx(# of para) @
   tini = tini + diffT;

   Qu = deptmp - indtmp*Fpara;    @ Tx1 @
   temp = indtmp'Qu*Qu'indtmp;    @ (# of para)x(# of para) @
   sum = sum + temp;

   ei = ei + 1;
endo;

RQvar = Qtemp*sum*Qtemp;
RQse  = sqrt(diag(RQvar));
RQt   = Fpara./RQse;

RDFERES = Fpara~RQse~RQt;

@ estimation of fixed effects @

Qind = zeros(n,1);
Qindse = zeros(n,1);
ei = 1;
do until ei > n;
   Qind[ei] = meandep[ei] - meanind[ei,.]*Fpara;
   Qindse[ei] = sqrt(Qsig2/T[ei]);
   ei = ei + 1;
endo;
Qindt = Qind./Qindse;

DFEIND = Qind~Qindse~Qindt;   @ Results for the individual fixed effects @

@  computation of the LR parameters @

if pp eq 0;     @ no lagged dep vble @

      beta1 = Fpara[1:k1];        @ beta1(i) @
      beta2 = Fpara[k1+1:k1+k2];  @ beta2(i) @
      Qtheta1 = beta1;
      Qtheta2 = beta2;
      Qtheta = Qtheta1|Qtheta2;

      setheta = Qse[1:k1+k2];         @ std. err. of theta1 and 2 @
      ttheta = Qtheta./setheta;       @ t-ratios of theta1 and 2 @
      Thres = Qtheta~setheta~ttheta;

      rsetheta = rQse[1:k1+k2];         @ std. err. of theta1 and 2 @
      rttheta = Qtheta./rsetheta;       @ t-ratios of theta1 and 2 @
      RThres = Qtheta~Rsetheta~Rttheta;

else;

      faihat = Fpara[1];            @ fai(i) @
      beta1 = Fpara[2:k1+1];        @ beta(i) @
      beta2 = Fpara[k1+2:k1+k2+1];  @ beta(i) @

      Qtheta1 = -beta1./faihat;     @ theta1(i) @
      Qtheta2 = -beta2./faihat;     @ theta2(i) @

      @  Computation of the standard errors and t-values
        for the LR parameters using the delta method @

      beta = beta1|beta2;
      Qtheta = Qtheta1|Qtheta2;
      dtheta1 = zeros(k1+k2,1);  @ first differential  w.r.t. phi @
      h = 1;
      do until h > k1+k2;
         dtheta1[h,1] = beta[h]/faihat^2;
         h = h + 1;
      endo;

      dtheta2 = eye(k1+k2);     @ first differential w.r.t beta @
      dtheta2 = -(1/faihat).*dtheta2;
      dtheta  = dtheta1~dtheta2;

    @ Using OLS @

      vtheta  = dtheta*Qvar[1:k1+k2+1,1:k1+k2+1]*dtheta'; @ variance of theta @
      setheta = sqrt(diag(vtheta));               @ std. err. of theta @
      ttheta  = Qtheta./setheta;                  @ t-ratios of theta @

    @ Using Robust @

      Rvtheta  = dtheta*RQvar[1:k1+k2+1,1:k1+k2+1]*dtheta';
      Rsetheta = sqrt(diag(Rvtheta));               @ std. err. of theta @
      Rttheta  = Qtheta./Rsetheta;                  @ t-ratios of theta @

      Thres = Qtheta~setheta~ttheta;
      RThres = Qtheta~Rsetheta~Rttheta;

endif;

@  computation of R2 and diagnostic statistics and so on @

if     wflag eq 0;  reg = QY1~QX10~QX20;
elseif wflag eq 1;  reg = QY1~QX10~QX20~Qww;
endif;
IF pp eq 0; Qpara = -1|Fpara;
ELSE;       Qpara = Fpara;
ENDIF;

dep = QdY;
yhat = reg*Qpara;
res = dep - yhat;

{ DFEDRES } = DIAGNOS(dep,reg,res,yhat,nreg);

@  Save the estimation results @

DFERES = Thres|DFERES;
@ in the order of theta1, theta2, phi, beta1, beta2
and other excluding intercepts @
RDFERES = RThres|RDFERES;

RETP(DFERES,DFEIND,RDFERES,DFEDRES,Qobs,QTHETA1,Qtheta2);ENDP;

/* ************************************************************************* */
   PROC (8) = FDGP2(ytemp,x1temp,x2temp,ztemp,plag,q1lag,q2lag,k1,k2,mxi,tis);
/* ***************************************************************************
*  Generation of the data used in estimation
*
*  y0: current value of dep vble
*  dy: first difference of dep vble
*  y1: lagged value of dep vble
*  x10: current value of indep vbles, x1
*  x20: current value of indep vbles, x2
*  ww: dy's, dx1's, dx2's and other fixed regressors excluding intercept
*  ********************************************************************* */

Local y1temp,dytemp,x11temp,dx1temp,x21temp,dx2temp,
      dy1temp,dy2temp,dy3temp,dystemp,wflag,
      dx10temp,dx11temp,dx12temp,dx13temp,dx1tmp,dx1flag,
      dx20temp,dx21temp,dx22temp,dx23temp,dx2tmp,dx2flag,
      inpt,ltemp,mxlag,sini,sig2,mxj,
      y0,y1,dy,x10,x20,ww;

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

 @ Including other fixed regressors excluding intercept @

   if wflag eq 0;  if cols(ztemp) eq 1;
                   else;  WW = ztemp[.,1:cols(ztemp)-1]; wflag = 1;
                   endif;
   else;           if cols(ztemp) eq 1;
                   else;  WW = WW~ztemp[.,1:cols(ztemp)-1];
                   endif;
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

   if     wflag eq 0; ww = 0;
   elseif wflag eq 1; ww = ww[sini:Tis,.];
          @ other variables excluding intercept @
   endif;

RETP(y0,dy,y1,x10,x20,ww,wflag,mxlag); ENDP;

