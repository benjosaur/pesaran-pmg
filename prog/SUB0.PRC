@  ********************************************************* @
@* Common sub-procedures for the main PMG estimation program @
@  ********************************************************* @

/*  ******************** */
    PROC(5) = BEG(n,nz);
/*  ******************** */

Local beg_op,k,yname,xname,zname;

@* ***************************** @
@* beg_op: the beginning options @
@* ***************************** @

@ ****** @
  BACK0:
@ ****** @

output off; cls;

"To compute the Pooled Mean Group Estimates (PMGE) you first need to
supply the program with information concerning your panel data
model such as the dependent variable, the regressors, the groups to be
included in the estimation, and so on.

Notice that it is assumed that you have already correctly loaded the
data, made the necessary transformations, and have correctly set up
the pooled regression model.

To simplify the specification when the same regression model is used,
you can skip this initial stage by typing 2, otherwise
type 1 to start with a new model.";"";

beg_op = con(1,1);"";
if beg_op eq 1;

   {k, Yname, Xname, Zname } = INITIAL(n,nz);
   save k Yname Xname Zname;

elseif beg_op eq 2;

   load k Yname Xname Zname;

else; "WARNING! Try again.";
      "Press any key to continue."; wait; cls;
      goto BACK0;
endif;

RETP(beg_op,k,Yname,Xname,Zname);
ENDP;

/*  ************************ */
    Proc(4) = INITIAL(n,nz);
/*  *********************************************************
*   Informations on the data process which will be used later
*
*   k: number of the regressors excluding intercepts
*   Yname: the name of dependent variable
*   Xname: the name of independent variables
*   ******************************************************* */

Local k,Yname,vname,Xname,Zname,kj,kn,kkn,vvname;

screen on;

@ ****** @
  BACK1:
@ ****** @
output off; cls;

"Type the number (maximum 10) of regressors (the number of columns of the X
matrix) specified in the data loading section of the main program.";"";

k = con(1,1);

if k > 10; "WARNING! Try again.";
      "Press any key to continue."; wait; cls;
      goto BACK1;
endif;

"Type the name (maximum 4 letters) you wish to give to the dependent
variable - this name will be used when reporting the results.";"";

vname = cons;"";

format /rd 3,0;
kj = 1;
do until kj > k;

   if     kj eq 1;
      "Type the name (maximum 4 letters) of the" kj "st regressor of X.";
   elseif kj eq 2;
      "Type the name (maximum 4 letters) of the" kj "nd regressor of X.";
   elseif kj eq 3;
      "Type the name (maximum 4 letters) of the" kj "rd regressor of X.";
   else;
      "Type the name (maximum 4 letters) of the" kj "th regressor of X.";
   endif;

   vvname = cons;"";  vname = vname|vvname;
   kj = kj + 1;
endo;
Yname = vname[1];
Xname = vname[2:rows(vname)];cls;

if     nz eq 1;  zname = "Inpt"|"";
elseif nz > 1;
   format /rd 3,0;
   kj = 1;
   do until kj > nz;
      "Type the name (maximum 8 letters) of the" kj "th regressor of Z.";
      vvname = cons;"";
      if kj eq 1; zname = vvname;
      else;       zname = zname|vvname;
      endif;
      kj = kj + 1;
   endo;
endif;

RETP(k,Yname,Xname,Zname);ENDP;

/*  ****************************** */
    Proc(4) = DETREGR(maxt,maxn);
/*  *********************************************************
*   Construction of deterministic regressors such as constants,
*   time trends and seasonal dummies
*
*   The dimension of each is maxt*maxn by 1.
*
*   inpt: ones
*   trend: time trends
*   qs1, qs2, qs3: quarterly seasnal dummies
*   ms1, ..., ms11: monthly seasnal dummies
*   ******************************************************* */

Local temp,inpt,trend,qs1tmp,qs2tmp,qs3tmp,qxi,qmaxt,ix,mmaxt,trtmp,
      qs1,qs2,qs3,qs,ms,mxi,
      ms1tmp,ms2tmp,ms3tmp,ms4tmp,ms5tmp,ms6tmp,
      ms7tmp,ms8tmp,ms9tmp,ms10tmp,ms11tmp,
      ms1,ms2,ms3,ms4,ms5,ms6, ms7,ms8,ms9,ms10,ms11;

temp = maxt*maxn;
inpt = ones(temp,1);    @ intercept @

qs1tmp = zeros(maxt,1); @ quarterly seasonal dummy1 @
qs2tmp = zeros(maxt,1); @ quarterly seasonal dummy2 @
qs3tmp = zeros(maxt,1); @ quarterly seasonal dummy3 @

qmaxt = ceil(maxt/4);
qxi = 1;
do until qxi > qmaxt;

    ix = 4*(qxi-1);
    if ix + 1 > maxt; goto QSout;
    else;             qs1tmp[ix+1] = 1;
    endif;
    if ix + 2 > maxt; goto QSout;
    else;             qs2tmp[ix+2] = 1;
    endif;
    if ix + 3 > maxt; goto QSout;
    else;             qs3tmp[ix+3] = 1;
    endif;

   qxi = qxi + 1;
endo;

@ ****** @
  QSout:
@ ****** @

ms1tmp = zeros(maxt,1);  @ monthly seasonal dummy1 @
ms2tmp = zeros(maxt,1);  @ monthly seasonal dummy2 @
ms3tmp = zeros(maxt,1);  @ monthly seasonal dummy3 @
ms4tmp = zeros(maxt,1);  @ monthly seasonal dummy4 @
ms5tmp = zeros(maxt,1);  @ monthly seasonal dummy5 @
ms6tmp = zeros(maxt,1);  @ monthly seasonal dummy6 @
ms7tmp = zeros(maxt,1);  @ monthly seasonal dummy7 @
ms8tmp = zeros(maxt,1);  @ monthly seasonal dummy8 @
ms9tmp = zeros(maxt,1);  @ monthly seasonal dummy9 @
ms10tmp = zeros(maxt,1); @ monthly seasonal dummy10 @
ms11tmp = zeros(maxt,1); @ monthly seasonal dummy11 @

mmaxt = ceil(maxt/12);
qxi = 1;
do until qxi > mmaxt;

   ix = 12*(qxi-1);
   if ix + 1 > maxt; goto MSout;
   else;             ms1tmp[ix+1] = 1;
   endif;
   if ix + 2 > maxt; goto MSout;
   else;             ms2tmp[ix+2] = 1;
   endif;
   if ix + 3 > maxt; goto MSout;
   else;             ms3tmp[ix+3] = 1;
   endif;
   if ix + 4 > maxt; goto MSout;
   else;             ms4tmp[ix+4] = 1;
   endif;
   if ix + 5 > maxt; goto MSout;
   else;             ms5tmp[ix+5] = 1;
   endif;
   if ix + 6 > maxt; goto MSout;
   else;             ms6tmp[ix+6] = 1;
   endif;
   if ix + 7 > maxt; goto MSout;
   else;             ms7tmp[ix+7] = 1;
   endif;
   if ix + 8 > maxt; goto MSout;
   else;             ms8tmp[ix+8] = 1;
   endif;
   if ix + 9 > maxt; goto MSout;
   else;             ms9tmp[ix+9] = 1;
   endif;
   if ix + 10 > maxt; goto MSout;
   else;              ms10tmp[ix+10] = 1;
   endif;
   if ix + 11 > maxt; goto MSout;
   else;              ms11tmp[ix+11] = 1;
   endif;

   qxi = qxi + 1;
endo;

@ ****** @
  MSout:
@ ****** @

trtmp = seqa(1,1,maxt); @ time trend @

mxi = 1;
do until mxi > maxn;

   if mxi eq 1; trend = trtmp;
                qs1 = qs1tmp; qs2 = qs2tmp; qs3 = qs3tmp;
                ms1 = ms1tmp; ms2 = ms2tmp; ms3 = ms3tmp;
                ms4 = ms4tmp; ms5 = ms5tmp; ms6 = ms6tmp;
                ms7 = ms7tmp; ms8 = ms8tmp; ms9 = ms9tmp;
                ms10 = ms10tmp; ms11 = ms11tmp;
   else;        trend = trend|trtmp;
                qs1 = qs1|qs1tmp; qs2 = qs2|qs2tmp; qs3 = qs3|qs3tmp;
                ms1 = ms1|ms1tmp; ms2 = ms2|ms2tmp; ms3 = ms3|ms3tmp;
                ms4 = ms4|ms4tmp; ms5 = ms5|ms5tmp; ms6 = ms6|ms6tmp;
                ms7 = ms7|ms7tmp; ms8 = ms8|ms8tmp; ms9 = ms9|ms9tmp;
                ms10 = ms10|ms10tmp; ms11 = ms11|ms11tmp;
   endif;

   mxi = mxi + 1;
endo;

QS = qs1~qs2~qs3;
MS = ms1~ms2~ms3~ms4~ms5~ms6~ms7~ms8~ms9~ms10~ms11;

RETP(inpt,trend,QS,MS);
ENDP;

/*  ******************************************* */
    PROC (6) = SUBGR(dat,k,n,maxt,grname,gn_op);
/*  *******************************************
*   Selecting the sub-group data if necessary
*   ******************************************* */

Local n_ex,n_pos,n_subgr,na_subgr,jjunk,jj0,jflag,jtemp,jj1;

cls;

"You may exclude a subset of the groups from the analysis.
Type the number of groups you wish to exclude.";"";

n_ex = con(1,1);"";
format /rd 3,0;
"Type the" n_ex "numeric position(s) of the group(s) you wish to exclude
from the whole data set.";"";

n_pos = con(n_ex,1);
n_pos = sortc(n_pos,1);

n_subgr  = 0;          @ group numbers included in the subset of data @
na_subgr = "ini";      @ group names included in the subset of data @
jjunk = zeros(1,k+1);  @ initialization @
jj0 = 1;
do until jj0 > n;

    jtemp = dat[(jj0-1)*MAXT+1:jj0*MAXT,.];

    jflag = 0;   @ check if the jj0-th group is excluded @
    jj1 = 1;
    do until jj1 > n_ex;
       if jj0 eq n_pos[jj1]; jflag = 1;
       endif;
       jj1 = jj1 + 1;
    endo;

    if jflag eq 0; jjunk = jjunk|jtemp;
                   n_subgr = n_subgr|jj0;
                   if gn_op eq 2; na_subgr = na_subgr|grname[jj0];
                   endif;
    endif;
    jj0 = jj0 + 1;
endo;

dat = jjunk[2:rows(jjunk),.];
n_subgr = n_subgr[2:rows(n_subgr)];
if gn_op eq 2;   na_subgr = na_subgr[2:rows(na_subgr)];
endif;
n = n - n_ex;   @ number of groups included in the subset of the data @

RETP(dat,n_subgr,na_subgr,n_ex,n_pos,n);ENDP;

/*  ****************************************  */
    PROC (1) = DEMEAN(odata,maxt,N,missing);
/*  ****************************************  */

@ Transforms the original data nto the (cross-sectional) demeaned ones,
*  and deals with missing values @

Local  dmdata,h3,h4,h5,h6,temp,meantemp,demdata,sumtemp,cnttemp,
       tdata,h1,h2,h2a,h2b,h2c;

tdata = zeros(n,maxt);
h1 = 1;            @ index for i = 1,...,N @
h2 = 1;            @ left index of tdata @
h2a = maxt;        @ right index of tdata @
do until h1 > n;
   tdata[h1,.] = odata[h2:h2a,1]';
   h2 = h2 + maxT;
   h2a = h2a + maxT;
   h1 = h1 + 1;
endo;

dmdata = zeros(n,maxt);
@ dmdata: N by maxT matrix of demeaned data matrix including missing values @

h3 = 1;
do until h3 > maxt;

   temp  = tdata[.,h3];     @ nx1 vector for each t @
   sumtemp = 0;
   cnttemp = 0;    @ number of groups for each t excluding missing values @
   h4 = 1;
   do until h4 > n;
      if temp[h4] eq missing;
      else;  sumtemp = sumtemp + temp[h4];
             cnttemp = cnttemp + 1;
      endif;
      h4 = h4 + 1;
   endo;
   if cnttemp eq 0;
   else; meantemp = sumtemp/cnttemp;  @ cross section mean for each t @
   endif;

   h5 = 1;
   do until h5 > n;
      if temp[h5] eq missing; dmdata[h5,h3] = missing;
      else;                   dmdata[h5,h3] = temp[h5] - meantemp;
      endif;
      h5 = h5 + 1;
   endo;

   h3 = h3 + 1;
endo;

@ vectoriztion of dmdata @

dmdata = vec(dmdata');

@ N*maxT-vector of demeaned data including missing values @

RETP(dmdata);ENDP;

/*  ***********************************************  */
    PROC (3) = DELMISS(odata,dmdat,maxt,N,missing);
/*  ***********************************************
*   Transforms the original data into ones without missing values.
*   which are sum(Ti)*N by (k+1) matrix to be used in the estimation
*   Also count Ti's for each group.
*
*   odata: orinal data: (maxT*N) by (1+k) matrix,
*   dmdat: demeaned data: (maxT*N) by (1+k) matrix,
*   maxT: maximum number of time periods,
*   N: number of cross-section units,
*   missing: missing values: default is 8934567,
*   tdata: sum(Ti) by (k+1) data matrix to be used in estimation */

Local  rdata,ddata,h1,h2,h3,h4,h5,cnttemp,TI,cflag,nvar,rtmp,dtmp;

nvar = cols(odata);    @ number of variables @
rdata = zeros(1,nvar);
ddata = zeros(1,nvar);
TI = zeros(n,1);

h1 = 1;            @ index for i = 1,...,N @
h2 = 1;            @ left index of tdata @
h3 = maxt;        @ right index of tdata @
do until h1 > n;

   rtmp = odata[h2:h3,.];  @ T by (k+1) data matrix for each group @
   dtmp = dmdat[h2:h3,.];  @ T by (k+1) data matrix for each group @
   cnttemp = 0;
   h4 = 1;
   do until h4 > maxt;

      cflag = 1;
      h5 = 1;
      do until h5 > nvar;

         if rtmp[h4,h5] eq missing; cflag = 0;
         endif;
         h5 = h5 + 1;
      endo;

      if cflag eq 1; cnttemp = cnttemp + 1;
                     rdata = rdata|rtmp[h4,.];
                     ddata = ddata|dtmp[h4,.];
      endif;

      h4 = h4 + 1;
   endo;

   TI[h1] = cnttemp;  @ number of T's for each group @

   h1 = h1 + 1;
   h2 = h2 + maxT;
   h3 = h3 + maxT;
endo;

if rows(rdata) > 1; rdata = rdata[2:rows(rdata),.];
                    ddata = ddata[2:rows(ddata),.];
else; "No data available";
endif;

RETP(rdata,ddata,TI); ENDP;

/*  ******************************************* */
    PROC (1) = DIAGNOS(dep,reg,res,yhat,nreg);
/*  ******************************************************************
*   This procedure computes the followings using the ECM specification
*
*        dy(i) = phi(i)*y(i,-1) + theta(i)*x(i) + W(i)*si(i) + e(i),
*
*  (1) R-squares and Rbar-squared
*  (2) Maximized log-likelihood, AIC, SBC, HQ
*  (3) Diagnostic statistics such as
*
*      Godfrey's test of residual serial correlation: chi(1)
*      Ramsey's RESET test of functional form: chi(1)
*      Jarque-Bera's test of the normality: chi(2)
*      Test of Heterosckedasticity: chi(1)
*
*  see Microfit manual
*  ********************************************** */

Local obs,sig2,sige,LL,AIC,SBC,HQ,Rbar,RSQ,ysq2,rbarsq,
      e,x,p,mx,w,p2,p1,num,den,mu1,e2,mu2,e3,mu3,e4,mu4,
      yhat2,re,ure,R2,tmpa,tmpb,tmpc,hdep,hconst,hindep,
      chiSC,chiFF,chiNO,chiHE,diagres,wtmp;

@ computation of the likelihood, AIC, SC, HQ and Rbar square @

obs = rows(res);
sig2 = res'res/obs;
sige = res'res/(obs-nreg);   @ variance estimate with dof @
sige = sqrt(sige);           @ standard error of regression @

LL = -(obs/2)*(1 + ln(2*pi*sig2));  @ log-likelihood @

AIC = LL - nreg;                   @ Akaike @
SBC = LL - (nreg/2)*ln(obs);       @ Schwarz @
HQ  = LL - nreg*ln(ln(obs));       @ Hannan-Quinn @

ysq2 = sumc((dep-meanc(dep))^2)/obs;         @ sum of squared dmeaned dep @
RSQ = 1 - (sig2/ysq2);                       @ R-squared @

RBARSQ = 1 - ((obs-1)/(obs-nreg))*(1 - RSQ); @ adjyusted R2 @

@ computation of diagnostic test statistics using the OLS residuals @

e = res;      @ Tx1 residual @
x = reg;      @ Txk indep vble @
p = 1;        @ lags used in SC test @

@ LM test for serial correlation @
@ computing W matrix for SC @

w = zeros(obs,p);   @ Txp @
p2 = 1;
do until p2 > p;

   w[obs,p+1-p2] = e[obs-p+p2-1,1];

   p1 = 1;
   do until p1 > obs-1;

      if obs-p+p2-1-p1 > 0;
         w[obs-p1,p+1-p2] = e[obs-p+p2-1-p1,1];
      endif;
      p1 = p1 + 1;

   endo;
   p2 = p2 + 1;

endo;

wtmp = w'w - w'x*invpd(x'x)*x'w;
num = e'w*invpd(wtmp)*w'e;
den = e'e/obs;
chiSC = num/den;

@ LM test for functional form @

w = yhat.^2;
wtmp = w'w - w'x*invpd(x'x)*x'w;
num = e'w*invpd(wtmp)*w'e;
chiFF = num/den;

@ LM test for normality @

mu1 = sumc(e)/obs;
e2 = e.^2;
mu2 = sumc(e2)/obs;
e3 = e.^3;
mu3 = sumc(e3)/obs;
e4 = e.^4;
mu4 = sumc(e4)/obs;

tmpa = mu3^2/(6*mu2^3);
tmpb = (mu4/mu2^2 - 3)^2;
tmpc = 3*mu1^2/(2*mu2) - (mu3*mu1/mu2^2);

chiNO = obs*(tmpa + (1/24)*tmpb) + obs*tmpc;

@ LM test for heterockedasticity @

hdep   = e2;
hconst = ones(obs,1);
yhat2 = yhat.^2;
hindep = hconst~yhat2;

ure = hdep - hindep*invpd(hindep'hindep)*hindep'hdep; @ unrestricted res @
re = hdep - hconst*invpd(hconst'hconst)*hconst'hdep;  @ restricted res @

R2 = (re're - ure'ure)/re're;
chiHE = obs*R2;

DIAGRES = RSQ|RBARSQ|SIGE|LL|AIC|SBC|HQ|chiSC|chiFF|chiNO|chiHE; @ 11x1 @

RETP(DIAGRES);ENDP;

/*  ********************************************* */
    PROC (4) = LAGSEL(lag_op,Y,X,Z,n,k,T,rel_op);
/*  *********************************************
*   Selection of the lag orders used in the dynamic regressors
*
*   PLAG: lag order of the dependent variable:  Nx1
*   QLAG: lag orders of the independent variables, X: Nxk
*   ***************************************************** */

Local plag,qlag,ms_op,maxlag,maxtmp,mtmp;

IF     lag_op eq 1;  ms_op = 0;
       { plag, qlag } = FIXLAG(n,t,k,rel_op,z);
       mtmp = plag|vec(qlag);   maxlag = maxc(mtmp);

ELSEIF lag_op eq 2;

@ ******** @
  SELBACK:
@ ******** @

   cls;
   "The program allows you to automatically select the lag orders
using any of the following three criteria:

Type 1 if you want to use the Akaike Information Criterion (AIC).
Type 2 if you want to use the Schwarz Baysian Criterion (SBC).
Type 3 if you want to use the Hannan-Quinn Criterion (HQ).";

@* ********************************* @
@* OPTION for lag selection criteria @
@* ********************************* @

  ms_op = con(1,1);

  IF     ms_op eq 1;
  ELSEIF ms_op eq 2;
  ELSEIF ms_op eq 3;
  ELSE; "WARNING! Try again.";
        "Press any key to continue."; wait; cls;
        goto SELBACK;
  ENDIF;

@ ******* @
  LBACK1:
@ ******* @

  "Type the maximum lag order. Notice that your choice must satisfy
the condition that the number of all parameters is less than
(minmum of T(i) - maxlag), that is,

(maxlag+1)*k + maxlag + number of Z < minmum of T(i) - maxlag.";

  maxlag = con(1,1);"";
  maxtmp = (maxlag+1)*k + maxlag + cols(z);  @ # of all parameters @
  if maxtmp < minc(t) - maxlag;
  else; "WARNING! Try again.";
        "Press any key to continue."; wait; cls;
        goto LBACK1;
  endif;

  { plag, qlag } = MSLAG(Y,X,Z,n,T,k,maxlag,ms_op);

ENDIF;

RETP(PLAG,QLAG,ms_op,maxlag);ENDP;

/* ********************************** */
   Proc (2) = FIXLAG(n,t,k,rel_op,z);
/* ********************************** */

Local fix_op,p,q,plag,qlag,qi,pli,samelag,sl,maxtmp,maxlag,sumlag,chklag;

@ ******* @
  FBACK1:
@ ******* @

output off; cls;

"Type 1 if you wish to impose the same lag order for all the variables
(dependent and independent). Otherwise type 2.";"";
samelag = con(1,1);

if     samelag eq 1;
elseif samelag eq 2;
else; "WARNING! Try again.";
      "Press any key to continue"; wait; cls;
      goto FBACK1;
endif;

if rel_op eq 2;
   "Notice that the lag order of the relative variable located in the last
column of X is set to 0 automatically by the program.";"";
endif;

"Notice that it must be that";"";
"[lag(y)+lag(x1)+...+lag(xk)] < minmum of T(i) - (maxlag+1), i = 1,...,N.";"";

if samelag eq 1;

@ ******* @
  SBACK1:
@ ******* @

   "Type the common lag order.";
   sl = con(1,1);

   maxtmp = (sl+1)*k + sl + cols(z);  @ # of all parameters @
   if maxtmp < minc(t) - sl;
   else; "WARNING! Too many parameters. Try again.";  goto SBACK1;
   endif;

   plag = sl.*ones(n,1);
   qlag = sl.*ones(n,k);

   if rel_op eq 2; qlag[.,k] = zeros(n,1);
   endif;

else;

@ ********* @
  FIXBACK1:
@ ********* @

   "Type the lag order for the dependent variable.";
   p = con(1,1);""; plag = p.*ones(n,1);

   format /rd 3,0;
   qlag = zeros(n,k);
   qi = 1;
   do until qi > k;

      if qi eq 1;
         "Type the lag order of the 1st regressor(X).";
      elseif qi eq 2;
         "Type the lag order of the 2nd regressor(X).";
      elseif qi eq 3;
         "Type the lag order of the 3rd regressor(X).";
      else;
         "Type the lag order of the" qi "th regressor(X).";
      endif;

      if qi eq k; if rel_op eq 2; q = 0;
                  else;           q = con(1,1);
                  endif;
                  @ The lag order for the relative variable in X
                    (located in the final column of X) MUST be 0. @
      else;  q = con(1,1);
      endif;

      qlag[.,qi] = q.*ones(n,1);
      qi = qi + 1;
   endo;

   chklag = plag[1]~qlag[1,.]; chklag = chklag';
   maxlag = maxc(chklag);
   sumlag = sumc(chklag);

   maxtmp = sumlag + k + cols(z);  @ # of all parameters @
   if maxtmp < minc(t) - maxlag;
   else; "WARNING! Too many parameters. Try again.";  goto FIXBACK1;
   endif;

endif;

RETP(plag,qlag);ENDP;

/* ******************************************* */
   Proc (2) = MSLAG(Y,X,Z,n,T,k,maxlag,ms_op);
/* ******************************************* */

Local ytemp,xtemp,ini,fini,mxi,tis,mdep,mreg,ztemp,
      y0tmp,y1tmp,y2tmp,y3tmp,y4tmp,chkcnt,
      x0tmp,x1tmp,x2tmp,x3tmp,x4tmp,
      preg,q1reg,q2reg,q3reg,q4reg,q5reg,
      q6reg,q7reg,q8reg,q9reg,q10reg,q6,q7,q8,q9,q10,
      p,q1,q2,q3,q4,q5,aic,sbc,hq,mscnt,msdep,msreg,PLAG,QLAG,
      oldaic,oldsbc,oldhq,aiclag,sbclag,hqlag,taiclag,tsbclag,thqlag;

TAICLAG = zeros(n,k+1);  @ Total lags by AIC @
TSBCLAG = zeros(n,k+1);  @ Total lags by SBC @
THQLAG  = zeros(n,k+1);  @ Total lags by HQ @

ini = 1;   @ initial position for ith group @
fini = 0;  @ end position for ith group @
mxi = 1;
do until mxi > N;


   Tis = T[mxi];             @ # of time periods for each group @
   fini = fini + Tis;

   ytemp = y[ini:fini];      @ y for ith group @
   xtemp = x[ini:fini,.];    @ x for ith group @
   ztemp = z[ini:fini,.];    @ z for ith group @
   ini = ini + Tis;

@ Generation y(-j), x1(-j), j = 1, ..., J, upto J = 4 @

   y0tmp  = ytemp;              @ y @
   y1tmp  = lag(y0tmp);       @ y(-1) @
   y2tmp  = lag(y1tmp);       @ y(-2) @
   y3tmp  = lag(y2tmp);       @ y(-3) @
   y4tmp  = lag(y3tmp);       @ y(-4) @

   x0tmp  = xtemp;          @ x @
   x1tmp  = lag(x0tmp);     @ x(-1) @
   x2tmp  = lag(x1tmp);     @ x(-2) @
   x3tmp  = lag(x2tmp);     @ x(-3) @
   x4tmp  = lag(x3tmp);     @ x(-4) @

   AICLAG = zeros(1,k+1);
   SBCLAG = zeros(1,k+1);
   HQLAG  = zeros(1,k+1);

   mscnt = 0;  @ number of total computations @

   mdep = y0tmp;
   msdep = mdep[maxlag+1:rows(mdep)];   @ dep vble used in MS regression @

   if     k eq 1; goto K1;
   elseif k eq 2; goto K2;
   elseif k eq 3; goto K3;
   elseif k eq 4; goto K4;
   elseif k eq 5; goto K5;
   elseif k eq 6; goto K6;
   elseif k eq 7; goto K7;
   elseif k eq 8; goto K8;
   elseif k eq 9; goto K9;
   elseif k eq 10; goto K10;
   endif;

K1:   @ Case of 1 regressor @

   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         mscnt = mscnt + 1;   @ number of computations @

         if p eq 0;
               mreg = ztemp~q1reg;
         else;
               mreg = ztemp~preg~q1reg;
         endif;
         msreg = mreg[maxlag+1:rows(mreg),.];

         { AIC, SBC, HQ } = MSC(msdep,msreg);

         if mscnt eq 1; OLDAIC = AIC; OLDSBC = SBC; OLDHQ = HQ;
                        AICLAG = p~q1;
                        SBCLAG = p~q1;
                        HQLAG  = p~q1;

         else; if AIC > OLDAIC; AICLAG = p~q1;
                        OLDAIC = AIC;
               endif;
                  if SBC > OLDSBC; SBCLAG = p~q1;
                           OLDSBC = SBC;
                  endif;
                  if HQ > OLDHQ; HQLAG = p~q1;
                          OLDHQ = HQ;
                  endif;
         endif;

         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

K2:   @ Case of 2 regressors @

   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         q2 = 0;
         do until q2 > maxlag;

            if     q2 eq 0; q2reg = x0tmp[.,2];
            elseif q2 eq 1; q2reg = x0tmp[.,2]~x1tmp[.,2];
            elseif q2 eq 2; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2];
            elseif q2 eq 3; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2];
            elseif q2 eq 4; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2]~x4tmp[.,2];
            endif;

            mscnt = mscnt + 1;   @ number of computations @

            if p eq 0;
               mreg = ztemp~q1reg~q2reg;
            else;
               mreg = ztemp~preg~q1reg~q2reg;
            endif;
            msreg = mreg[maxlag+1:rows(mreg),.];

            { AIC, SBC, HQ } = MSC(msdep,msreg);

            if mscnt eq 1; OLDAIC = AIC; OLDSBC = SBC; OLDHQ = HQ;
                           AICLAG = p~q1~q2;
                           SBCLAG = p~q1~q2;
                           HQLAG  = p~q1~q2;

            else; if AIC > OLDAIC; AICLAG = p~q1~q2;
                           OLDAIC = AIC;
                  endif;
                  if SBC > OLDSBC; SBCLAG = p~q1~q2;
                           OLDSBC = SBC;
                  endif;
                  if HQ > OLDHQ; HQLAG = p~q1~q2;
                          OLDHQ = HQ;
                  endif;
            endif;

            q2 = q2 + 1;
         endo;
         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

K3:   @ Case of 3 regressors @

   chkcnt = 0;
   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         q2 = 0;
         do until q2 > maxlag;

            if     q2 eq 0; q2reg = x0tmp[.,2];
            elseif q2 eq 1; q2reg = x0tmp[.,2]~x1tmp[.,2];
            elseif q2 eq 2; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2];
            elseif q2 eq 3; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2];
            elseif q2 eq 4; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                 ~x3tmp[.,2]~x4tmp[.,2];
            endif;

            q3 = 0;
            do until q3 > maxlag;

               if     q3 eq 0; q3reg = x0tmp[.,3];
               elseif q3 eq 1; q3reg = x0tmp[.,3]~x1tmp[.,3];
               elseif q3 eq 2; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3];
               elseif q3 eq 3; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3];
               elseif q3 eq 4; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3]~x4tmp[.,3];
               endif;

               mscnt = mscnt + 1;   @ number of computations @

               if p eq 0;
                  mreg = ztemp~q1reg~q2reg~q3reg;
               else;
                  mreg = ztemp~preg~q1reg~q2reg~q3reg;
               endif;
               msreg = mreg[maxlag+1:rows(mreg),.];

               { AIC, SBC, HQ } = MSC(msdep,msreg);

               if mscnt eq 1; OLDAIC = AIC; OLDSBC = SBC; OLDHQ = HQ;
                              AICLAG = p~q1~q2~q3;
                              SBCLAG = p~q1~q2~q3;
                              HQLAG  = p~q1~q2~q3;

               else; if AIC > OLDAIC; AICLAG = p~q1~q2~q3;
                              OLDAIC = AIC;
                     endif;
                     if SBC > OLDSBC; SBCLAG = p~q1~q2~q3;
                              OLDSBC = SBC;
                     endif;
                     if HQ > OLDHQ; HQLAG = p~q1~q2~q3;
                             OLDHQ = HQ;
                     endif;
               endif;

               q3 = q3 + 1;
            endo;
            q2 = q2 + 1;
         endo;
         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

K4:   @ Case of 4 regressors @

   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         q2 = 0;
         do until q2 > maxlag;

            if     q2 eq 0; q2reg = x0tmp[.,2];
            elseif q2 eq 1; q2reg = x0tmp[.,2]~x1tmp[.,2];
            elseif q2 eq 2; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2];
            elseif q2 eq 3; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2];
            elseif q2 eq 4; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                 ~x3tmp[.,2]~x4tmp[.,2];
            endif;

            q3 = 0;
            do until q3 > maxlag;

               if     q3 eq 0; q3reg = x0tmp[.,3];
               elseif q3 eq 1; q3reg = x0tmp[.,3]~x1tmp[.,3];
               elseif q3 eq 2; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3];
               elseif q3 eq 3; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3];
               elseif q3 eq 4; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3]~x4tmp[.,3];
               endif;

               q4 = 0;
               do until q4 > maxlag;

                  if     q4 eq 0; q4reg = x0tmp[.,4];
                  elseif q4 eq 1; q4reg = x0tmp[.,4]~x1tmp[.,4];
                  elseif q4 eq 2; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4];
                  elseif q4 eq 3; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4];
                  elseif q4 eq 4; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4]~x4tmp[.,4];
                  endif;

                  mscnt = mscnt + 1;   @ number of computations @

                  if p eq 0;
                     mreg = ztemp~q1reg~q2reg~q3reg~q4reg;
                  else;
                     mreg = ztemp~preg~q1reg~q2reg~q3reg~q4reg;
                  endif;
                  msreg = mreg[maxlag+1:rows(mreg),.];

                  { AIC, SBC, HQ } = MSC(msdep,msreg);

                  if mscnt eq 1; OLDAIC = AIC; OLDSBC = SBC; OLDHQ = HQ;
                                 AICLAG = p~q1~q2~q3~q4;
                                 SBCLAG = p~q1~q2~q3~q4;
                                 HQLAG = p~q1~q2~q3~q4;

                  else; if AIC > OLDAIC; AICLAG = p~q1~q2~q3~q4;
                                 OLDAIC = AIC;
                        endif;
                        if SBC > OLDSBC; SBCLAG = p~q1~q2~q3~q4;
                                 OLDSBC = SBC;
                        endif;
                        if HQ > OLDHQ; HQLAG = p~q1~q2~q3~q4;
                                OLDHQ = HQ;
                        endif;
                  endif;

                  q4 = q4 + 1;
               endo;
               q3 = q3 + 1;
            endo;
            q2 = q2 + 1;
         endo;
         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

K5:   @ Case of 5 regressors @

   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         q2 = 0;
         do until q2 > maxlag;

            if     q2 eq 0; q2reg = x0tmp[.,2];
            elseif q2 eq 1; q2reg = x0tmp[.,2]~x1tmp[.,2];
            elseif q2 eq 2; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2];
            elseif q2 eq 3; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2];
            elseif q2 eq 4; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                 ~x3tmp[.,2]~x4tmp[.,2];
            endif;

            q3 = 0;
            do until q3 > maxlag;

               if     q3 eq 0; q3reg = x0tmp[.,3];
               elseif q3 eq 1; q3reg = x0tmp[.,3]~x1tmp[.,3];
               elseif q3 eq 2; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3];
               elseif q3 eq 3; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3];
               elseif q3 eq 4; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3]~x4tmp[.,3];
               endif;

               q4 = 0;
               do until q4 > maxlag;

                  if     q4 eq 0; q4reg = x0tmp[.,4];
                  elseif q4 eq 1; q4reg = x0tmp[.,4]~x1tmp[.,4];
                  elseif q4 eq 2; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4];
                  elseif q4 eq 3; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4];
                  elseif q4 eq 4; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4]~x4tmp[.,4];
                  endif;

                  q5 = 0;
                  do until q5 > maxlag;

                     if     q5 eq 0; q5reg = x0tmp[.,5];
                     elseif q5 eq 1; q5reg = x0tmp[.,5]~x1tmp[.,5];
                     elseif q5 eq 2; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5];
                     elseif q5 eq 3; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5]~x3tmp[.,5];
                     elseif q5 eq 4; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                       ~x2tmp[.,5]~x3tmp[.,5]~x4tmp[.,5];
                     endif;

                     mscnt = mscnt + 1;   @ number of computations @

                     if p eq 0;
                        mreg = ztemp~q1reg~q2reg~q3reg~q4reg~q5reg;
                     else;
                        mreg = ztemp~preg~q1reg~q2reg~q3reg~q4reg~q5reg;
                     endif;
                     msreg = mreg[maxlag+1:rows(mreg),.];

                     { AIC, SBC, HQ } = MSC(msdep,msreg);

                     if mscnt eq 1; OLDAIC = AIC; OLDSBC = SBC; OLDHQ = HQ;
                                    AICLAG = p~q1~q2~q3~q4~q5;
                                    SBCLAG = p~q1~q2~q3~q4~q5;
                                    HQLAG = p~q1~q2~q3~q4~q5;

                     else; if AIC > OLDAIC; AICLAG = p~q1~q2~q3~q4~q5;
                                    OLDAIC = AIC;
                           endif;
                           if SBC > OLDSBC; SBCLAG = p~q1~q2~q3~q4~q5;
                                    OLDSBC = SBC;
                           endif;
                           if HQ > OLDHQ; HQLAG = p~q1~q2~q3~q4~q5;
                                   OLDHQ = HQ;
                           endif;
                     endif;

                     q5 = q5 + 1;
                  endo;
                  q4 = q4 + 1;
               endo;
               q3 = q3 + 1;
            endo;
            q2 = q2 + 1;
         endo;
         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

K6:   @ Case of 6 regressors @

   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         q2 = 0;
         do until q2 > maxlag;

            if     q2 eq 0; q2reg = x0tmp[.,2];
            elseif q2 eq 1; q2reg = x0tmp[.,2]~x1tmp[.,2];
            elseif q2 eq 2; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2];
            elseif q2 eq 3; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2];
            elseif q2 eq 4; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                 ~x3tmp[.,2]~x4tmp[.,2];
            endif;

            q3 = 0;
            do until q3 > maxlag;

               if     q3 eq 0; q3reg = x0tmp[.,3];
               elseif q3 eq 1; q3reg = x0tmp[.,3]~x1tmp[.,3];
               elseif q3 eq 2; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3];
               elseif q3 eq 3; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3];
               elseif q3 eq 4; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3]~x4tmp[.,3];
               endif;

               q4 = 0;
               do until q4 > maxlag;

                  if     q4 eq 0; q4reg = x0tmp[.,4];
                  elseif q4 eq 1; q4reg = x0tmp[.,4]~x1tmp[.,4];
                  elseif q4 eq 2; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4];
                  elseif q4 eq 3; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4];
                  elseif q4 eq 4; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4]~x4tmp[.,4];
                  endif;

                  q5 = 0;
                  do until q5 > maxlag;

                     if     q5 eq 0; q5reg = x0tmp[.,5];
                     elseif q5 eq 1; q5reg = x0tmp[.,5]~x1tmp[.,5];
                     elseif q5 eq 2; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5];
                     elseif q5 eq 3; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5]~x3tmp[.,5];
                     elseif q5 eq 4; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                       ~x2tmp[.,5]~x3tmp[.,5]~x4tmp[.,5];
                     endif;

                     q6 = 0;
                     do until q6 > maxlag;

                        if     q6 eq 0; q6reg = x0tmp[.,6];
                        elseif q6 eq 1; q6reg = x0tmp[.,6]~x1tmp[.,6];
                        elseif q6 eq 2; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6];
                        elseif q6 eq 3; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6]~x3tmp[.,6];
                        elseif q6 eq 4; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                       ~x2tmp[.,6]~x3tmp[.,6]~x4tmp[.,6];
                        endif;

                        mscnt = mscnt + 1;   @ number of computations @

                        if p eq 0;
                           mreg = ztemp~q1reg~q2reg~q3reg~q4reg~q5reg~q6reg;
                        else;
                           mreg = ztemp~preg~q1reg~q2reg~q3reg~q4reg~q5reg
                                 ~q6reg;
                        endif;
                        msreg = mreg[maxlag+1:rows(mreg),.];

                        { AIC, SBC, HQ } = MSC(msdep,msreg);

                        if mscnt eq 1; OLDAIC = AIC; OLDSBC = SBC; OLDHQ = HQ;
                                       AICLAG = p~q1~q2~q3~q4~q5~q6;
                                       SBCLAG = p~q1~q2~q3~q4~q5~q6;
                                       HQLAG = p~q1~q2~q3~q4~q5~q6;

                        else; if AIC > OLDAIC; AICLAG = p~q1~q2~q3~q4~q5~q6;
                                       OLDAIC = AIC;
                              endif;
                              if SBC > OLDSBC; SBCLAG = p~q1~q2~q3~q4~q5~q6;
                                       OLDSBC = SBC;
                              endif;
                              if HQ > OLDHQ; HQLAG = p~q1~q2~q3~q4~q5~q6;
                                      OLDHQ = HQ;
                              endif;
                        endif;

                        q6 = q6 + 1;
                     endo;
                     q5 = q5 + 1;
                  endo;
                  q4 = q4 + 1;
               endo;
               q3 = q3 + 1;
            endo;
            q2 = q2 + 1;
         endo;
         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

K7:   @ Case of 7 regressors @

   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         q2 = 0;
         do until q2 > maxlag;

            if     q2 eq 0; q2reg = x0tmp[.,2];
            elseif q2 eq 1; q2reg = x0tmp[.,2]~x1tmp[.,2];
            elseif q2 eq 2; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2];
            elseif q2 eq 3; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2];
            elseif q2 eq 4; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                 ~x3tmp[.,2]~x4tmp[.,2];
            endif;

            q3 = 0;
            do until q3 > maxlag;

               if     q3 eq 0; q3reg = x0tmp[.,3];
               elseif q3 eq 1; q3reg = x0tmp[.,3]~x1tmp[.,3];
               elseif q3 eq 2; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3];
               elseif q3 eq 3; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3];
               elseif q3 eq 4; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3]~x4tmp[.,3];
               endif;

               q4 = 0;
               do until q4 > maxlag;

                  if     q4 eq 0; q4reg = x0tmp[.,4];
                  elseif q4 eq 1; q4reg = x0tmp[.,4]~x1tmp[.,4];
                  elseif q4 eq 2; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4];
                  elseif q4 eq 3; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4];
                  elseif q4 eq 4; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4]~x4tmp[.,4];
                  endif;

                  q5 = 0;
                  do until q5 > maxlag;

                     if     q5 eq 0; q5reg = x0tmp[.,5];
                     elseif q5 eq 1; q5reg = x0tmp[.,5]~x1tmp[.,5];
                     elseif q5 eq 2; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5];
                     elseif q5 eq 3; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5]~x3tmp[.,5];
                     elseif q5 eq 4; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                       ~x2tmp[.,5]~x3tmp[.,5]~x4tmp[.,5];
                     endif;

                     q6 = 0;
                     do until q6 > maxlag;

                        if     q6 eq 0; q6reg = x0tmp[.,6];
                        elseif q6 eq 1; q6reg = x0tmp[.,6]~x1tmp[.,6];
                        elseif q6 eq 2; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6];
                        elseif q6 eq 3; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6]~x3tmp[.,6];
                        elseif q6 eq 4; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                       ~x2tmp[.,6]~x3tmp[.,6]~x4tmp[.,6];
                        endif;

                        q7 = 0;
                        do until q7 > maxlag;

                           if     q7 eq 0; q7reg = x0tmp[.,7];
                           elseif q7 eq 1; q7reg = x0tmp[.,7]~x1tmp[.,7];
                           elseif q7 eq 2; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                                  ~x2tmp[.,7];
                           elseif q7 eq 3; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                                  ~x2tmp[.,7]~x3tmp[.,7];
                           elseif q7 eq 4; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                       ~x2tmp[.,7]~x3tmp[.,7]~x4tmp[.,7];
                           endif;

                           mscnt = mscnt + 1;   @ number of computations @

                           if p eq 0;
                              mreg = ztemp~q1reg~q2reg~q3reg~q4reg~q5reg
                                    ~q6reg~q7reg;
                           else;
                              mreg = ztemp~preg~q1reg~q2reg~q3reg~q4reg~q5reg
                                    ~q6reg~q7reg;
                           endif;
                           msreg = mreg[maxlag+1:rows(mreg),.];

                           { AIC, SBC, HQ } = MSC(msdep,msreg);

                           if mscnt eq 1; OLDAIC = AIC; OLDSBC = SBC;
                                          OLDHQ = HQ;
                                          AICLAG = p~q1~q2~q3~q4~q5~q6~q7;
                                          SBCLAG = p~q1~q2~q3~q4~q5~q6~q7;
                                          HQLAG = p~q1~q2~q3~q4~q5~q6~q7;

                           else; if AIC > OLDAIC;
                                          AICLAG = p~q1~q2~q3~q4~q5~q6~q7;
                                          OLDAIC = AIC;
                                 endif;
                                 if SBC > OLDSBC;
                                          SBCLAG = p~q1~q2~q3~q4~q5~q6~q7;
                                          OLDSBC = SBC;
                                 endif;
                                 if HQ > OLDHQ;
                                         HQLAG = p~q1~q2~q3~q4~q5~q6~q7;
                                         OLDHQ = HQ;
                                 endif;
                           endif;

                           q7 = q7 + 1;
                        endo;
                        q6 = q6 + 1;
                     endo;
                     q5 = q5 + 1;
                  endo;
                  q4 = q4 + 1;
               endo;
               q3 = q3 + 1;
            endo;
            q2 = q2 + 1;
         endo;
         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

K8:   @ Case of 8 regressors @

   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         q2 = 0;
         do until q2 > maxlag;

            if     q2 eq 0; q2reg = x0tmp[.,2];
            elseif q2 eq 1; q2reg = x0tmp[.,2]~x1tmp[.,2];
            elseif q2 eq 2; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2];
            elseif q2 eq 3; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2];
            elseif q2 eq 4; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                 ~x3tmp[.,2]~x4tmp[.,2];
            endif;

            q3 = 0;
            do until q3 > maxlag;

               if     q3 eq 0; q3reg = x0tmp[.,3];
               elseif q3 eq 1; q3reg = x0tmp[.,3]~x1tmp[.,3];
               elseif q3 eq 2; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3];
               elseif q3 eq 3; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3];
               elseif q3 eq 4; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3]~x4tmp[.,3];
               endif;

               q4 = 0;
               do until q4 > maxlag;

                  if     q4 eq 0; q4reg = x0tmp[.,4];
                  elseif q4 eq 1; q4reg = x0tmp[.,4]~x1tmp[.,4];
                  elseif q4 eq 2; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4];
                  elseif q4 eq 3; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4];
                  elseif q4 eq 4; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4]~x4tmp[.,4];
                  endif;

                  q5 = 0;
                  do until q5 > maxlag;

                     if     q5 eq 0; q5reg = x0tmp[.,5];
                     elseif q5 eq 1; q5reg = x0tmp[.,5]~x1tmp[.,5];
                     elseif q5 eq 2; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5];
                     elseif q5 eq 3; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5]~x3tmp[.,5];
                     elseif q5 eq 4; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                       ~x2tmp[.,5]~x3tmp[.,5]~x4tmp[.,5];
                     endif;

                     q6 = 0;
                     do until q6 > maxlag;

                        if     q6 eq 0; q6reg = x0tmp[.,6];
                        elseif q6 eq 1; q6reg = x0tmp[.,6]~x1tmp[.,6];
                        elseif q6 eq 2; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6];
                        elseif q6 eq 3; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6]~x3tmp[.,6];
                        elseif q6 eq 4; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                       ~x2tmp[.,6]~x3tmp[.,6]~x4tmp[.,6];
                        endif;

                        q7 = 0;
                        do until q7 > maxlag;

                           if     q7 eq 0; q7reg = x0tmp[.,7];
                           elseif q7 eq 1; q7reg = x0tmp[.,7]~x1tmp[.,7];
                           elseif q7 eq 2; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                                  ~x2tmp[.,7];
                           elseif q7 eq 3; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                                  ~x2tmp[.,7]~x3tmp[.,7];
                           elseif q7 eq 4; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                       ~x2tmp[.,7]~x3tmp[.,7]~x4tmp[.,7];
                           endif;

                           q8 = 0;
                           do until q8 > maxlag;

                              if     q8 eq 0; q8reg = x0tmp[.,8];
                              elseif q8 eq 1; q8reg = x0tmp[.,8]~x1tmp[.,8];
                              elseif q8 eq 2; q8reg = x0tmp[.,8]~x1tmp[.,8]
                                                     ~x2tmp[.,8];
                              elseif q8 eq 3; q8reg = x0tmp[.,8]~x1tmp[.,8]
                                                     ~x2tmp[.,8]~x3tmp[.,8];
                              elseif q8 eq 4; q8reg = x0tmp[.,8]~x1tmp[.,8]
                                       ~x2tmp[.,7]~x3tmp[.,8]~x4tmp[.,8];
                              endif;

                              mscnt = mscnt + 1;   @ number of computations @

                              if p eq 0;
                                 mreg = ztemp~q1reg~q2reg~q3reg~q4reg~q5reg
                                       ~q6reg~q7reg~q8reg;
                              else;
                                mreg = ztemp~preg~q1reg~q2reg~q3reg~q4reg~q5reg
                                       ~q6reg~q7reg~q8reg;
                              endif;
                              msreg = mreg[maxlag+1:rows(mreg),.];

                              { AIC, SBC, HQ } = MSC(msdep,msreg);

                              if mscnt eq 1;
                                 OLDAIC = AIC; OLDSBC = SBC; OLDHQ = HQ;
                                 AICLAG = p~q1~q2~q3~q4~q5~q6~q7~q8;
                                 SBCLAG = p~q1~q2~q3~q4~q5~q6~q7~q8;
                                 HQLAG = p~q1~q2~q3~q4~q5~q6~q7~q8;

                              else; if AIC > OLDAIC;
                                       AICLAG = p~q1~q2~q3~q4~q5~q6~q7~q8;
                                       OLDAIC = AIC;
                                    endif;
                                    if SBC > OLDSBC;
                                       SBCLAG = p~q1~q2~q3~q4~q5~q6~q7~q8;
                                       OLDSBC = SBC;
                                    endif;
                                    if HQ > OLDHQ;
                                       HQLAG = p~q1~q2~q3~q4~q5~q6~q7~q8;
                                       OLDHQ = HQ;
                                    endif;
                              endif;

                              q8 = q8 + 1;
                           endo;
                           q7 = q7 + 1;
                        endo;
                        q6 = q6 + 1;
                     endo;
                     q5 = q5 + 1;
                  endo;
                  q4 = q4 + 1;
               endo;
               q3 = q3 + 1;
            endo;
            q2 = q2 + 1;
         endo;
         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

K9:   @ Case of 9 regressors @

   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         q2 = 0;
         do until q2 > maxlag;

            if     q2 eq 0; q2reg = x0tmp[.,2];
            elseif q2 eq 1; q2reg = x0tmp[.,2]~x1tmp[.,2];
            elseif q2 eq 2; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2];
            elseif q2 eq 3; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2];
            elseif q2 eq 4; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                 ~x3tmp[.,2]~x4tmp[.,2];
            endif;

            q3 = 0;
            do until q3 > maxlag;

               if     q3 eq 0; q3reg = x0tmp[.,3];
               elseif q3 eq 1; q3reg = x0tmp[.,3]~x1tmp[.,3];
               elseif q3 eq 2; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3];
               elseif q3 eq 3; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3];
               elseif q3 eq 4; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3]~x4tmp[.,3];
               endif;

               q4 = 0;
               do until q4 > maxlag;

                  if     q4 eq 0; q4reg = x0tmp[.,4];
                  elseif q4 eq 1; q4reg = x0tmp[.,4]~x1tmp[.,4];
                  elseif q4 eq 2; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4];
                  elseif q4 eq 3; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4];
                  elseif q4 eq 4; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4]~x4tmp[.,4];
                  endif;

                  q5 = 0;
                  do until q5 > maxlag;

                     if     q5 eq 0; q5reg = x0tmp[.,5];
                     elseif q5 eq 1; q5reg = x0tmp[.,5]~x1tmp[.,5];
                     elseif q5 eq 2; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5];
                     elseif q5 eq 3; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5]~x3tmp[.,5];
                     elseif q5 eq 4; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                       ~x2tmp[.,5]~x3tmp[.,5]~x4tmp[.,5];
                     endif;

                     q6 = 0;
                     do until q6 > maxlag;

                        if     q6 eq 0; q6reg = x0tmp[.,6];
                        elseif q6 eq 1; q6reg = x0tmp[.,6]~x1tmp[.,6];
                        elseif q6 eq 2; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6];
                        elseif q6 eq 3; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6]~x3tmp[.,6];
                        elseif q6 eq 4; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                       ~x2tmp[.,6]~x3tmp[.,6]~x4tmp[.,6];
                        endif;

                        q7 = 0;
                        do until q7 > maxlag;

                           if     q7 eq 0; q7reg = x0tmp[.,7];
                           elseif q7 eq 1; q7reg = x0tmp[.,7]~x1tmp[.,7];
                           elseif q7 eq 2; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                                  ~x2tmp[.,7];
                           elseif q7 eq 3; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                                  ~x2tmp[.,7]~x3tmp[.,7];
                           elseif q7 eq 4; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                       ~x2tmp[.,7]~x3tmp[.,7]~x4tmp[.,7];
                           endif;

                           q8 = 0;
                           do until q8 > maxlag;

                              if     q8 eq 0; q8reg = x0tmp[.,8];
                              elseif q8 eq 1; q8reg = x0tmp[.,8]~x1tmp[.,8];
                              elseif q8 eq 2; q8reg = x0tmp[.,8]~x1tmp[.,8]
                                                     ~x2tmp[.,8];
                              elseif q8 eq 3; q8reg = x0tmp[.,8]~x1tmp[.,8]
                                                     ~x2tmp[.,8]~x3tmp[.,8];
                              elseif q8 eq 4; q8reg = x0tmp[.,8]~x1tmp[.,8]
                                       ~x2tmp[.,7]~x3tmp[.,8]~x4tmp[.,8];
                              endif;

                              q9 = 0;
                              do until q9 > maxlag;

                                 if     q9 eq 0; q9reg = x0tmp[.,9];
                                 elseif q9 eq 1; q9reg = x0tmp[.,9]~x1tmp[.,9];
                                 elseif q9 eq 2; q9reg = x0tmp[.,9]~x1tmp[.,9]
                                                        ~x2tmp[.,9];
                                 elseif q9 eq 3; q9reg = x0tmp[.,9]~x1tmp[.,9]
                                                        ~x2tmp[.,9]~x3tmp[.,9];
                                 elseif q9 eq 4; q9reg = x0tmp[.,9]~x1tmp[.,9]
                                             ~x2tmp[.,9]~x3tmp[.,9]~x4tmp[.,9];
                                 endif;

                                 mscnt = mscnt + 1; @ number of computations @

                                 if p eq 0;
                                    mreg = ztemp~q1reg~q2reg~q3reg~q4reg~q5reg
                                          ~q6reg~q7reg~q8reg~q9reg;
                                 else;
                                    mreg = ztemp~preg~q1reg~q2reg~q3reg~q4reg
                                          ~q5reg~q6reg~q7reg~q8reg~q9reg;
                                 endif;
                                 msreg = mreg[maxlag+1:rows(mreg),.];

                                 { AIC, SBC, HQ } = MSC(msdep,msreg);

                                 if mscnt eq 1;
                                    OLDAIC = AIC; OLDSBC = SBC; OLDHQ = HQ;
                                    AICLAG = p~q1~q2~q3~q4~q5~q6~q7~q8~q9;
                                    SBCLAG = p~q1~q2~q3~q4~q5~q6~q7~q8~q9;
                                    HQLAG = p~q1~q2~q3~q4~q5~q6~q7~q8~q9;

                                 else; if AIC > OLDAIC;
                                         AICLAG = p~q1~q2~q3~q4~q5~q6~q7~q8~q9;
                                         OLDAIC = AIC;
                                       endif;
                                       if SBC > OLDSBC;
                                         SBCLAG = p~q1~q2~q3~q4~q5~q6~q7~q8~q9;
                                          OLDSBC = SBC;
                                       endif;
                                       if HQ > OLDHQ;
                                          HQLAG = p~q1~q2~q3~q4~q5~q6~q7~q8~q9;
                                          OLDHQ = HQ;
                                       endif;
                                 endif;

                                 q9 = q9 + 1;
                              endo;
                              q8 = q8 + 1;
                           endo;
                           q7 = q7 + 1;
                        endo;
                        q6 = q6 + 1;
                     endo;
                     q5 = q5 + 1;
                  endo;
                  q4 = q4 + 1;
               endo;
               q3 = q3 + 1;
            endo;
            q2 = q2 + 1;
         endo;
         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

K10:  @ Case of 10 regressors @

   p = 0;
   do until p > maxlag;

      if     p eq 0;
      elseif p eq 1; preg = y1tmp;
      elseif p eq 2; preg = y1tmp~y2tmp;
      elseif p eq 3; preg = y1tmp~y2tmp~y3tmp;
      elseif p eq 4; preg = y1tmp~y2tmp~y3tmp~y4tmp;
      endif;

      q1 = 0;
      do until q1 > maxlag;

         if     q1 eq 0; q1reg = x0tmp[.,1];
         elseif q1 eq 1; q1reg = x0tmp[.,1]~x1tmp[.,1];
         elseif q1 eq 2; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1];
         elseif q1 eq 3; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]~x3tmp[.,1];
         elseif q1 eq 4; q1reg = x0tmp[.,1]~x1tmp[.,1]~x2tmp[.,1]
                                ~x3tmp[.,1]~x4tmp[.,1];
         endif;

         q2 = 0;
         do until q2 > maxlag;

            if     q2 eq 0; q2reg = x0tmp[.,2];
            elseif q2 eq 1; q2reg = x0tmp[.,2]~x1tmp[.,2];
            elseif q2 eq 2; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2];
            elseif q2 eq 3; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                   ~x3tmp[.,2];
            elseif q2 eq 4; q2reg = x0tmp[.,2]~x1tmp[.,2]~x2tmp[.,2]
                                 ~x3tmp[.,2]~x4tmp[.,2];
            endif;

            q3 = 0;
            do until q3 > maxlag;

               if     q3 eq 0; q3reg = x0tmp[.,3];
               elseif q3 eq 1; q3reg = x0tmp[.,3]~x1tmp[.,3];
               elseif q3 eq 2; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3];
               elseif q3 eq 3; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3];
               elseif q3 eq 4; q3reg = x0tmp[.,3]~x1tmp[.,3]~x2tmp[.,3]
                                    ~x3tmp[.,3]~x4tmp[.,3];
               endif;

               q4 = 0;
               do until q4 > maxlag;

                  if     q4 eq 0; q4reg = x0tmp[.,4];
                  elseif q4 eq 1; q4reg = x0tmp[.,4]~x1tmp[.,4];
                  elseif q4 eq 2; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4];
                  elseif q4 eq 3; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4];
                  elseif q4 eq 4; q4reg = x0tmp[.,4]~x1tmp[.,4]
                                       ~x2tmp[.,4]~x3tmp[.,4]~x4tmp[.,4];
                  endif;

                  q5 = 0;
                  do until q5 > maxlag;

                     if     q5 eq 0; q5reg = x0tmp[.,5];
                     elseif q5 eq 1; q5reg = x0tmp[.,5]~x1tmp[.,5];
                     elseif q5 eq 2; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5];
                     elseif q5 eq 3; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                          ~x2tmp[.,5]~x3tmp[.,5];
                     elseif q5 eq 4; q5reg = x0tmp[.,5]~x1tmp[.,5]
                                       ~x2tmp[.,5]~x3tmp[.,5]~x4tmp[.,5];
                     endif;

                     q6 = 0;
                     do until q6 > maxlag;

                        if     q6 eq 0; q6reg = x0tmp[.,6];
                        elseif q6 eq 1; q6reg = x0tmp[.,6]~x1tmp[.,6];
                        elseif q6 eq 2; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6];
                        elseif q6 eq 3; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                               ~x2tmp[.,6]~x3tmp[.,6];
                        elseif q6 eq 4; q6reg = x0tmp[.,6]~x1tmp[.,6]
                                       ~x2tmp[.,6]~x3tmp[.,6]~x4tmp[.,6];
                        endif;

                        q7 = 0;
                        do until q7 > maxlag;

                           if     q7 eq 0; q7reg = x0tmp[.,7];
                           elseif q7 eq 1; q7reg = x0tmp[.,7]~x1tmp[.,7];
                           elseif q7 eq 2; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                                  ~x2tmp[.,7];
                           elseif q7 eq 3; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                                  ~x2tmp[.,7]~x3tmp[.,7];
                           elseif q7 eq 4; q7reg = x0tmp[.,7]~x1tmp[.,7]
                                       ~x2tmp[.,7]~x3tmp[.,7]~x4tmp[.,7];
                           endif;

                           q8 = 0;
                           do until q8 > maxlag;

                              if     q8 eq 0; q8reg = x0tmp[.,8];
                              elseif q8 eq 1; q8reg = x0tmp[.,8]~x1tmp[.,8];
                              elseif q8 eq 2; q8reg = x0tmp[.,8]~x1tmp[.,8]
                                                     ~x2tmp[.,8];
                              elseif q8 eq 3; q8reg = x0tmp[.,8]~x1tmp[.,8]
                                                     ~x2tmp[.,8]~x3tmp[.,8];
                              elseif q8 eq 4; q8reg = x0tmp[.,8]~x1tmp[.,8]
                                       ~x2tmp[.,7]~x3tmp[.,8]~x4tmp[.,8];
                              endif;

                              q9 = 0;
                              do until q9 > maxlag;

                                 if     q9 eq 0; q9reg = x0tmp[.,9];
                                 elseif q9 eq 1; q9reg = x0tmp[.,9]~x1tmp[.,9];
                                 elseif q9 eq 2; q9reg = x0tmp[.,9]~x1tmp[.,9]
                                                        ~x2tmp[.,9];
                                 elseif q9 eq 3; q9reg = x0tmp[.,9]~x1tmp[.,9]
                                                        ~x2tmp[.,9]~x3tmp[.,9];
                                 elseif q9 eq 4; q9reg = x0tmp[.,9]~x1tmp[.,9]
                                             ~x2tmp[.,9]~x3tmp[.,9]~x4tmp[.,9];
                                 endif;

                                 q10 = 0;
                                 do until q10> maxlag;

                                    if     q10 eq 0; q10reg = x0tmp[.,10];
                                    elseif q10 eq 1;
                                           q10reg = x0tmp[.,9]~x1tmp[.,10];
                                    elseif q10 eq 2;
                                           q10reg = x0tmp[.,10]~x1tmp[.,10]
                                                   ~x2tmp[.,10];
                                    elseif q10 eq 3;
                                           q10reg = x0tmp[.,10]~x1tmp[.,10]
                                                   ~x2tmp[.,10]~x3tmp[.,10];
                                    elseif q10 eq 4;
                                           q10reg = x0tmp[.,10]~x1tmp[.,10]
                                                   ~x2tmp[.,10]~x3tmp[.,10]
                                                   ~x4tmp[.,10];
                                    endif;

                                    mscnt = mscnt + 1;

                                    if p eq 0;
                                       mreg = ztemp~q1reg~q2reg~q3reg~q4reg
                                        ~q5reg~q6reg~q7reg~q8reg~q9reg~q10reg;
                                    else;
                                      mreg = ztemp~preg~q1reg~q2reg~q3reg~q4reg
                                        ~q5reg~q6reg~q7reg~q8reg~q9reg~q10reg;
                                    endif;
                                    msreg = mreg[maxlag+1:rows(mreg),.];

                                    { AIC, SBC, HQ } = MSC(msdep,msreg);

                                    if mscnt eq 1;
                                       OLDAIC = AIC; OLDSBC = SBC; OLDHQ = HQ;
                                    AICLAG = p~q1~q2~q3~q4~q5~q6~q7~q8~q9~q10;
                                    SBCLAG = p~q1~q2~q3~q4~q5~q6~q7~q8~q9~q10;
                                    HQLAG = p~q1~q2~q3~q4~q5~q6~q7~q8~q9~q10;

                                    else; if AIC > OLDAIC;
                                             AICLAG = p~q1~q2~q3~q4~q5~q6
                                                     ~q7~q8~q9~q10;
                                             OLDAIC = AIC;
                                          endif;
                                          if SBC > OLDSBC;
                                             SBCLAG = p~q1~q2~q3~q4~q5~q6
                                                     ~q7~q8~q9~q10;
                                             OLDSBC = SBC;
                                          endif;
                                          if HQ > OLDHQ;
                                             HQLAG = p~q1~q2~q3~q4~q5~q6
                                                    ~q7~q8~q9~q10;
                                             OLDHQ = HQ;
                                          endif;
                                    endif;

                                    q10 = q10 + 1;
                                 endo;
                                 q9 = q9 + 1;
                              endo;
                              q8 = q8 + 1;
                           endo;
                           q7 = q7 + 1;
                        endo;
                        q6 = q6 + 1;
                     endo;
                     q5 = q5 + 1;
                  endo;
                  q4 = q4 + 1;
               endo;
               q3 = q3 + 1;
            endo;
            q2 = q2 + 1;
         endo;
         q1 = q1 + 1;
      endo;
      p = p + 1;
   endo;
   goto MOUT;

MOUT:

   TAICLAG[mxi,.] = AICLAG;
   TSBCLAG[mxi,.] = SBCLAG;
   THQLAG[mxi,.]  = HQLAG;

   mxi = mxi + 1;
endo;

if     ms_op eq 1;   @ AIC @
       plag = TAICLAG[.,1];
       qlag = TAICLAG[.,2:k+1];

elseif ms_op eq 2;   @ SBC @
       plag = TSBCLAG[.,1];
       qlag = TSBCLAG[.,2:k+1];

elseif ms_op eq 3;   @ HQ @
       plag = THQLAG[.,1];
       qlag = THQLAG[.,2:k+1];
endif;

RETP(PLAG,QLAG);ENDP;

/* **************************** */
   PROC (3) = MSC(msdep,msreg);
/* **************************** */

Local nreg,mobs,mtmp,mpara,mres,msig2,mll,maic,msbc,mhq;

@ computation of the likelihood, AIC, SC, HQ and Rbar square @

nreg = cols(msreg);
mobs = rows(msdep);
mtmp = invpd(msreg'msreg);
mpara = mtmp*msreg'msdep;
mres = msdep - msreg*mpara;
msig2 = mres'mres/mobs;

mLL = -(mobs/2)*(1 + ln(2*pi*msig2));  @ log-likelihood @

mAIC = mLL - nreg;                    @ Akaike @
mSBC = mLL - (nreg/2)*ln(mobs);       @ Schwarz @
mHQ  = mLL - nreg*ln(ln(mobs));       @ Hannan-Quinn @

RETP(mAIC,MSBC,mHQ);ENDP;

/*  ********************************************* */
    PROC (8) = MLEOP(X,Xname,k,sumt,QLAG,rel_op);
/*  ********************************************* */

Local k1,k2,xpos,x1pos,x2pos,ki,k1tmp,k2tmp,x1,x2,x1name,x2name,q1lag,q2lag,
      junk,cnt,kj,tmp,x1pj;

xpos = seqa(1,1,k);    @ order of the stochastic regressors @

@ ******** @
  MLBACK2:
@ ******** @

output off; cls;

"Type the number of regressors (X) with identical long-run coeffcients
across the groups.";

k1 = con(1,1);"";
if k1 eq k; "WARNING! Try again.";
            "Press any key to continue."; wait; cls;
            goto MLBACK2;
endif;
k2 = k - k1;"";

"Specify the numeric position(s) of the regressors whose long-run coefficients
are assumed to be identical across the groups.

For example, if they are included in the X matrix in its 1st and 3rd columns
then type 1, space, and then 3. Then press the ENTER key.";"";

if     rel_op eq 1; x1pos = con(k1,1);
elseif rel_op eq 2;

X1POSBAK:

   "Note: When including a relative variable, its long-run coefficent MUST
be restricted to be the same across groups. Otherwise the coefficient
is not identified.

You only need to specify " k1-1 " such numeric position(s) since the program
automatically imposes the restriction on the relative variable located in
the last column of the X matrix.";

   if k1 eq 1; x1pos = k; "Press any key to continue.";""; wait;
   else;       x1pos = con(k1-1,1);
               x1pj = 1;
               do until x1pj > k1 -1;

                  if x1pos[x1pj] eq k;
                     "WARNING! The relative variable is included automatically
by the program, so do not include it yourself. Try again.";
                     "Press any key to continue."; wait; cls;
                     goto X1POSBAK;
                  endif;
                  x1pj = x1pj + 1;
               endo;
               x1pos = x1pos|k;  @ including the relative variable @
   endif;
endif;

x1pos = sortc(x1pos,1);

junk = xpos;
ki = 1;
do until ki > k1;

   tmp = x1pos[ki];
   cnt = 0;
   kj = 1;
   do until kj > rows(junk);

      if tmp eq junk[kj];
      else;
             if cnt eq 0; x2pos = junk[kj];
             else;        x2pos = x2pos|junk[kj];
             endif;
             cnt = cnt + 1;
      endif;

      kj = kj + 1;
   endo;

   junk = x2pos;
   ki = ki + 1;
endo;
x2pos = sortc(x2pos,1);

ki = 1;
do until ki > k1;

   k1tmp = x1pos[ki];
   if ki eq 1; X1 = X[.,k1tmp];
               X1name = Xname[k1tmp];
               Q1LAG = QLAG[.,k1tmp];
   else;       X1 = X1~X[.,k1tmp];
               X1name = X1name|Xname[k1tmp];
               Q1LAG = Q1LAG~QLAG[.,k1tmp];
   endif;
   ki = ki + 1;
endo;

X2 = zeros(sumT,1);
X2name = "";
ki = 1;
do until ki > k2;

   k2tmp = x2pos[ki];
   if ki eq 1; X2 = X[.,k2tmp];
               X2name = Xname[k2tmp];
               Q2LAG = QLAG[.,k2tmp];
   else;       X2 = X2~X[.,k2tmp];
               X2name = X2name|Xname[k2tmp];
               Q2LAG = Q2LAG~QLAG[.,k2tmp];
   endif;
   ki = ki + 1;
endo;

RETP(X1,X2,X1name,X2name,k1,k2,Q1LAG,Q2LAG);ENDP;

/*  ***********************  */
    PROC (2) = AVG(tmpvec);
/*  ****************************************************************
    Compute the average, its standard error in a non-parameteric way
    **************************************************************** */

Local dimvec,mvec,devvec,sumvec,sevec,tvec,mresult;

dimvec = rows(tmpvec);                        @ dimension of the vector @
mvec = meanc(tmpvec);                         @ mean of vector @
devvec = tmpvec - mvec.*ones(dimvec,1);       @ demeaned vector @
sumvec = sumc(devvec.^2)/(dimvec*(dimvec-1)); @ NP variance of vector @
sevec = sqrt(sumvec);                         @ NP SE of vector @
tvec = mvec/sevec;                            @ t-ratio @

mresult = mvec~sevec~tvec;

RETP(mresult,devvec);ENDP;

