@ ******************* @
@ Printing procedures @
@ ******************* @

/* ********************************************************************* */
   PROC (1) = PRSUM(cflag,iter,sub_op,n_pos,dem_op,ml_op,ini_op,
       RSLL,URSLL,n_LR,n,T,ms_op,Alag,COMPRES,maxlag,nu_op,theta,sm_op);
/* ********************************************************************* */

Local LR,i,tobs,maxalag;

"          **********************************************************";
"          Alternative Estimates for Dynamic Heterogenous Panel Model";
"          **********************************************************";"";

format /rd 6,0;
if cflag eq 1;
   "";"Computations converged after" iter-1 "iterations.";
else;
   "Computations did not converge. Increase the number of iterations.";
endif;"";

format /rd 4,0;
if sub_op eq 1;
@   "All groups are included.";     @
   "The number of groups is N = ";;n;
else;
   "The following groups are excluded (by position in the data set):";""; n_pos'
;"";
   "The number of included groups is N = ";;n;
endif;"";

"The number of whole time periods by groups are:";""; t';"";

if     ms_op eq 0;
   "Orders of fixed lags in the ARDL model which you have specified:";
elseif ms_op eq 1;
   "Orders of lags in the ARDL model which are selected by AIC:";
elseif ms_op eq 2;
   "Orders of lags in the ARDL model which are selected by SBC:";
elseif ms_op eq 3;
   "Orders of lags in the ARDL model which are selected by HQ:";
endif;"";

format /rd 3,0;
if ms_op eq 0;  Alag[1,.];

   if maxlag eq 0; tobs = T - 1;
   else;           tobs = T - maxlag;
   endif;

else;
   tobs = zeros(n,1);
   i=1;
   do until i > n;
      "Group";;i;;": ";;Alag[i,.];
      maxalag = maxc(alag[i,.]');
      if maxalag eq 0; tobs[i] = T[i] - 1;
      else;            tobs[i] = T[i] - maxalag;
      endif;
      i = i + 1;
   endo;
endif;"";

format /rd 4,0;
"The number of time periods by groups used in estimation are:";""; tobs';"";

@ IF     dem_op eq 1;  "The raw data is used.";   @
@ elseif dem_op eq 2;  "The demeaned data is used.";   @
@ endif;"";   @

@IF     ml_op eq 1;  @
@"All the long-run parameters are restricted to be the same across groups.";  @
@elseif ml_op eq 2;  @
@"A subset of the long-run parameters is restricted to be the same across groups
.  @
@";  @
@endif;"";  @

@ format /rd 6,3;  @
@ if     ini_op eq 1;  @
@   "Mean group estimates are used as initial values.";  @
@ elseif ini_op eq 2;  @
@   "Static fixed effects estimates are used as initial values.";  @
@ elseif ini_op eq 3;  @
@   "Dynamic fixed effects estimates are used as initial values.";  @
@ elseif ini_op eq 4;  @
@   "** Your own values are used as initial values.";  @
@ endif;"";  @

@ if     nu_op eq 1;  @
@    "Back-Substitution algorithm is used.";  @
@ elseif nu_op eq 2;  @
@    "Newton-Raphson algorithm is used.";  @
@ endif;"";  @

format /rd 15,4;
"Restricted log likelihood:  ";;RSLL;
"Unrestricted log likelihood:";;URSLL;"";
LR = 2*(URSLL - RSLL);
"LR statistic testing for equal long-run parameters:";;LR;

format /rd 4,0;
"LR statistic is distributed as Chi-Squared with";;(n-1)*n_LR;;
format /rd 7,4;
"degrees of freedom,";
"and its p-value is";;cdfchic(LR,(n-1)*n_LR);"";

"       --------------------------------------------------";
"       Alternative Initial Estimates of the LR Parameters";
"       --------------------------------------------------";

if ms_op eq 0;

"         Mean Group Estimates    Static Fixed Effects    Dynamic Fixed Effects"
;
"         --------------------    --------------------    ---------------------"
;

"        Coef.  St. Er. t-ratio  Coef.  St. Er. t-ratio  Coef. St. Er. t-ratio";
"";
call printfm(compres,0~1~1~1~1~1~1~1~1~1,
             "s"~4~4|"lf"~8~2|"lf"~8~2|"lf"~8~2|"lf"~8~2|"lf"~8~2
             |"lf"~8~2|"lf"~8~2|"lf"~8~2|"lf"~8~2);"";

else;

"          Mean Group Estimates    Static Fixed Effects";
"          --------------------    --------------------";
"        Coef.  St. Er. t-ratio  Coef.  St. Er. t-ratio";"";

call printfm(compres,0~1~1~1~1~1~1,
             "s"~4~4|"lf"~8~2|"lf"~8~2|"lf"~8~2|"lf"~8~2|"lf"~8~2|"lf"~8~2);"";
endif;"";

if sm_op eq 2;
   format /rd 12,6;
   "PMG Estimates of the long-run coefficients are (without standard errors):";
   ""; theta;"";
endif;

RETP(Tobs);ENDP;

/* ************************************************************************** */
   PROC (2) = PRGROUP1(OLSRES,PMLERES,CFNAME,OLSDRES,PMLEDRES,maxlag,
        k,n,T,Tobs,numsi,sub_op,n_gr,grname,na_gr,gn_op1,gn_op2,yname,page);
/* ************************************************************************** */

@ ML ESTIMATION OPTION 1
* ----------------------
* OLSRES obtained for each group by OLS
* ------
* Result in the order of theta1, theta2, phi, beta1, beta2, other para
* for group 1, 2, ..., N:
* # of parameters in each group is k+1+k+numsi(i)
*
* PMLERES
* -------
* Allpara parameters in the order of theta1, theta2(1),...theta2(N),
* phi(1), ..., phi(N), and si(1), ..., si(N):
* The dimension is k+n+k+sumnumsi-vector @

Local ini,fini,pini,pfini,bini,bfini,oini,ofini,ni,otmp,otemp,ctemp,nx,
      temp1,temp2,temp3,temp4,ptemp,prout,drtmp1,drtmp2,drout1,drout2,drout3,
      PMLESUM;

ini  = 1;                fini  = 0;
Pini = k + 1;            Pfini = k;
Bini  = k + n + 1;       Bfini = k + n;
Oini  = k + n + k*n + 1; Ofini = k + n + k*n;

ni = 1;
do until ni > n;

   ""; page = page + 1; format /rd 3,0;
"PAGE "  page;"";

      otmp  = k+1+k+numsi[ni];
      fini  = fini + otmp;
      otemp = OLSRES[ini:fini,.];
      ctemp = CFname[ini:fini,.];
      ini   = ini + otmp;

      temp1 = PMLERES[1:k,.];

      Pfini = Pfini + 1;
      temp2 = PMLERES[Pini:Pfini,.];
      Pini  = Pini + 1;

      Bfini = Bfini + k;
      temp3 = PMLERES[Bini:Bfini,.];
      Bini  = Bini + k;

      Ofini = Ofini + numsi[ni];
      temp4 = PMLERES[Oini:Ofini,.];
      Oini  = Oini + numsi[ni];

      ptemp = temp1|temp2|temp3|temp4;

@ Matching the form of PMLE estimation results with that of OLSRES @

      if ni eq 1; PMLESUM = ptemp;
      else;       PMLESUM = PMLESUM|ptemp;
      endif;

      PRout = ctemp~ptemp~otemp;
      PRHEAD1(k,Tobs,sub_op,n_gr,grname,na_gr,gn_op1,gn_op2,yname,prout,ni);

      drtmp1 = pmledres[.,ni];
      drtmp2 = olsdres[.,ni];

      PRDIAG(drtmp1,drtmp2);

   ni = ni + 1;
endo;

RETP(page,PMLESUM);ENDP;

/* ************************************************************************ */
   PROC (2) = PRGROUP2(OLSRES,PMLERES,CFNAME,OLSDRES,PMLEDRES,maxlag,k1,k2,
        n,T,Tobs,numsi,sub_op,n_gr,grname,na_gr,gn_op1,gn_op2,yname,page);
/* ************************************************************************ */

@ ML ESTIMATION OPTION 2
* ----------------------
* OLSRES obtained for each group by OLS
* ------
* Result in the order of theta1, theta2, phi, beta1, beta2, other para
* for group 1, 2, ..., N:
* # of parameters in each group is k1+k2+1+k1+k2+numsi(i)
*
* PMLERES
* -------
* Allpara parameters in the order of theta1, theta2(1),...theta2(N),
* phi(1), ..., phi(N), and si(1), ..., si(N):
* The dimension is k1+k2*n+n+k1+k2*n+sumnumsi-vector @

Local ini,fini,pini,pfini,b1ini,b1fini,oini,ofini,ni,otmp,otemp,ctemp,
      temp1,temp2,temp3,temp4,ptemp,prout,drtmp1,drtmp2,drout1,drout2,
      temp5,temp6,b2ini,b2fini,t2ini,t2fini,nx,pmlesum;

ini  = 1;                                fini  = 0;
T2ini = k1 + 1;                          T2fini = k1;
Pini  = k1 + k2*n + 1;                   Pfini  = k1 + k2*n;
B1ini = k1 + k2*n + n + 1;               B1fini = k1 + k2*n + n;
B2ini = k1 + k2*n + n + k1*n + 1;        B2fini = k1 + k2*n + n + k1*n;
Oini  = k1 + k2*n + n + k1*n + k2*n + 1; Ofini  = k1 + k2*n + n + k1*n + k2*n;

ni = 1;
do until ni > n;

   ""; page = page + 1; format /rd 3,0;
"PAGE "  page;"";

      otmp = k1+k2+1+k1+k2+numsi[ni];
      fini = fini + otmp;
      otemp = OLSRES[ini:fini,.];
      ctemp = CFname[ini:fini];
      ini = ini + otmp;

      temp1 = PMLERES[1:k1,.];

      T2fini = T2fini + k2;
      temp2  = PMLERES[T2ini:T2fini,.];
      T2ini  = T2ini + k2;

      Pfini = Pfini + 1;
      temp3 = PMLERES[Pini:Pfini,.];
      Pini  = Pini + 1;

      B1fini = B1fini + k1;
      temp4  = PMLERES[B1ini:B1fini,.];
      B1ini  = B1ini + k1;

      B2fini = B2fini + k2;
      temp5  = PMLERES[B2ini:B2fini,.];
      B2ini  = B2ini + k2;

      Ofini = Ofini + numsi[ni];
      temp6 = PMLERES[Oini:Ofini,.];
      Oini  = Oini + numsi[ni];

      ptemp = temp1|temp2|temp3|temp4|temp5|temp6;

@ Matching the form of PMLE estimation results with that of OLSRES @

      if ni eq 1; PMLESUM = ptemp;
      else;       PMLESUM = PMLESUM|ptemp;
      endif;

      prout = ctemp~ptemp~otemp;
      PRHEAD2(k1,k2,Tobs,sub_op,n_gr,grname,na_gr,gn_op1,gn_op2,
              yname,prout,ni);

      drtmp1 = pmledres[.,ni];
      drtmp2 = olsdres[.,ni];

      PRDIAG(drtmp1,drtmp2);

   ni = ni + 1;
endo;

RETP(page,PMLESUM);ENDP;

/* ************ */
   PROC (0) =
   PRHEAD1(k,Tobs,sub_op,n_gr,grname,na_gr,gn_op1,gn_op2,yname,prout,ni);
/* ********************************************************************** */

@ Print the ML and OLS estimation results for each group for ML option 1 @

Local grout,prj,prtmp,prtmp1,chklag;

format /rd 3,0;
"*********  ***************************     ****************************";
IF     sub_op eq 1;
"Group" ni
"  Pooled Mean Group Estimates       Group-Specific Estimates";

ELSEIF sub_op eq 2;
"Group" n_gr[ni]
"  Pooled Mean Group Estimates       Group-Specific Estimates";
ENDIF;
"*********  ***************************     ****************************";

IF     gn_op1 eq 1;
       "             Dependent variable: ";;
       call printfm(Yname,0,"s"~4~4);
       "                         T = " Tobs[ni];

ELSEIF gn_op1 eq 2;

    if     sub_op eq 1; grout = grname[ni];
    elseif sub_op eq 2; grout = na_gr[ni];
    endif;

    if     gn_op2 eq 1;
           call printfm(grout,1,"lf"~8~0);
    elseif gn_op2 eq 2;
           call printfm(grout,0,"s"~8~8);
    endif;

    "    Dependent variable: ";;
    call printfm(Yname,0,"s"~4~4);
       "                         T = " Tobs[ni];

ENDIF;"";

"            Coef.     St. Er.  t-ratio      Coef.     St. Er.  t-ratio";"";

prtmp = prout[1:k,.];
"Long-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
            |"lf"~12~4|"lf"~10~4|"lf"~10~4);
if rows(prtmp) eq 1; "";
endif;"";

prtmp = prout[k+1,.];
"Error Correction Coefficients";"";
if prtmp[2] eq -1;  @ if phi = -1 @
   prtmp1 = prtmp[1]~prtmp[2]~"NA"~"NA"
           ~prtmp[5]~prtmp[6]~prtmp[7];
   call printfm(PRtmp1,0~1~0~0~1~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10
               |"lf"~12~4|"s"~10~10|"s"~10~10);"";
else;
   call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
               |"lf"~12~4|"lf"~10~4|"lf"~10~4);"";
endif;"";

prtmp = prout[k+2:rows(prout),.];
"Short-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
            |"lf"~12~4|"lf"~10~4|"lf"~10~4);
if rows(prtmp) eq 1; "";
endif;"";

RETP;ENDP;

/* **************************************************************** */
   PROC (0) = PRHEAD2
   (k1,k2,Tobs,sub_op,n_gr,grname,na_gr,gn_op1,gn_op2,yname,prout,ni);
/* ****************************************************************** */

@ Print the ML and OLS estimation results for each group for ML option 2 @

Local grout,prj,prtmp,prtmp1;

format /rd 3,0;
"*********  ***************************     ****************************";
IF     sub_op eq 1;
"Group" ni
"  Pooled Mean Group Estimates       Group-Specific Estimates";

ELSEIF sub_op eq 2;
"Group" n_gr[ni]
"  Pooled Mean Group Estimates       Group-Specific Estimates";
ENDIF;
"*********  ***************************     ****************************";

IF     gn_op1 eq 1;
       "             Dependent variable: ";;
       call printfm(Yname,0,"s"~4~4);
       "                         T = " Tobs[ni];

ELSEIF gn_op1 eq 2;

    if     sub_op eq 1; grout = grname[ni];
    elseif sub_op eq 2; grout = na_gr[ni];
    endif;

    if     gn_op2 eq 1;
           call printfm(grout,1,"lf"~8~0);
    elseif gn_op2 eq 2;
           call printfm(grout,0,"s"~8~8);
    endif;

    "    Dependent variable: ";;
    call printfm(Yname,0,"s"~4~4);
       "                         T = " Tobs[ni];

ENDIF;"";

"            Coef.     St. Er.  t-ratio      Coef.     St. Er.  t-ratio";"";

prtmp = prout[1:k1,.];
"Long-Run Coefficients Restricted to be the Same Across All Groups";"";
call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
            |"lf"~12~4|"lf"~10~4|"lf"~10~4);
if rows(prtmp) eq 1; "";
endif;"";

prtmp = prout[k1+1:k1+k2,.];
"Unrestricted Long-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
            |"lf"~12~4|"lf"~10~4|"lf"~10~4);
if rows(prtmp) eq 1; "";
endif;"";

prtmp = prout[k1+k2+1,.];
"Error Correction Coefficients";"";
if prtmp[2] eq -1;  @ if phi = -1 @
   prtmp1 = prtmp[1]~prtmp[2]~"NA"~"NA"
           ~prtmp[5]~prtmp[6]~prtmp[7];
   call printfm(PRtmp1,0~1~0~0~1~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10
               |"lf"~12~4|"s"~10~10|"s"~10~10);"";
else;
   call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
               |"lf"~12~4|"lf"~10~4|"lf"~10~4);"";
endif;"";

prtmp = prout[k1+k2+2:rows(prout),.];
"Short-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
            |"lf"~12~4|"lf"~10~4|"lf"~10~4);
if rows(prtmp) eq 1; "";
endif;"";

RETP;ENDP;

/* ***************************************************************** */
   PROC (0) = PRDIAG(drtmp1,drtmp2);
/* ***************************************************************** */

@ Print the summary statistics and diagnostics @

Local drout1,drout2,drout3;

drout1 = drtmp1[1:3]'~drtmp2[1:3]';
drout2 = drtmp1[4:6]'~drtmp2[4:6]';
drout3 = drtmp1[8:11]'~drtmp2[8:11]';

"                      Summary statistics and diagnostics";
"                      ----------------------------------";"";

"              RSQ     RBARSQ     SIGMA        RSQ     RBARSQ     SIGMA";
call printfm(drout1,1~1~1~1~1~1,
           "lf"~18~3|"lf"~10~3|"lf"~10~3|"lf"~12~3|"lf"~10~3|"lf"~10~3);"";
"";

"               LL       AIC        SC          LL       AIC       SC";
call printfm(drout2,1~1~1~1~1~1,
           "lf"~18~2|"lf"~10~2|"lf"~10~2|"lf"~12~2|"lf"~10~2|"lf"~10~2);"";"";

"         CH-SC   CH-FF   CH-NO   CH-HE       CH-SC   CH-FF   CH-NO   CH-HE";
call printfm(drout3,1~1~1~1~1~1~1~1,
           "lf"~14~2|"lf"~8~2|"lf"~8~2|"lf"~8~2
          |"lf"~12~2|"lf"~8~2|"lf"~8~2|"lf"~8~2);"";"";

RETP;ENDP;

/* ******************************************************************* */
   PROC (0) = PRFIX(Yname,XXname,ZNAME,SFERES,RSFERES,SFEDRES,SFEOBS);
/* *******************************************************************
*  Printing the fixed effects estimation results
*  ********************************************* */

Local pout;

"                  ******************************";
"                  Static Fixed Effects Estimates";
"                  ******************************";
"          Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);
"  Sample size: T(1) + ... + T(N) = " SFEobs;"";
"                            FE                 Robust";"";

"            Coef.     St. Er.  t-ratio    St. Er.  t-ratio";"";

if rows(zname) eq 1;
else; xxname = xxname|zname[1:rows(zname)-1];
endif;

pout = XXname~SFERES~RSFERES[.,2:3];
call printfm(pout,0~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
             |"lf"~10~4|"lf"~10~4);"";"";

"                  Summary statistics and diagnostics";
"                  ----------------------------------";"";
pout = SFEDRES[2:6]';
"            RBARSQ     SIGMA      LL        AIC       SC";"";
call printfm(POUT,1~1~1~1~1,"lf"~18~3|"lf"~10~3|"lf"~10~2|"lf"~10~2|"lf"~10~2);
"";"";

pout = SFEDRES[8:11]';
"             CH-SC     CH-FF     CH-NO     CH-HE";"";
call printfm(POUT,1~1~1~1,"lf"~18~2|"lf"~10~2|"lf"~10~2|"lf"~10~2);"";"";

RETP;ENDP;

/* ****************************************************************** */
   PROC (0) =
   PRFIX1(CFNAME,YNAME,DFERES,RDFERES,DFEDRES,DFEOBS,k,numot,pp);
/* ****************************************************************** */

Local pout,dimind,ctemp,ptemp,prj,prtmp;

"                  *******************************";
"                  Dynamic Fixed Effects Estimates";
"                  *******************************";
"          Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);
"  Sample size: T(1) + ... + T(N) = " DFEobs;"";
"                            FE                 Robust";"";
"            Coef.     St. Er.  t-ratio    St. Er.  t-ratio";"";

pout = DFERES~RDFERES[.,2:3];
if pp eq 0; @ no lagged dep vble @
   ptemp = -1~"NA"~"NA"~"NA"~"NA";
   pout = pout[1:k,.]|ptemp|pout[k+1:rows(pout),.];
endif;
dimind = k + 1 + k + numot;
ctemp = CFname[1:dimind,1];
pout = ctemp~pout;

prtmp = pout[1:k,.];
"Long-Run Coefficients";"";
call printfm(prtmp,0~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
            |"lf"~10~4|"lf"~10~4);"";
if rows(prtmp) eq 1; "";
endif;"";

prtmp = pout[k+1,.];
"Error Correction Coefficients";"";
if prtmp[2] eq -1;  @ if phi = -1 @
   prtmp = prtmp[1]~prtmp[2]~"NA"~"NA"~"NA"~"NA";
   call printfm(prtmp,0~1~0~0~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10
               |"s"~10~10|"s"~10~10);"";
else;
   call printfm(prtmp,0~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
               |"lf"~10~4|"lf"~10~4);"";
endif;"";

prtmp = pout[k+2:rows(pout),.];
"Short-Run Coefficients";"";
call printfm(prtmp,0~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
            |"lf"~10~4|"lf"~10~4);"";
if rows(prtmp) eq 1; "";
endif;"";

"                  Summary statistics and diagnostics";
"                  ----------------------------------";"";
pout = DFEDRES[2:6]';
"            RBARSQ     SIGMA      LL        AIC       SC";"";
call printfm(POUT,1~1~1~1~1,"lf"~18~3|"lf"~10~3|"lf"~10~2|"lf"~10~2|"lf"~10~2);
"";"";

pout = DFEDRES[8:11]';
"             CH-SC     CH-FF     CH-NO     CH-HE";"";
call printfm(POUT,1~1~1~1,"lf"~18~2|"lf"~10~2|"lf"~10~2|"lf"~10~2);"";"";

RETP;ENDP;

/* ******************************************************************* */
   PROC (0) =
   PRFIX2(CFNAME,YNAME,DFERES,RDFERES,DFEDRES,DFEOBS,kk,numot,pp,k1,k2);
/* ******************************************************************** */

Local pout,dimind,ctemp,ptemp,prj,prtmp;

"                  *******************************";
"                  Dynamic Fixed Effects Estimates";
"                  *******************************";
"          Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);
"  Sample size: T(1) + ... + T(N) = " DFEobs;"";
"                            Fe                 Robust";"";
"            Coef.     St. Er.  t-ratio    St. Er.  t-ratio";"";

pout = DFERES~RDFERES[.,2:3];
if pp eq 0; @ no lagged dep vble @
   ptemp = -1~"NA"~"NA"~"NA"~"NA";
   pout = pout[1:kk,.]|ptemp|pout[kk+1:rows(pout),.];
endif;
dimind = kk + 1 + kk + numot;
ctemp = CFname[1:dimind,1];
pout = ctemp~pout;

prtmp = pout[1:k1+k2,.];
"Long-Run Coefficients";"";
call printfm(prtmp,0~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
            |"lf"~10~4|"lf"~10~4);"";
if rows(prtmp) eq 1; "";
endif;"";

prtmp = pout[k1+k2+1,.];
"Error Correction Coefficients";"";
if prtmp[2] eq -1;  @ if phi = -1 @
   prtmp = prtmp[1]~prtmp[2]~"NA"~"NA"~"NA"~"NA";
   call printfm(prtmp,0~1~0~0~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10
               |"s"~10~10|"s"~10~10);"";
else;
   call printfm(prtmp,0~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
               |"lf"~10~4|"lf"~10~4);"";
endif;"";

prtmp = pout[k1+k2+2:rows(pout),.];
"Short-Run Coefficients";"";
call printfm(prtmp,0~1~1~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4
            |"lf"~10~4|"lf"~10~4);"";
if rows(prtmp) eq 1; "";
endif;"";

"                  Summary statistics and diagnostics";
"                  ----------------------------------";"";
pout = DFEDRES[2:6]';
"            RBARSQ     SIGMA      LL        AIC       SC";"";
call printfm(POUT,1~1~1~1~1,"lf"~18~3|"lf"~10~3|"lf"~10~2|"lf"~10~2|"lf"~10~2);
"";"";

pout = DFEDRES[8:11]';
"             CH-SC     CH-FF     CH-NO     CH-HE";"";
call printfm(POUT,1~1~1~1,"lf"~18~2|"lf"~10~2|"lf"~10~2|"lf"~10~2);"";"";

RETP;ENDP;

/* ************************************************************* */
   PROC (0) = PRSUMLR(YNAME,XXname,ms_op,OLSRES,OLSDRES,n,numsi,
              gn_op1,gn_op2,kk,grname,sub_op,n_gr,na_gr,page);
/* ************************************************************* */

Local  iii,c1temp,c2temp,c1tmp,c2tmp,ini,fini,ntmp,c1out,c2out,heading,
       c1aout,c1bout,c1cout,headinga,headingb,ctmp,HA1,HA2;

if     ms_op eq 0;
"      Group-Specific Estimates of the Long-Run Coefficients
       Based on ARDL Specifications With Pre-Specified Lags";
elseif ms_op eq 1;
"      Group-Specific Estimates of the Long-Run Coefficients
       Based on ARDL Secifications Selected Using the Akaike Criterion";
elseif ms_op eq 2;
"      Group-Specific Estimates of the Long-Run Coefficients
       Based on ARDL Specifications Selected Using the Schwarz Criterion";
elseif ms_op eq 3;
"      Group-Specific Estimates of the Long-Run Coefficients
       Based on ARDL Specifications Selected Using the Hannan-Quinn Criterion";
endif;"";

Heading = "Phi"~XXname'~"SIGMA"~"Ch-SC"~"CH-FF"~"CH-NO"~"CH-HE"~"RBARSQ"~"LL";

if kk > 2; HeadingA = "Group"~Heading[1:kk+1];
           HeadingB = Heading[kk+2:cols(Heading)];

"Dependent variable: ";; call printfm(Yname,0,"s"~4~4);"";"";

endif;

IF kk eq 1;

"********   *******************  **************************";
" Group        Estimates           Diagnostic Statistics";
"********   *******************  **************************";
"Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";
call printfm(Heading,0~0~0~0~0~0~0~0~0,
     "s"~16~16|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7
    |"s"~7~7|"s"~7~7);"";"";

ELSEIF kk eq 2;

"********   **************************  **************************";
" Group             Estimates             Diagnostic Statistics";
"********   **************************  **************************";
"Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";

call printfm(Heading,0~0~0~0~0~0~0~0~0~0,
     "s"~16~16|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7
    |"s"~7~7|"s"~7~7|"s"~7~7);"";"";

ELSEIF kk eq 3;
   call printfm(HeadingA,0~0~0~0~0,
        "s"~8~8|"s"~8~8|"s"~7~7|"s"~7~7|"s"~7~7);"";"";

ELSEIF kk eq 4;
   call printfm(HeadingA,0~0~0~0~0~0,
        "s"~8~8|"s"~8~8|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7);"";"";

ELSEIF kk eq 5;
   call printfm(HeadingA,0~0~0~0~0~0~0,
        "s"~8~8|"s"~8~8|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7);"";"";

ELSEIF kk eq 6;
   call printfm(HeadingA,0~0~0~0~0~0~0~0,
        "s"~8~8|"s"~8~8|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7);"";"";

ELSEIF kk eq 7;
   call printfm(HeadingA,0~0~0~0~0~0~0~0~0,
        "s"~8~8|"s"~8~8|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7
       |"s"~7~7);"";"";

ELSEIF kk eq 8;
   call printfm(HeadingA,0~0~0~0~0~0~0~0~0~0,
        "s"~8~8|"s"~8~8|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7|"s"~7~7
       |"s"~7~7|"s"~7~7);"";"";

ELSEIF kk eq 9;
   call printfm(HeadingA,0~0~0~0~0~0~0~0~0~0~0,
        "s"~8~8|"s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6
       |"s"~6~6|"s"~6~6|"s"~6~6);"";"";

ELSEIF kk eq 10;
   call printfm(HeadingA,0~0~0~0~0~0~0~0~0~0~0~0,
        "s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6
       |"s"~6~6|"s"~6~6|"s"~6~6|"s"~6~6);"";"";

ENDIF;

ini = 1;
fini = 0;
iii = 1;
do until iii > n;

   ntmp = kk + 1 + kk + numsi[iii];
   fini = fini + ntmp;

   c1tmp = OLSRES[ini:fini,1];
   c1tmp = c1tmp[kk+1]|c1tmp[1:kk];
   c1tmp = c1tmp';
   c2tmp = OLSDRES[3,iii]~OLSDRES[8:11,iii]'~OLSDRES[2,iii]~OLSDRES[4,iii];
   c1temp = c1tmp~c2tmp;

   @ in the order of phi, theta(1), ..., theta(k), SIGMA,
     chi-SC, chi-FF, chi-NO, chi-HE, RBARSQ, LL @

   c2temp = OLSRES[ini:fini,2];
   c2temp = c2temp[kk+1]|c2temp[1:kk];
   c2temp = c2temp';                   @ standard errors @

   ini = ini + ntmp;

   if gn_op1 eq 1;
      if     sub_op eq 1; c1out = iii~c1temp;
      elseif sub_op eq 2; c1out = n_gr[iii]~c1temp;
      endif;
   elseif gn_op1 eq 2;
      if     sub_op eq 1; c1out = grname[iii]~c1temp;
      elseif sub_op eq 2; c1out = na_gr[iii]~c1temp;
      endif;
   endif;

   if kk > 2;  c1aout = c1out[.,1:2+kk];                      @ estimates @
               c1bout = c1out[.,1]~c1out[.,3+kk:cols(c1out)]; @ diagnostics @
               if iii eq 1;   c1cout = c1bout;
               else;          c1cout = c1cout|c1bout;
               endif;

   endif;

   if kk eq 1;

      if gn_op2 eq 1;
         call printfm(c1out,1~1~1~1~1~1~1~1~1~1,
              "lf"~8~0|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~2|"lf"~7~2
              |"lf"~7~2|"lf"~7~2|"lf"~7~2|"lf"~7~2);
      elseif gn_op2 eq 2;
         call printfm(c1out,0~1~1~1~1~1~1~1~1~1,
              "s"~8~8|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~2|"lf"~7~2
              |"lf"~7~2|"lf"~7~2|"lf"~7~2|"lf"~7~2);
      endif;"";

   elseif kk eq 2;

      if gn_op2 eq 1;
         call printfm(c1out,1~1~1~1~1~1~1~1~1~1~1,
              "lf"~8~0|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~2|"lf"~7~2
             |"lf"~7~2|"lf"~7~2|"lf"~7~2|"lf"~7~2);
      elseif gn_op2 eq 2;
         call printfm(c1out,0~1~1~1~1~1~1~1~1~1~1,
              "s"~8~8|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~2|"lf"~7~2
             |"lf"~7~2|"lf"~7~2|"lf"~7~2|"lf"~7~2);
      endif;"";

   elseif kk eq 3;

      if gn_op2 eq 1;
         call printfm(c1aout,1~1~1~1~1,
              "lf"~8~0|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3);
      elseif gn_op2 eq 2;
         call printfm(c1aout,0~1~1~1~1,
              "s"~8~8|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3);
      endif;"";

   elseif kk eq 4;

      if gn_op2 eq 1;
         call printfm(c1aout,1~1~1~1~1~1,
              "lf"~8~0|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3);
      elseif gn_op2 eq 2;
         call printfm(c1aout,0~1~1~1~1~1,
              "s"~8~8|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3);
      endif;"";

   elseif kk eq 5;

      if gn_op2 eq 1;
         call printfm(c1aout,1~1~1~1~1~1~1,
              "lf"~8~0|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3
             |"lf"~7~3);
      elseif gn_op2 eq 2;
         call printfm(c1aout,0~1~1~1~1~1~1,
              "s"~8~8|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3
             |"lf"~7~3);
      endif;"";

   elseif kk eq 6;

      if gn_op2 eq 1;
         call printfm(c1aout,1~1~1~1~1~1~1~1,
              "lf"~8~0|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3
              |"lf"~7~3);
      elseif gn_op2 eq 2;
         call printfm(c1aout,0~1~1~1~1~1~1~1,
              "s"~8~8|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3
              |"lf"~7~3);
      endif;"";

/*      IF c1temp[1] eq -1;"";   @ if phi = -1 @
      endif;
*/

   elseif kk eq 7;

      if gn_op2 eq 1;
         call printfm(c1aout,1~1~1~1~1~1~1~1~1,
              "lf"~8~0|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3
             |"lf"~7~3|"lf"~7~3);
      elseif gn_op2 eq 2;
         call printfm(c1aout,0~1~1~1~1~1~1~1~1,
              "s"~8~8|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3
             |"lf"~7~3|"lf"~7~3);
      endif;"";

   elseif kk eq 8;

      if gn_op2 eq 1;
         call printfm(c1aout,1~1~1~1~1~1~1~1~1~1,
           "lf"~8~0|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3
          |"lf"~7~3|"lf"~7~3|"lf"~7~3);
      elseif gn_op2 eq 2;
         call printfm(c1aout,0~1~1~1~1~1~1~1~1~1,
           "s"~8~8|"lf"~8~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3|"lf"~7~3
          |"lf"~7~3|"lf"~7~3|"lf"~7~3);
      endif;"";

   elseif kk eq 9;

      if gn_op2 eq 1;
         call printfm(c1aout,1~1~1~1~1~1~1~1~1~1~1,
           "lf"~8~0|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2
          |"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2);
      elseif gn_op2 eq 2;
         call printfm(c1aout,0~1~1~1~1~1~1~1~1~1~1,
           "s"~8~8|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2
          |"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2);
      endif;"";

   elseif kk eq 10;

      if gn_op2 eq 1;
         call printfm(c1aout,1~1~1~1~1~1~1~1~1~1~1~1,
           "lf"~6~0|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2
          |"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2);
      elseif gn_op2 eq 2;
         call printfm(c1aout,0~1~1~1~1~1~1~1~1~1~1~1,
           "s"~6~6|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2
          |"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2|"lf"~6~2);
      endif;"";

   endif;

   IF c1temp[1] eq -1;   @ if phi = -1 @
      if     kk eq 1;
         c2out = "(NA)  ("~c2temp[2]~")";
         call printfm(c2out,0~1~0,
              "s"~18~18|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 2;
         c2out = "(NA)  ("~c2temp[2]~")("~c2temp[3]~")";
         call printfm(c2out,0~1~0~1~0,
              "s"~18~18|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 3;
         c2out = "(NA)  ("~c2temp[2]~")("~c2temp[3]~")("~c2temp[4]~")";
         call printfm(c2out,0~1~0~1~0~1~0,
              "s"~18~18|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 4;
         c2out = "(NA)  ("~c2temp[2]~")("~c2temp[3]~")("~c2temp[4]~")("
                 ~c2temp[5]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0,
              "s"~18~18|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
              |"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 5;
         c2out = "(NA)  ("~c2temp[2]~")("~c2temp[3]~")("~c2temp[4]~")("
                 ~c2temp[5]~")("~c2temp[6]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0,
              "s"~18~18|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
              |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 6;
         c2out = "(NA)  ("~c2temp[2]~")("~c2temp[3]~")("~c2temp[4]~")("
                 ~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~18~18|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
              |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 7;
         c2out = "(NA)  ("~c2temp[2]~")("~c2temp[3]~")("~c2temp[4]~")("
                 ~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]~")("~c2temp[8]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~18~18|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
              |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
              |"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 8;
         c2out = "(NA)  ("~c2temp[2]~")("~c2temp[3]~")("~c2temp[4]~")("
                 ~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]~")("~c2temp[8]~")("
                 ~c2temp[9]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~18~18|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
              |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
              |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 9;
         c2out = "(NA)  ("~c2temp[2]~")("~c2temp[3]~")("~c2temp[4]~")("
                 ~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]~")("~c2temp[8]~")("
                 ~c2temp[9]~")("~c2temp[10]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~16~16|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
              |"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
              |"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~1~1);"";
      elseif kk eq 10;
         c2out = "(NA)  ("~c2temp[2]~")("~c2temp[3]~")("~c2temp[4]~")("
                 ~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]~")("~c2temp[8]~")("
                 ~c2temp[9]~")("~c2temp[10]~")("~c2temp[11]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~14~14|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
              |"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
              |"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
              |"s"~2~2|"lf"~4~2|"s"~1~1);"";
      endif;

   ELSE;

      if     kk eq 1;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")";
         call printfm(c2out,0~1~0~1~0,
              "s"~11~11|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 2;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")("~c2temp[3]~")";
         call printfm(c2out,0~1~0~1~0~1~0,
              "s"~11~11|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 3;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")("~c2temp[3]~")("
                ~c2temp[4]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0,
              "s"~11~11|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
             |"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 4;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")("~c2temp[3]~")("
                ~c2temp[4]~")("~c2temp[5]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0,
              "s"~11~11|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
             |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 5;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")("~c2temp[3]~")("
                ~c2temp[4]~")("~c2temp[5]~")("~c2temp[6]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~11~11|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
             |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 6;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")("~c2temp[3]~")("
                ~c2temp[4]~")("~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~11~11|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
             |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
             |"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 7;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")("~c2temp[3]~")("
                ~c2temp[4]~")("~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]
                ~")("~c2temp[8]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~11~11|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
             |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
             |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 8;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")("~c2temp[3]~")("
                ~c2temp[4]~")("~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]
                ~")("~c2temp[8]~")("~c2temp[9]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~11~11|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
             |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3
             |"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~2~2|"lf"~5~3|"s"~1~1);"";
      elseif kk eq 9;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")("~c2temp[3]~")("
                ~c2temp[4]~")("~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]
                ~")("~c2temp[8]~")("~c2temp[9]~")("~c2temp[10]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~9~9|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
             |"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
             |"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
             |"s"~2~2|"lf"~4~2|"s"~1~1);"";
      elseif kk eq 10;
         c2out = "("~c2temp[1]~")("~c2temp[2]~")("~c2temp[3]~")("
                ~c2temp[4]~")("~c2temp[5]~")("~c2temp[6]~")("~c2temp[7]
                ~")("~c2temp[8]~")("~c2temp[9]~")("~c2temp[10]~")("
                ~c2temp[11]~")";
         call printfm(c2out,0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0~1~0,
              "s"~7~7|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
             |"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
             |"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2
             |"s"~2~2|"lf"~4~2|"s"~2~2|"lf"~4~2|"s"~1~1);"";
      endif;
   ENDIF;

   iii = iii + 1;
endo;"";

"Figures in brackets are the standard errors.";"";

if kk > 2;  ""; page = page + 1; format /rd 3,0;
            "PAGE "  page;"";

"********           *****************************";
" Group                 Diagnostic Statistics";
"********           *****************************";

   call printfm(Headingb,0~0~0~0~0~0~0,
        "s"~16~16|"s"~8~8|"s"~8~8|"s"~8~8|"s"~8~8|"s"~8~8|"s"~8~8);"";"";
   iii = 1;
   do until iii > n;

            ctmp = c1cout[iii,.];
            if gn_op2 eq 1;
               call printfm(ctmp,1~1~1~1~1~1~1~1,
                   "lf"~8~0|"lf"~8~3|"lf"~8~2|"lf"~8~2
                  |"lf"~8~2|"lf"~8~2|"lf"~8~2|"lf"~8~2);"";
            elseif gn_op2 eq 2;
               call printfm(ctmp,0~1~1~1~1~1~1~1,
                   "s"~8~8|"lf"~8~3|"lf"~8~2|"lf"~8~2
                  |"lf"~8~2|"lf"~8~2|"lf"~8~2|"lf"~8~2);"";
            endif;
            iii = iii + 1;
   endo;
endif;"";


RETP;ENDP;

/* ************************************************************************** */
   PROC (0) = PRAVG1(ms_op,OLSRES,PMLERES,n,numsi,kk,plag,qlag,CFname,maxlag,
                     yname,xname,zname,MGECOV,PMLECOV);
/* *************************************************************
   Print the mean group estimates of the LR and SR parameters
   ********************************************************** */

Local ini,fini,oini,ofini,msi,otmp,omatrix,pmatrix,Otemp,Ptemp,
      pini,pfini,bini,bfini,temp1,temp2,temp3,temp4,matrix1,matrix2,matrix3,
      matrix4,OARES,PARES,npara,pp,coname,prout,maxp,maxq,maxpara,sumq,
      otmp0,otmp1,otmp2,otmp3,otmp4,temp5,temp6,temp7,maxot,n_ot,pot,qot,
      junk,thetares,prtmp1,prtmp2,prtmp3,prtmp4,hresult,jresult,rflag,hflag;

IF ms_op eq 0;   @ the case of fixed lags:
                   same number of parameters for each group @

   @ Computing MGE of the OLS parameters @

   ini = 1;   fini = 0;
   msi = 1;
   do until msi > n;

      otmp = kk + 1 + kk + numsi[msi];
      fini = fini + otmp;
      if msi eq 1; Omatrix = OLSRES[ini:fini,1];
      else;        Omatrix = Omatrix~OLSRES[ini:fini,1];  @ all para by N @
      endif;
      ini = ini + otmp;

      msi = msi + 1;
   endo;

   npara = rows(Omatrix);  @ total # of parameters @
   OARES = zeros(npara,3); @ avg, se and t @
   pp = plag[1];           @ same p lag order for each group @
   msi = 1;
   do until msi > npara;

      if pp eq 0;
         if msi eq kk+1;  OARES[msi,.] = Omatrix[msi,1]~"NA"~"NA";
         else;            Otemp = Omatrix[msi,.]';
                          { OARES[msi,.] } = AVG1(Otemp);
         endif;
      else;
         Otemp = Omatrix[msi,.]';
         { OARES[msi,.] } = AVG1(Otemp);
      endif;

      msi = msi + 1;
   endo;

   @ Computing MGE of the Pooled MLE parameters @

   Pini = kk + 1;            Pfini = kk;
   Bini = kk + n + 1;        Bfini = kk + n;
   Oini = kk + n + kk*n + 1; Ofini = kk + n + kk*n;

   msi = 1;
   do until msi > n;

      temp1 = PMLERES[1:kk,1];             @ theta: kx1 @

      Pfini = Pfini + 1;
      temp2 = PMLERES[Pini:Pfini,1];       @ phi(i): 1x1 @
      Pini  = Pini + 1;

      Bfini = Bfini + kk;
      temp3 = PMLERES[Bini:Bfini,1];       @ beta(1i), ..., beta(ki): kx1 @
      Bini  = Bini + kk;

      Ofini = Ofini + numsi[msi];
      temp4 = PMLERES[Oini:Ofini,1];       @ other: numsi[i]x1 @
      Oini  = Oini + numsi[msi];

      if msi eq 1; Matrix1 = temp1;
                   Matrix2 = temp2;
                   Matrix3 = temp3;
                   Matrix4 = temp4;
      else;        Matrix1 = Matrix1~temp1;       @ kxN @
                   Matrix2 = Matrix2~temp2;       @ 1xN @
                   Matrix3 = Matrix3~temp3;       @ kxN @
                   Matrix4 = Matrix4~temp4;       @ otherxN @
      endif;

      msi = msi + 1;
   endo;

   PMatrix = Matrix1|Matrix2|Matrix3|matrix4;

   npara = rows(Pmatrix);  @ total # of parameters @
   PARES = zeros(npara,3); @ avg, se and t @
   msi = 1;
   do until msi > npara;
         Ptemp = Pmatrix[msi,.]';
         { PARES[msi,.] } = AVG1(Ptemp);
         msi = msi + 1;
   endo;

   @ printing the results @

   otmp = kk + 1 + kk + numsi[1];
   coname = cfname[1:otmp];
   prout = coname~PARES~OARES;

   @ the case of the same theta @

   thetares = coname[1:kk]~pmleres[1:kk,.]~prout[1:kk,5:7];
   prout = thetares|prout[kk+1:rows(prout),.];

"          **********************    **********************";
"           Pooled MGE Estimates         MGE Estimates";
"          **********************    **********************";
"             Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";

"           Coef.  St. Er. t-ratio    Coef. St. Er. t-ratio h-test  p-val";"";

   if plag eq 0;  @ if all phi's are zero @

      "Long-Run Coefficients";"";

      prtmp1 = prout[1:kk,.];
      { Hresult, hflag } = HAUSMAN(prtmp1); @ Hausman test @
      if hflag eq 0;
         prtmp1 = prtmp1~hresult;
        call printfm(prtmp1,0~1~1~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3|"lf"~7~2|"lf"~7~2);"";
      elseif hflag eq 1;
         prtmp1 = prtmp1;
        call printfm(prtmp1,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      endif;
      if     rows(prtmp1) eq 1; "";
      elseif rows(prtmp1) > 1;
         { Jresult, Rflag }
         = JHAUSMAN(prtmp1,MGECOV,PMLECOV); @ Joint Hausman test @
         if rflag eq 0;
"                                      Joint Hausman test: ";;
            call printfm(Jresult,1~1,"lf"~7~2|"lf"~7~2);"";
         elseif rflag eq 1;
"        ** V(MGE) - V(PMLE) is not positive definite. **";
         endif;

      endif;

      "Error Correction Coefficients";"";
      prtmp2 = prout[kk+1,1:2]~"NA"~"NA"~prout[kk+1,5]~"NA"~"NA";
      call printfm(prtmp2,0~1~0~0~1~0~0,"s"~8~8|"lf"~8~3|"s"~8~8|"s"~8~8
                   |"lf"~8~3|"s"~8~8|"s"~8~8);"";"";

      "Short-Run Coefficients";"";
      prtmp3 = prout[kk+2:kk+1+kk,1]~pmleres[1:kk,.]~prout[kk+2:kk+1+kk,5:7];
      call printfm(prtmp3,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);
      if rows(prtmp3) eq 1; "";
      endif;

      prtmp4 = prout[kk+1+kk+1:rows(prout),.];
      call printfm(prtmp4,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      if rows(prtmp1) eq 1; "";
      endif;

   else;
      PRAVG1S(prout,kk,MGECOV,PMLECOV);
   endif;

ELSE;   @ different number of parameters for each group @

   @ Computing MGE of the OLS parameters @

   if maxlag eq 0; maxot = rows(zname);
   else;           maxot = maxlag*(kk+1) - 1 + rows(zname);
   endif;

   @ maximum # of other parameters @

   maxpara = kk + 1 + kk + maxot;  @ maximum number of parameters @
   Omatrix = zeros(maxpara,n);
   ini = 1;   fini = 0;
   msi = 1;
   do until msi > n;

      otmp = kk + 1 + kk + numsi[msi];
      fini = fini + otmp;
      otmp0 = OLSRES[ini:fini,1];
      ini = ini + otmp;

      if numsi[msi] eq maxot;  @ if there are all parameters @
         Omatrix[.,msi] = otmp0;

      else;  @ if less than all parameters @
           otmp1 = otmp0[1:kk+1+kk];    @ theta, phi, beta @

           if numsi[msi] eq 1;   @ if only intercept @
              /* "CHK" maxot; */
              otmp3 = zeros(maxot-1,1);
           else;
              otmp2 = otmp0[kk+1+kk+1:rows(otmp0)-1];  @ other parameters @
              pot = plag[msi]; qot = qlag[msi,.];
             { otmp3 } = OTHERP(otmp2,pot,qot,maxlag,maxot,kk);
           endif;

           @ fixed regressors @
           otmp4 = otmp0[rows(otmp0)-rows(zname)+1:rows(otmp0)];

           Omatrix[.,msi] = otmp1|otmp3|otmp4;
      endif;

      msi = msi + 1;
   endo;

   OARES = zeros(maxpara,3); @ avg, se and t @
   msi = 1;
   do until msi > maxpara;
      otemp = omatrix[msi,.]';
      { OARES[msi,.] } = AVG1(Otemp);
      msi = msi + 1;
   endo;

   @ Computing MGE of the Pooled MLE parameters @

   Pmatrix = zeros(maxpara,n);
   Pini = kk + 1;            Pfini = kk;
   Bini = kk + n + 1;        Bfini = kk + n;
   Oini = kk + n + kk*n + 1; Ofini = kk + n + kk*n;

   msi = 1;
   do until msi > n;

      temp1 = PMLERES[1:kk,1];             @ theta: kx1 @

      Pfini = Pfini + 1;
      temp2 = PMLERES[Pini:Pfini,1];       @ phi(i): 1x1 @
      Pini  = Pini + 1;

      Bfini = Bfini + kk;
      temp3 = PMLERES[Bini:Bfini,1];       @ beta(1i), ..., beta(ki): kx1 @
      Bini  = Bini + kk;

      Ofini = Ofini + numsi[msi];
      temp4 = PMLERES[Oini:Ofini,1];       @ other: numsi[i]x1 @
      Oini  = Oini + numsi[msi];

      if numsi[msi] eq maxot;  @ if there are all parameters @
         Pmatrix[.,msi] = temp1|temp2|temp3|temp4;

      else;  @ if less than all parameters @

         if numsi[msi] eq 1;   @ if only intercept @
            temp6 = zeros(maxot-1,1);
         else;
            temp5 = temp4[1:rows(temp4)-1];
            pot = plag[msi]; qot = qlag[msi,.];
            { temp6 } = OTHERP(temp5,pot,qot,maxlag,maxot,kk);
         endif;

         @ fixed regressors @
         temp7 = temp4[rows(temp4)-rows(zname)+1:rows(temp4)];
         Pmatrix[.,msi] = temp1|temp2|temp3|temp6|temp7;

      endif;

      msi = msi + 1;
   endo;

   PARES = zeros(maxpara,3); @ avg, se and t @
   msi = 1;
   do until msi > maxpara;
         Ptemp = Pmatrix[msi,.]';
         { PARES[msi,.] } = AVG1(Ptemp);
         msi = msi + 1;
   endo;

   @ printing the results @

    {COname} = MXCNAME1(maxlag,kk,yname,xname,zname);
   prout = coname~PARES~OARES;

   @ the case of the same theta @

   thetares = coname[1:kk]~pmleres[1:kk,.]~prout[1:kk,5:7];
   prout = thetares|prout[kk+1:rows(prout),.];

"         ************************    **********************";
"           Pooled MGE Estimates         MGE Estimates";
"         ************************    **********************";
"             Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";"";

"           Coef.  St. Er. t-ratio    Coef. St. Er. t-ratio h-test p-val";"";

   if plag eq 0;  @ if all phi's are zero @

      "Long-Run Coefficients";"";
      prtmp1 = prout[1:kk,.];
      { Hresult, hflag } = HAUSMAN(prtmp1); @ Hausman test @
      if hflag eq 0;
         prtmp1 = prtmp1~hresult;
        call printfm(prtmp1,0~1~1~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3|"lf"~7~2|"lf"~7~2);"";
      elseif hflag eq 1;
         prtmp1 = prtmp1;
        call printfm(prtmp1,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      endif;
      if rows(prtmp1) eq 1; "";
      elseif rows(prtmp1) > 1;
         { Jresult, Rflag }
         = JHAUSMAN(prtmp1,MGECOV,PMLECOV); @ Joint Hausman test @
         if rflag eq 0;
"                                      Joint Hausman test: ";;
            call printfm(Jresult,1~1,"lf"~7~2|"lf"~7~2);"";
         elseif rflag eq 1;
"        ** V(MGE) - V(PMLE) is not positive definite. **";
         endif;
      endif;

      "Error Correction Coefficients";"";
      prtmp2 = prout[kk+1,1:2]~"NA"~"NA"~prout[kk+1,5]~"NA"~"NA";
      call printfm(prtmp2,0~1~0~0~1~0~0,"s"~8~8|"lf"~8~3|"s"~8~8|"s"~8~8
                   |"lf"~10~3|"s"~8~8|"s"~8~8);"";"";

      "Short-Run Coefficients";"";
      prtmp3 = prout[kk+2:kk+1+kk,1]~pmleres[1:kk,.]~prout[kk+2:kk+1+kk,5:7];
      call printfm(prtmp3,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";

      prtmp4 = prout[kk+1+kk+1:rows(prout),.];
      call printfm(prtmp4,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      if rows(prtmp4) eq 1; "";
      endif;

   else;
      PRAVG1S(prout,kk,MGECOV,PMLECOV);
   endif;

ENDIF;

"";
"* h-test is the Hausman test statistic with its associated p-value.";"";

RETP;ENDP;

/* ******************************************* */
   PROC (0) = PRAVG1S(prout,k,MGECOV,PMLECOV);
/* ******************************************* */

Local prtmp,hresult,jresult,rflag,hflag;

"Long-Run Coefficients";"";
prtmp = prout[1:k,.];
{ Hresult, hflag } = HAUSMAN(prtmp); @ Hausman test @
if hflag eq 0;
   prtmp = prtmp~hresult;
   call printfm(prtmp,0~1~1~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3|"lf"~7~2|"lf"~7~2);"";
elseif hflag eq 1;
   prtmp = prtmp;
   call printfm(prtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
 endif;

if rows(prtmp) eq 1; "";
elseif rows(prtmp) > 1;
      { Jresult, Rflag }
      = JHAUSMAN(prtmp,MGECOV,PMLECOV); @ Joint Hausman test @
      if rflag eq 0;
"                                      Joint Hausman test: ";;
         call printfm(Jresult,1~1,"lf"~7~2|"lf"~7~2);"";
      elseif rflag eq 1;
"        ** V(MGE) - V(PMLE) is not positive definite. **";
      endif;
endif;"";

prtmp = prout[k+1,.];
"Error Correction Coefficients";"";
if prtmp[2] eq -1;  @ if phi = -1 @
   prtmp = prtmp[1]~prtmp[2]~"NA"~"NA"
           ~prtmp[5]~prtmp[6]~prtmp[7];
   call printfm(PRtmp,0~1~0~0~1~0~0,"s"~8~8|"lf"~8~3|"s"~8~8|"s"~8~8
               |"lf"~10~3|"s"~8~8|"s"~8~8);"";
else;
   call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
               |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
endif;"";

prtmp = prout[k+2:rows(prout),.];
"Short-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
            |"lf"~10~3|"lf"~8~3|"lf"~8~3);
if rows(prtmp) eq 1; "";
endif;"";

RETP;ENDP;

/* ************************* */
   Proc (2) = HAUSMAN(htmp);
/* ************************* */

   Local nh,htest,pval,hi,htmp0,htmp1,htmp2,hresult,hflag;

   nh = rows(htmp);
   htest = zeros(nh,1);   @ Hausman test @
   pval = zeros(nh,1);    @ p-value @
   hflag = 0;
   hi = 1;
   do until hi > nh;

         htmp0 = htmp[hi,.];
         htmp1 = (htmp0[5] - htmp0[2])^2;  @ theta(c) - theta(e) @
         htmp2 = htmp0[6]^2 - htmp0[3]^2;  @ V(t_c) - V(t_e) @
         if htmp2 < 0.00000001; hflag = 1; goto Hout;
         endif;
         htest[hi] = htmp1/htmp2;
         pval[hi] = cdfchic(htest[hi],1);

         hi = hi + 1;
   endo;

HOUT:

   hresult = htest~pval;

   RETP(hresult,hflag); ENDP;

/* ****************************************************************** */
   Proc (2) = JHAUSMAN(jtmp,MGECOV,PMLECOV); @ Joint Hausman test @
/* ****************************************************************** */
   Local nj,jtest,pval,ji,jtmp1,jresult,chk,RK,Vdiff,rflag;

   nj = rows(jtmp);
   jtmp1 = zeros(nj,1);   @ kx1 vector of theta(c) - theta(e) @
   ji = 1;
   do until ji > nj;

      jtmp1[ji] = (jtmp[ji,5] - jtmp[ji,2]);  @ theta(c) - theta(e) @
      ji = ji + 1;
   endo;

   Vdiff = mgecov - pmlecov;
   RK = det(Vdiff);
   if rk > 0.00000001; /* "CHK: Y"; */
      Rflag = 0;
      jtest = jtmp1'invpd(Vdiff)*jtmp1;
      pval = cdfchic(jtest,nj);
   else; Rflag =1; jtest = 0; pval = 0;  /* "CHK: N"; */
   endif;
   jresult = jtest~pval;

   RETP(Jresult,Rflag); ENDP;

/* ************************************************************************** */
   PROC (0) =
   PRAVG2(ms_op,OLSRES,PMLERES,n,numsi,k1,k2,plag,q1lag,q2lag,CFname,maxlag,
          yname,x1name,x2name,zname,MGECOV,PMLECOV);
/* ***************************************************************************
   Print the mean group estimates of the LR and SR parameters
   for the ML option case 2
   ********************************************************** */

Local ini,fini,oini,ofini,msi,otmp,omatrix,pmatrix,Otemp,Ptemp,
      pini,pfini,b1ini,b1fini,temp1,temp2,temp3,temp4,matrix1,matrix2,matrix3,
      t2ini,t2fini,b2ini,b2fini,matrix4,matrix5,matrix6,temp5,temp6,temp7,
      temp8,temp9,OARES,PARES,npara,pp,coname,prout,maxp,maxq,maxpara,sumq,
      otmp0,otmp1,otmp2,otmp3,otmp4,maxot,n_ot,pot,qot,
      junk,thetares,kk,prtmp1,prtmp2,prtmp3,prtmp4,hresult,jresult,rflag,hflag;

kk = k1 + k2;
IF ms_op eq 0;   @ same number of parameters for each group @

   @ Computing MGE of the OLS parameters @

   ini = 1;   fini = 0;
   msi = 1;
   do until msi > n;

      otmp = kk + 1 + kk + numsi[msi];
      fini = fini + otmp;
      if msi eq 1; Omatrix = OLSRES[ini:fini,1];
      else;        Omatrix = Omatrix~OLSRES[ini:fini,1];  @ all para by N @
      endif;
      ini = ini + otmp;

      msi = msi + 1;
   endo;

   npara = rows(Omatrix);  @ total # of parameters @
   OARES = zeros(npara,3); @ avg, se and t @
   pp = plag[1];           @ same p lag order for each group @
   msi = 1;
   do until msi > npara;

      if pp eq 0;
         if msi eq kk+1;  OARES[msi,.] = Omatrix[msi,1]~"NA"~"NA";
         else;            Otemp = Omatrix[msi,.]';
                          { OARES[msi,.] } = AVG1(Otemp);
         endif;
      else;
         Otemp = Omatrix[msi,.]';
         { OARES[msi,.] } = AVG1(Otemp);
      endif;

      msi = msi + 1;
   endo;

   @ Computing MGE of the Pooled MLE parameters @

   T2ini = k1 + 1;                   T2fini = k1;
   Pini  = k1 + k2*n + 1;            Pfini  = k1 + k2*n;
   B1ini = k1 + k2*n + n + 1;        B1fini = k1 + k2*n + n;
   B2ini = k1 + k2*n + n + k1*n + 1; B2fini = k1 + k2*n + n + k1*n;
   Oini  = k1 + k2*n + n + kk*n + 1; Ofini  = k1 + k2*n + n + kk*n;

   msi = 1;
   do until msi > n;

      temp1 = PMLERES[1:k1,1];         @ theta1: k1x1 @

      T2fini = T2fini + k2;
      temp2 = PMLERES[T2ini:T2fini,1]; @ theta2(i): k2x1 @
      T2ini  = T2ini + k2;

      Pfini = Pfini + 1;
      temp3 = PMLERES[Pini:Pfini,1];   @ phi(i): 1x1 @
      Pini  = Pini + 1;

      B1fini = B1fini + k1;
      temp4  = PMLERES[B1ini:B1fini,1];   @ beta1(1i), ..., beta1(k1i): k1x1 @
      B1ini  = B1ini + k1;

      B2fini = B2fini + k2;
      temp5  = PMLERES[B2ini:B2fini,1];   @ beta2(1i), ..., beta2(k2i): k2x1 @
      B2ini  = B2ini + k2;

      Ofini = Ofini + numsi[msi];
      temp6 = PMLERES[Oini:Ofini,1];       @ other: numsi[i]x1 @
      Oini  = Oini + numsi[msi];

      if msi eq 1; Matrix1 = temp1;
                   Matrix2 = temp2;
                   Matrix3 = temp3;
                   Matrix4 = temp4;
                   Matrix5 = temp5;
                   Matrix6 = temp6;
      else;        Matrix1 = Matrix1~temp1;       @ k1xN @
                   Matrix2 = Matrix2~temp2;       @ k2xN @
                   Matrix3 = Matrix3~temp3;       @ 1xN @
                   Matrix4 = Matrix4~temp4;       @ k1xN @
                   Matrix5 = Matrix5~temp5;       @ k2xN @
                   Matrix6 = Matrix6~temp6;       @ otherxN @
      endif;

      msi = msi + 1;
   endo;

   PMatrix = Matrix1|Matrix2|Matrix3|matrix4|matrix5|matrix6;

   npara = rows(Pmatrix);  @ total # of parameters @
   PARES = zeros(npara,3); @ avg, se and t @
   msi = 1;
   do until msi > npara;
         Ptemp = Pmatrix[msi,.]';
         { PARES[msi,.] } = AVG1(Ptemp);
         msi = msi + 1;
   endo;

   @ printing the results @

   otmp = kk + 1 + kk + numsi[1];
   coname = cfname[1:otmp];
   prout = coname~PARES~OARES;

   @ the case of the same theta @

   thetares = coname[1:k1]~pmleres[1:k1,.]~prout[1:k1,5:7];
   prout = thetares|prout[k1+1:rows(prout),.];

"         ************************     **********************";
"           Pooled MGE Estimates          MGE Estimates";
"         ************************     **********************";
"             Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";
"           Coef.  St. Er. t-ratio    Coef. St. Er. t-ratio h-test  p-val";"";

   if plag eq 0;  @ if all phi's are zero @

      "Long-Run Coefficients Restricted to be the Same Across All Groups";"";
      prtmp1 = prout[1:k1,.];
      { Hresult, hflag } = HAUSMAN(prtmp1); @ Hausman test @
      if hflag eq 0;
         prtmp1 = prtmp1~hresult;
        call printfm(prtmp1,0~1~1~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3|"lf"~7~2|"lf"~7~2);"";
      elseif hflag eq 1;
         prtmp1 = prtmp1;
        call printfm(prtmp1,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      endif;

      if rows(prtmp1) eq 1; "";
      elseif rows(prtmp1) > 1;
         { Jresult, rflag }
         = JHAUSMAN(prtmp1,MGECOV,PMLECOV); @ Joint Hausman test @
         if rflag eq 0;
"                                      Joint Hausman test: ";;
            call printfm(Jresult,1~1,"lf"~7~2|"lf"~7~2);"";
         elseif rflag eq 1;
"        ** V(MGE) - V(PMLE) is not positive definite. **";
         endif;
      endif;

     "Unrestricted Long-Run Coefficients";"";
      prtmp1 = prout[k1+1:k1+k2,.];
      call printfm(prtmp1,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      if rows(prtmp1) eq 1; "";
      endif;

      "Error Correction Coefficients";"";
      prtmp2 = prout[kk+1,1:2]~"NA"~"NA"~prout[kk+1,5]~"NA"~"NA";
      call printfm(prtmp2,0~1~0~0~1~0~0,"s"~8~8|"lf"~8~3|"s"~8~8|"s"~8~8
                   |"lf"~10~3|"s"~8~8|"s"~8~8);"";"";

      "Short-Run Coefficients";"";
      prtmp3 = prout[kk+2:kk+1+k1,1]~pmleres[1:k1,.]~prout[kk+2:kk+1+k1,5:7];
      call printfm(prtmp3,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";

      prtmp4 = prout[kk+1+k1+1:rows(prout),.];
      call printfm(prtmp4,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      if rows(prtmp1) eq 1; "";
      endif;

   else;
      PRAVG2S(prout,k1,k2,MGECOV,PMLECOV);
   endif;

ELSE;   @ different number of parameters for each group @

   @ Computing MGE of the OLS parameters @

   if maxlag eq 0; maxot = rows(zname);
   else;           maxot = maxlag*(kk+1) - 1 + rows(zname);
   endif;
                   @ maximum # of other parameters @
   maxpara = kk + 1 + kk + maxot;  @ maximum number of parameters @
   Omatrix = zeros(maxpara,n);
   ini = 1;   fini = 0;
   msi = 1;
   do until msi > n;

      otmp = kk + 1 + kk + numsi[msi];
      fini = fini + otmp;
      otmp0 = OLSRES[ini:fini,1];
      ini = ini + otmp;

      if numsi[msi] eq maxot;  @ if there are all parameters @
         Omatrix[.,msi] = otmp0;

      else;  @ if less than all parameters @
           otmp1 = otmp0[1:kk+1+kk];    @ theta, phi, beta @

           if numsi[msi] eq 1;   @ if only intercept @
              otmp3 = zeros(maxot-1,1);
           else;
              otmp2 = otmp0[kk+1+kk+1:rows(otmp0)-1];  @ other parameters @
              pot = plag[msi]; qot = q1lag[msi,.]~q2lag[msi,.];
             { otmp3 } = OTHERP(otmp2,pot,qot,maxlag,maxot,kk);
           endif;

           otmp4 = otmp0[rows(otmp0)-rows(zname)+1:rows(otmp0)];
           @ fixed regressors @

           Omatrix[.,msi] = otmp1|otmp3|otmp4;
      endif;

      msi = msi + 1;
   endo;

   OARES = zeros(maxpara,3); @ avg, se and t @
   msi = 1;
   do until msi > maxpara;
      otemp = omatrix[msi,.]';
      { OARES[msi,.] } = AVG1(Otemp);
      msi = msi + 1;
   endo;

   @ Computing MGE of the Pooled MLE parameters @

   Pmatrix = zeros(maxpara,n);

   T2ini = k1 + 1;                   T2fini = k1;
   Pini  = k1 + k2*n + 1;            Pfini  = k1 + k2*n;
   B1ini = k1 + k2*n + n + 1;        B1fini = k1 + k2*n + n;
   B2ini = k1 + k2*n + n + k1*n + 1; B2fini = k1 + k2*n + n + k1*n;
   Oini  = k1 + k2*n + n + kk*n + 1; Ofini  = k1 + k2*n + n + kk*n;

   msi = 1;
   do until msi > n;

      temp1 = PMLERES[1:k1,1];         @ theta1: k1x1 @

      T2fini = T2fini + k2;
      temp2 = PMLERES[T2ini:T2fini,1]; @ theta2(i): k2x1 @
      T2ini  = T2ini + k2;

      Pfini = Pfini + 1;
      temp3 = PMLERES[Pini:Pfini,1];   @ phi(i): 1x1 @
      Pini  = Pini + 1;

      B1fini = B1fini + k1;
      temp4  = PMLERES[B1ini:B1fini,1];   @ beta1(1i), ..., beta1(k1i): k1x1 @
      B1ini  = B1ini + k1;

      B2fini = B2fini + k2;
      temp5  = PMLERES[B2ini:B2fini,1];   @ beta2(1i), ..., beta2(k2i): k2x1 @
      B2ini  = B2ini + k2;

      Ofini = Ofini + numsi[msi];
      temp6 = PMLERES[Oini:Ofini,1];       @ other: numsi[i]x1 @
      Oini  = Oini + numsi[msi];

      if numsi[msi] eq maxot;  @ if there are all parameters @
         Pmatrix[.,msi] = temp1|temp2|temp3|temp4|temp5|temp6;

      else;  @ if less than all parameters @

         if numsi[msi] eq 1;   @ if only intercept @
            temp8 = zeros(maxot-1,1);
         else;
            temp7 = temp6[1:rows(temp6)-1];
            pot = plag[msi]; qot = q1lag[msi,.]~q2lag[msi,.];
            { temp8 } = OTHERP(temp7,pot,qot,maxlag,maxot,kk);
         endif;

         temp9 = temp6[rows(temp6)-rows(zname)+1:rows(temp6)];
         @ fixed regressors @
         Pmatrix[.,msi] = temp1|temp2|temp3|temp4|temp5|temp8|temp9;

      endif;

      msi = msi + 1;
   endo;

   PARES = zeros(maxpara,3); @ avg, se and t @
   msi = 1;
   do until msi > maxpara;
         Ptemp = Pmatrix[msi,.]';
         { PARES[msi,.] } = AVG1(Ptemp);
         msi = msi + 1;
   endo;

   @ printing the results @

    {COname} = MXCNAME2(maxlag,k1,k2,yname,x1name,x2name,zname);
    prout = coname~PARES~OARES;

   @ the case of the same theta @

    thetares = coname[1:k1]~pmleres[1:k1,.]~prout[1:k1,5:7];
    prout = thetares|prout[k1+1:rows(prout),.];

"         ************************    **********************";
"           Pooled MGE Estimates         MGE Estimates";
"         ************************    **********************";
"             Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";
"           Coef.  St. Er. t-ratio    Coef. St. Er. t-ratio h-test  p-val";"";

   if plag eq 0;  @ if all phi's are zero @

      "Long-Run Coefficients Restricted to be the Same Across All Groups";"";
      prtmp1 = prout[1:k1,.];
      { Hresult, hflag } = HAUSMAN(prtmp1); @ Hausman test @
      if hflag eq 0;
         prtmp1 = prtmp1~hresult;
        call printfm(prtmp1,0~1~1~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3|"lf"~7~2|"lf"~7~2);"";
      elseif hflag eq 1;
         prtmp1 = prtmp1;
        call printfm(prtmp1,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      endif;
      if rows(prtmp1) eq 1; "";
      elseif rows(prtmp1) > 1;
         { Jresult, Rflag }
         = JHAUSMAN(prtmp1,MGECOV,PMLECOV); @ Joint Hausman test @
         if rflag eq 0;
"                                      Joint Hausman test: ";;
         call printfm(Jresult,1~1,"lf"~7~2|"lf"~7~2);"";
         elseif rflag eq 1;
"        ** V(MGE) - V(PMLE) is not positive definite. **";
         endif;
      endif;"";

      "Unrestricted Long-Run Coefficients";"";
      prtmp1 = prout[k1+1:k1+k2,.];
      call printfm(prtmp1,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      if rows(prtmp1) eq 1; "";
      endif;"";

      "Error Correction Coefficients";"";
      prtmp2 = prout[kk+1,1:2]~"NA"~"NA"~prout[kk+1,5]~"NA"~"NA";
      call printfm(prtmp2,0~1~0~0~1~0~0,"s"~8~8|"lf"~8~3|"s"~8~8|"s"~8~8
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";"";

      "Short-Run Coefficients";"";
      prtmp3 = prout[kk+2:kk+1+k1,1]~pmleres[1:k1,.]~prout[kk+2:kk+1+k1,5:7];
      call printfm(prtmp3,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";

      prtmp4 = prout[kk+1+k1+1:rows(prout),.];
      call printfm(prtmp4,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
      if rows(prtmp4) eq 1; "";
      endif;"";

   else;
      PRAVG2S(prout,k1,k2,MGECOV,PMLECOV);
   endif;

ENDIF;

"";
"* h-test is the Hausman test statistic with its associated p-value.";"";

RETP;ENDP;

/* *********************************************** */
   PROC (0) = PRAVG2S(prout,k1,k2,MGECOV,PMLECOV);
/* *********************************************** */

Local prtmp,hresult,jresult,rflag,hflag;

"Long-Run Coefficients Restricted to be the Same Across all Groups";"";
prtmp = prout[1:k1,.];
{ Hresult, hflag } = HAUSMAN(prtmp); @ Hausman test @
if hflag eq 0;
   prtmp = prtmp~hresult;
   call printfm(prtmp,0~1~1~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3|"lf"~7~2|"lf"~7~2);"";
elseif hflag eq 1;
   prtmp = prtmp;
   call printfm(prtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
                   |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
endif;
if rows(prtmp) eq 1; "";
elseif rows(prtmp) > 1;
      { Jresult, Rflag }
      = JHAUSMAN(prtmp,MGECOV,PMLECOV); @ Joint Hausman test @
      if rflag eq 0;
"                                      Joint Hausman test: ";;
      call printfm(Jresult,1~1,"lf"~7~2|"lf"~7~2);"";
      elseif rflag eq 1;
"        ** V(MGE) - V(PMLE) is not positive definite. **";
      endif;
endif;"";

prtmp = prout[k1+1:k1+k2,.];
"Unrestricted Long-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
            |"lf"~10~3|"lf"~8~3|"lf"~8~3);
if rows(prtmp) eq 1; "";
endif;"";

prtmp = prout[k1+k2+1,.];
"Error Correction Coefficients";"";
if prtmp[2] eq -1;  @ if phi = -1 @
   prtmp = prtmp[1]~prtmp[2]~"NA"~"NA"
           ~prtmp[5]~prtmp[6]~prtmp[7];
   call printfm(PRtmp,0~1~0~0~1~0~0,"s"~8~8|"lf"~8~3|"s"~8~8|"s"~8~8
               |"lf"~10~3|"s"~8~8|"s"~8~8);"";
else;
   call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
               |"lf"~10~3|"lf"~8~3|"lf"~8~3);"";
endif;"";

prtmp = prout[k1+k2+2:rows(prout),.];
"Short-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1~1~1~1,"s"~8~8|"lf"~8~3|"lf"~8~3|"lf"~8~3
            |"lf"~10~3|"lf"~8~3|"lf"~8~3);
if rows(prtmp) eq 1; "";
endif;"";

RETP;ENDP;

/*  *************************************************  */
    Proc (1) = MXCNAME1(maxlag,kk,yname,xname,zname);
/*  *************************************************************
*   Give all the coefficents name for the case of model selection
*   ************************************************************* */

Local ttmp1,ttmp2,btmp1,btmp2,srname,mxi,s,dyname,dxname,dzname,dxtmp0,dxtmp1,
      dztmp0,dztmp1,mxj,coname,srname1;

   SRname = "ini";
   if     maxlag eq 2;
          s = seqa(1,1,1);
          dyname = 0 $+ "d" $+ yname $+ "(-" $+ ftocv(s,1,0) $+ ")";
          SRname = SRname|dyname;
   elseif maxlag eq 3;
          s = seqa(1,1,2);
          dyname = 0 $+ "d" $+ yname $+ "(-" $+ ftocv(s,1,0) $+ ")";
          SRname = SRname|dyname;
   elseif maxlag eq 4;
          s = seqa(1,1,3);
          dyname = 0 $+ "d" $+ yname $+ "(-" $+ ftocv(s,1,0) $+ ")";
          SRname = SRname|dyname;
   endif;

   mxj = 1;
   do until mxj > kk;

      if     maxlag eq 1;
             dxname = 0 $+ "d" $+ xname[mxj];
             SRname = SRname|dxname;

      elseif maxlag eq 2;
             dxtmp0 = 0 $+ "d" $+ xname[mxj];
             s = seqa(1,1,1);
             dxtmp1 = 0 $+ "d" $+ xname[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      elseif maxlag eq 3;
             dxtmp0 = 0 $+ "d" $+ xname[mxj];
             s = seqa(1,1,2);
             dxtmp1 = 0 $+ "d" $+ xname[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      elseif maxlag eq 4;
             dxtmp0 = 0 $+ "d" $+ xname[mxj];
             s = seqa(1,1,3);
             dxtmp1 = 0 $+ "d" $+ xname[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      endif;

      mxj = mxj + 1;
   endo;

   SRname = SRname|zname;
   SRname = SRname[2:rows(SRname)];

coname = xname|"Phi"|xname|SRname;

RETP(COname);ENDP;

/*  ************************************************************  */
    Proc (1) = MXCNAME2(maxlag,k1,k2,yname,x1name,x2name,zname);
/*  *************************************************************
*   Give all the coefficents name for the case of model selection
*   ************************************************************* */

Local ttmp1,ttmp2,btmp1,btmp2,srname,mxi,s,dyname,dxname,dzname,dxtmp0,dxtmp1,
      dztmp0,dztmp1,mxj,coname,srname1,srname2;

   SRname = "ini";
   if     maxlag eq 2;
          s = seqa(1,1,1);
          dyname = 0 $+ "d" $+ yname $+ "(-" $+ ftocv(s,1,0) $+ ")";
          SRname = SRname|dyname;
   elseif maxlag eq 3;
          s = seqa(1,1,2);
          dyname = 0 $+ "d" $+ yname $+ "(-" $+ ftocv(s,1,0) $+ ")";
          SRname = SRname|dyname;
   elseif maxlag eq 4;
          s = seqa(1,1,3);
          dyname = 0 $+ "d" $+ yname $+ "(-" $+ ftocv(s,1,0) $+ ")";
          SRname = SRname|dyname;
   endif;

   mxj = 1;
   do until mxj > k1;

      if     maxlag eq 1;
             dxname = 0 $+ "d" $+ x1name[mxj];
             SRname = SRname|dxname;

      elseif maxlag eq 2;
             dxtmp0 = 0 $+ "d" $+ x1name[mxj];
             s = seqa(1,1,1);
             dxtmp1 = 0 $+ "d" $+ x1name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      elseif maxlag eq 3;
             dxtmp0 = 0 $+ "d" $+ x1name[mxj];
             s = seqa(1,1,2);
             dxtmp1 = 0 $+ "d" $+ x1name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      elseif maxlag eq 4;
             dxtmp0 = 0 $+ "d" $+ x1name[mxj];
             s = seqa(1,1,3);
             dxtmp1 = 0 $+ "d" $+ x1name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      endif;

      mxj = mxj + 1;
   endo;

   mxj = 1;
   do until mxj > k2;

      if     maxlag eq 1;
             dxname = 0 $+ "d" $+ x2name[mxj];
             SRname = SRname|dxname;

      elseif maxlag eq 2;
             dxtmp0 = 0 $+ "d" $+ x2name[mxj];
             s = seqa(1,1,1);
             dxtmp1 = 0 $+ "d" $+ x2name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      elseif maxlag eq 3;
             dxtmp0 = 0 $+ "d" $+ x2name[mxj];
             s = seqa(1,1,2);
             dxtmp1 = 0 $+ "d" $+ x2name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      elseif maxlag eq 4;
             dxtmp0 = 0 $+ "d" $+ x2name[mxj];
             s = seqa(1,1,3);
             dxtmp1 = 0 $+ "d" $+ x2name[mxj] $+ "(-" $+ ftocv(s,1,0) $+ ")";
             dxname = dxtmp0|dxtmp1;
             SRname = SRname|dxname;

      endif;

      mxj = mxj + 1;
   endo;

   SRname = SRname|zname;
   SRname = SRname[2:rows(SRname)];

coname = x1name|x2name|"Phi"|x1name|x2name|SRname;

RETP(COname);ENDP;

/*  ***********************  */
    PROC (1) = AVG1(tmpvec);
/*  ****************************************************************
*   Compute the average, its standard error in a non-parameteric way
*   **************************************************************** */

Local dimvec,mvec,devvec,sumvec,sevec,tvec,mresult;

dimvec = rows(tmpvec);                        @ dimension of the vector @
mvec = meanc(tmpvec);                         @ mean of vector @
devvec = tmpvec - mvec.*ones(dimvec,1);       @ demeaned vector @
sumvec = sumc(devvec.^2)/(dimvec*(dimvec-1)); @ NP variance of vector @
sevec = sqrt(sumvec);                         @ NP SE of vector @

if sevec < 0.0000000001; tvec = "NA";
else;                   tvec = mvec/sevec;             @ t-ratio @
endif;

if sevec < 0.000000001; mresult = mvec~"NA"~"NA";
else;          mresult = mvec~sevec~tvec;
endif;

RETP(mresult);ENDP;

/*  *************************************************  */
    PROC (1) = OTHERP(ottmp,pot,qot,maxlag,maxot,kk);
/*  ****************************************************************
*   Construct the matrix of all possible other parameters
*   including the nil values
*
*   Input: numsi(i)-1 by 1
*   Output is n_ot-1 by 1.
*   **************************************************************** */

Local otresp,otresq,otres,ti,oi,oj,qi;

ti = 1;   @ index for inpput: numsi(i)-1 by 1 @

if maxlag < 2;  goto nextq;
else; otresp = zeros(maxlag-1,1);
endif;

if pot < 2;
else;
      oj = 1;
      do until oj > pot-1;
         otresp[oj] = ottmp[ti];
         ti = ti + 1;
         oj = oj + 1;
      endo;
endif;

NEXTQ:

otresq = zeros(maxlag,kk);
qi = 1;
do until qi > kk;

   if qot[qi] eq 0;
   else; oj = 1;
         do until oj > qot[qi];
            otresq[oj,qi] = ottmp[ti];
            ti = ti + 1;
            oj = oj + 1;
         endo;
   endif;

   qi = qi + 1;
endo;

if maxlag < 2; otres = vec(otresq);
else;          otres = otresp|vec(otresq);
endif;

RETP(otres);ENDP;

