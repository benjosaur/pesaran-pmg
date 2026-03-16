/* ********************************************************************* */
   PROC (0) = PRAVG1SM(ms_op,OLSRES,n,numsi,kk,plag,qlag,CFname,maxlag,
                       yname,xname,zname);
/* **********************************************************
   Print the mean group estimates of the LR and SR parameters
   ********************************************************** */

Local ini,fini,oini,ofini,msi,otmp,omatrix,pmatrix,Otemp,Ptemp,
      pini,pfini,bini,bfini,temp1,temp2,temp3,temp4,matrix1,matrix2,matrix3,
      matrix4,OARES,PARES,npara,pp,coname,prout,maxp,maxq,maxpara,sumq,
      otmp0,otmp1,otmp2,otmp3,otmp4,temp5,temp6,temp7,maxot,n_ot,pot,qot,
      junk,thetares,prtmp1,prtmp2,prtmp3,prtmp4;

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

   @ printing the results @

   otmp = kk + 1 + kk + numsi[1];
   coname = cfname[1:otmp];
   prout = coname~OARES;

"          ****************************";
"                  MGE Estimates";
"          ****************************";
"             Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";

"            Coef.     St. Er.  t-ratio";"";

   if plag eq 0;  @ if all phi's are zero @

      "Long-Run Coefficients";"";
      prtmp1 = prout[1:kk,.];
      call printfm(prtmp1,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";
      if rows(prtmp1) eq 1; "";
      endif;

      "Error Correction Coefficients";"";
      prtmp2 = prout[kk+1,1:2]~"NA"~"NA";
      call printfm(prtmp2,0~1~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10);"";"";

      "Short-Run Coefficients";"";
      prtmp3 = prout[kk+2:rows(prout),.];
      call printfm(prtmp3,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";"";

   else;
      PR1SM(prout,kk);
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

   @ printing the results @

    {COname} = MXCNAME1(maxlag,kk,yname,xname,zname);
   prout = coname~OARES;

"          ****************************";
"                  MGE Estimates";
"          ****************************";
"             Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";"";

"            Coef.     St. Er.  t-ratio";"";

   if plag eq 0;  @ if all phi's are zero @

      "Long-Run Coefficients";"";
      prtmp1 = prout[1:kk,.];
      call printfm(prtmp1,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";
      if rows(prtmp1) eq 1; "";
      endif;

      "Error Correction Coefficients";"";
      prtmp2 = prout[kk+1,1:2]~"NA"~"NA";
      call printfm(prtmp2,0~1~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10);
                   "";"";

      "Short-Run Coefficients";"";
      prtmp3 = prout[kk+2:rows(prout),.];
      call printfm(prtmp3,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";

   else;
      PR1SM(prout,kk);
   endif;

ENDIF;

RETP;ENDP;

/* **************************** */
   PROC (0) = PR1SM(prout,k);
/* **************************** */

Local prtmp;

prtmp = prout[1:k,.];
"Long-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";
if rows(prtmp) eq 1; "";
endif;"";

prtmp = prout[k+1,.];
"Error Correction Coefficients";"";
if prtmp[2] eq -1;  @ if phi = -1 @
   prtmp = prtmp[1]~prtmp[2]~"NA"~"NA";
   call printfm(PRtmp,0~1~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10);"";"";
else;
   call printfm(PRtmp,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";
endif;"";

prtmp = prout[k+2:rows(prout),.];
"Short-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);
if rows(prtmp) eq 1; "";
endif;"";

RETP;ENDP;

/* ************************************************************************** */
   PROC (0) =
   PRAVG2SM(ms_op,OLSRES,n,numsi,k1,k2,plag,q1lag,q2lag,CFname,maxlag,
          yname,x1name,x2name,zname);
/* ***************************************************************************
   Print the mean group estimates of the LR and SR parameters
   for the ML option case 2
   ********************************************************** */

Local ini,fini,oini,ofini,msi,otmp,omatrix,pmatrix,Otemp,Ptemp,
      pini,pfini,b1ini,b1fini,temp1,temp2,temp3,temp4,matrix1,matrix2,matrix3,
      t2ini,t2fini,b2ini,b2fini,matrix4,matrix5,matrix6,temp5,temp6,temp7,
      temp8,temp9,OARES,PARES,npara,pp,coname,prout,maxp,maxq,maxpara,sumq,
      otmp0,otmp1,otmp2,otmp3,otmp4,maxot,n_ot,pot,qot,
      junk,thetares,kk,prtmp1,prtmp2,prtmp3,prtmp4;

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

   @ printing the results @

   otmp = kk + 1 + kk + numsi[1];
   coname = cfname[1:otmp];
   prout = coname~OARES;

"          ****************************";
"               MGE Estimates";
"          ****************************";
"             Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";
"            Coef.     St. Er.  t-ratio";"";

   if plag eq 0;  @ if all phi's are zero @

      "Long-Run Coefficients Restricted to be the Same Across All Groups";"";
      prtmp1 = prout[1:k1,.];
      call printfm(prtmp1,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";
      if rows(prtmp1) eq 1; "";
      endif;

     "Unrestricted Long-Run Coefficients";"";
      prtmp1 = prout[k1+1:k1+k2,.];
      call printfm(prtmp1,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";
      if rows(prtmp1) eq 1; "";
      endif;

      "Error Correction Coefficients";"";
      prtmp2 = prout[kk+1,1:2]~"NA"~"NA";
      call printfm(prtmp2,0~1~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10);
                   "";"";

      "Short-Run Coefficients";"";
      prtmp3 = prout[kk+2:rows(prout),.];
      call printfm(prtmp3,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";

   else;
      PR2SM(prout,k1,k2);
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

   @ printing the results @

    {COname} = MXCNAME2(maxlag,k1,k2,yname,x1name,x2name,zname);
    prout = coname~OARES;

"          ****************************";
"                   MGE Estimates";
"          ****************************";
"             Dependent variable: ";;
call printfm(Yname,0,"s"~4~4);"";"";
"            Coef.     St. Er.  t-ratio      Coef.     St. Er.  t-ratio";"";

   if plag eq 0;  @ if all phi's are zero @

      "Long-Run Coefficients Restricted to be the Same Across All Groups";"";
      prtmp1 = prout[1:k1,.];
      call printfm(prtmp1,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";
      if rows(prtmp1) eq 1; "";
      endif;"";

      "Unrestricted Long-Run Coefficients";"";
      prtmp1 = prout[k1+1:k1+k2,.];
      call printfm(prtmp1,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);
      if rows(prtmp1) eq 1; "";
      endif;"";

      "Error Correction Coefficients";"";
      prtmp2 = prout[kk+1,1:2]~"NA"~"NA";
      call printfm(prtmp2,0~1~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10);
                   "";"";

      "Short-Run Coefficients";"";
      prtmp3 = prout[kk+2:rows(prout),.];
      call printfm(prtmp3,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";

   else;
      PR2SM(prout,k1,k2);
   endif;

ENDIF;

RETP;ENDP;

/* ******************************** */
   PROC (0) = PR2SM(prout,k1,k2);
/* ******************************** */

Local prtmp;

prtmp = prout[1:k1,.];
"Long-Run Coefficients Restricted to be the Same Across All Groups";"";
call printfm(PRtmp,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);
if rows(prtmp) eq 1; "";
endif;"";

prtmp = prout[k1+1:k1+k2,.];
"Unrestricted Long-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);
if rows(prtmp) eq 1; "";
endif;"";

prtmp = prout[k1+k2+1,.];
"Error Correction Coefficients";"";
if prtmp[2] eq -1;  @ if phi = -1 @
   prtmp = prtmp[1]~prtmp[2]~"NA"~"NA";
   call printfm(PRtmp,0~1~0~0,"s"~8~8|"lf"~10~4|"s"~10~10|"s"~10~10);"";
else;
   call printfm(PRtmp,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);"";
endif;"";

prtmp = prout[k1+k2+2:rows(prout),.];
"Short-Run Coefficients";"";
call printfm(PRtmp,0~1~1~1,"s"~8~8|"lf"~10~4|"lf"~10~4|"lf"~10~4);
if rows(prtmp) eq 1; "";
endif;"";

RETP;ENDP;

