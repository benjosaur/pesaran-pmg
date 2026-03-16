/* ****************************************************************************
*  Pooled Maximum Likelihood & Mean Group Estimation of
*  Dynamic Heterogeneous Panel Data Models
*  -----------------------------------------------------------
*  Written by Yongcheol Shin
*  First version: January 1997
*  Current version: November 1998
*  If you encounter any problems in running this program, send e-mail to
*  Yongcheol.Shin@ed.ac.uk
*  -----------------------------------------------------------
*  See also Pesaran, Shin and Smith, Pooled Mean Group Estimation of
*  Dynamic Heterogeneous Panels, 1998, JASA (forthcoming). This paper
*  is also available on the internet at:
*  http://www.econ.cam.ac.uk/faculty/pesaran/
*  or
*  http://www.ed.ac.uk/~shiny/
*  -----------------------------------------------------------
*
* 1. This program estimates the dynamic heterogenous panel data model,
*    allowing for short-run dynamic heterogeneity but restricting all or
*    a subset of the long-run coefficients to be the same across groups.
*
* 2. Step 1 of the program specifies the subroutines, loads the data,
*    sets the output file, and provides further information to the
*    program such as selection of the dependent variable, selection of
*    the regressors, specification of the group names/group codes, etc.
*    This part of the program needs to be adapted by the user for his/her
*    application.
*
*    Steps 2-10 of the program execute the estimation of the model
*    specified in Step 1. These and all subroutines called by the program
*    (sub0.prc, sub1.prc, sub2.prc, fix.prc, pr.prc, prsm.prc) should
*    not be modified by the user.
*
*    Upon execution of the program you will be prompted for further input
*    on the screen.
*
* 3. This program can also deal with the case where one of the regressors
*    is the relative or reference variable, for example, income relative to
*    U.S. income. Upon execution of the program you will be prompted to
*    choose the appropriate options.
*
* 4. It is recommended to run the program in its full version. But it may
*    also be run in a small version if the user's computer memory or
*    workspace is not large enough.
**************************************************************************** */

NEW;

/* ****************************************************************************
*  The following six files containing sub-procedures MUST be included
*  for running this program.
*  As written, the program assumes that this program and all the
*  subroutines are placed in the same directory, and that the GAUSS
*  configuration file (gauss.cfg) has been modified so that the
*  src_path specification includes the directory in which you
*  have placed the main program and the program files.
*  If you do not want to modify the GAUSS configuration file, you
*  will need to modify the path statements following the INCLUDE
*  commands to indicate in which directory you have placed the
*  subroutines.
**************************************************************************** */

@* common sub-procedures @

#INCLUDE sub0.prc;

@* sub-procedures necessary for PMG estimation when all long-run parameters
** are restricted to be the same across groups @

#INCLUDE sub1.prc;

@* sub-procedures necessary for PMG estimation when a subset of long-run
** parameters is restricted to be the same across groups @

#INCLUDE sub2.prc;

@* sub-procedures necessary for fixed effects estimation @

#INCLUDE fix.prc;

@* sub-procedures necessary for printing the full-version results @

#INCLUDE pr.prc;

@* sub-procedures necessary for printing the results of a small version @

#INCLUDE prsm.prc;

/* ****************************************************************************
*  Step 1: Loading the data and setting the output file
*  ****************************************************************************
*  IMPORTANT!! PLEASE READ BEFORE RUNNING THE PROGRAM
*  ****************************************************************************
*
*  Before running the program, you need to go through the following steps
*  right after these comments:
*
*  1. Choose the way your data is read into GAUSS.
*
*  2. Specify the the output file(s).
*
*  3. Make the necessary transformations of the data or variables.
*
*  4. Select the dependent (Y) and independent variables (X), which vary over
*     time periods and over groups.
*
*  5. IMPORTANT!!
*     If there is a relative or reference variable in X, you MUST put
*     the data for the reference group in the last MAXT rows of this variable,
*     and when selecting the X regressors, you MUST put these data
*     in the last column of the X matrix.
*
*  6. (Optional) Select the independent fixed variables (Z), such as constants
*     and time trend, which vary (at most) over time periods. In most cases
*     you should include at least a constant.
*
*  7. Supply the necessary information on the (maximum) number of time periods,
*     groups, and iterations, and specify a value for missing data.
*
*  The data are assumed to be stored in ASCII files, and are simple
*  matrices containing no alphabetics, but can include missing values.
*
*  [1] Data Option 1 (default and recommended)
*
*  Choose this if you have the data files in a form such that each file
*  contains a MAXT*MAXN-vector of the data of one variable at a time,
*  where MAXT is the maximum of T(1), ..., T(N), and T and N = MAXN refer
*  to the number of time periods and groups, respectively.
*
*  For example, if you have the three data files named energy.dat, output.dat,
*  and price.dat, then the files are to be read into GAUSS as follows:
*
*               load data1[MAXT*N,1] = energy.dat;
*               load data2[MAXT*N,1] = output.dat;
*               load data3[MAXT*N,1] = price.dat;
*
*  [2] Data Option 2 (currently not available)
*
*  This case is relevant if each data file contains the whole data matrix
*  of each group such that there are N files in total. For the i-th group,
*  the data file should be arrayed as [Y(i),X(i)], where Y(i) is a
*  MAXT*1-vector containing the dependent variable, and X(i) is MAXT*v
*  matrix of independent variables excluding intercept. In this case you must
*  load the N data files sequentially.
*
*  For example, if N = 3, there are v indpendent variables, and
*  the data file names are given by ban.dat, ind.dat, and kor.dat, then
*  the files are to be read into GAUSS as follows:
*
*               load dat1[MAXT,v] = ban.dat;
*               load dat2[MAXT,v] = ind.dat;
*               load dat3[MAXT,v] = kor.dat;
**************************************************************************** */

@* Specify the name of output file @

OUTPUT FILE = output.out RESET;

@*  Set the maximum number of iterations for the pooled ML estimation.
**  In case of non-convergence problems, try increasing this number. @

MAXITER = 1000;

@*  Set the value for the missing observations - the default value, 8934567,
**  is the one used in Microfit. But this can be changed by the user. @

MISSING = 8934567;

@* Specify the maximum number of time periods (MAXT) and groups (MAXN). @

MAXT = 34;  @ 1960 to 1993 @
MAXN = 24;  @ OECD countries @

@* Set the data option: the default is data_op = 1. If you want
** to use the data option 2, reset this to 2. @

Data_op = 1;

IF     Data_op eq 1;

   temp = maxt*maxn;

@*  Loading the (raw) data, each variable on N different groups in separate
**  files. There are no restrictions on the number of variables that can be
**  read in at this stage. @

@ * ***************************************** @
@ * JASA: Application 1: Consumption function @
@ * *************************************************
*   For the raw data and the relevant transformations
*   see the attached program called JASADAT1.PRG
*   ************************************************* @

@* Load the variables. @

load lpc[temp,1] = lpc.dat;
   @ natural logarithm of per capita real consumption @

load lndi[temp,1] = lndi.dat;
   @ natural logarithm of per capita real disposable income @

load lpdi[temp,1] = lpdi.dat;
   @ natural logarithm of per capita private real disposable income @

load dp[temp,1] = dp.dat;    @ inflation (deflator) @

@*  Load group names (alphanumeric) or codes (numeric) if any.
**  Its dimension MUST be MAXN*1. If there are no group names/codes,
**  leave this space blank. You MUST call this variable "GRNAME". @

load GRNAME[maxn,1] = J1name.dat; @  Country codes: alphanumeric @

@* Note that there are three types of the variables.
**
** (1) Dependent variabe (Y),
**
** (2) "X" regressors, which are the part of ARDL lag order search process,
**
** (3) "Z" Regressors, which are not subject to ARDL lag search process. @

@* Selecting dependent variable (Y). @

   Y = lpc;

@* Selecting X regressors. @

   X = lndi~dp;

@* Selecting deterministic or fixed regressors (Z).
**
**  In the program deterministic regressors such as constants, time trends,
**  and seasonal dummies are provided as below:
**
**  inpt: constant terms
**  trend: time trends
**  QS = qs1~qs2~qs3: 3 quarterly seasonal dummies
**  MS = ms1~ms2~...~ms11: 11 monthly seasonal dummies
**
**  You may include one or several of the above deterministic regressors
**  in specifying Z.
**
**  IMPORTANT!!
**  Note that the intercept MUST always be placed in the last column of Z.
**
**  If you wish to include another fixed regressor, say one called other,
**  which has been loaded as
**
**       load other[temp,1] = other.dat;
**
**  and trend, quarterly seasonal dummies, and intercept in Z, then put,
**  for example,
**
**       Z = other~trend~QS~inpt;
**
**  and note that the number of variables in Z is 6 in the above case. @

@* Creating the intercept, the linear time trend, and seasonal dummies @

   { inpt, trend, qs, ms } = DETREGR(maxt,maxn);  @ Do not modify this line! @

@* Selecting Z regressors. @

    Z = inpt;

ELSEIF Data_op eq 2;

ENDIF;

/* ****************************************************************************
*  Be sure to check if you have completed the data loading process correctly.
*                              !!IMPORTANT!!
*  From this point onwards, PLEASE do not try to change or modify the program!!
**************************************************************************** */

n = maxn;
nz = cols(z);   @ number of fixed regressors including intercept @

@ Data matrix used in the data transformations @

DAT = Y~X;

/* ****************************************************************************
*  Step 2: Initial information which is needed
*          for the remainder of the program.
**************************************************************************** */

{beg_op, k, Yname, Xname, Zname } = BEG(N,nz);

if nz eq 1; zname = zname[1];
endif;
if beg_op eq 2;
   if nz eq 1; zname = zname[1];
   endif;
   goto SKIP1;
endif;

/* ****************************************************************************
*  STEP3: Some further data handling, if necessary.
*******************************************************************************
*  You may select a subset of the group, transform the raw data to
*  the (cross-section) demeaned ones, and delete missing values.
*  ************************************************************************* */

@*  *********************** @
@*  OPTION for group naming @
@*  *************************
*   grname: nx1 vector of group names
*   gn_op1: option whether or not there is group name or code
*   gn_op2: option whether group name is numeric or alphanumeric @

@ ******* @
  BACKG4:
@ ******* @

"Type 2 if you have already loaded the group names or their numeric codes.
Otherwise type 1.";"";

gn_op1 = con(1,1);

if     gn_op1 eq 1;  gn_op2 = 1;
elseif gn_op1 eq 2;

@ ******* @
  BACKG4Z:
@ ******* @

   "Type 1 if you have loaded the groups' numeric codes.";
   "Type 2 if you have loaded the group names (alphanumeric).";"";

    gn_op2 = con(1,1);
    if     gn_op2 eq 1;
    elseif gn_op2 eq 2;
    else; "WARNING! Try again.";
          "Press any key to continue."; wait; cls;
          goto BACKG4Z;
    endif;

else; "WARNING! Try again.";
      "Press any key to continue."; wait; cls;
      goto BACKG4;
endif;

@* ************************************************ @
@* OPTION for selecting the whole or sub-group data @
@* ************************************************ @

output off;

@ ******* @
  BACKS1:
@ ******* @

"Type 1 if you want to use the whole group data.
Type 2 if you want to use only a sub-group of the data.";"";

sub_op = con(1,1);"";

IF     sub_op eq 1;
ELSEIF sub_op eq 2;
       {dat,n_subgr,na_subgr,n_ex,n_pos,n} = SUBGR(dat,k,n,maxt,grname,gn_op1);
ELSE; "WARNING! Try again.";
      "Press any key to continue."; wait; cls;
      goto BACKS1;
ENDIF;

@* *************************************************************************** @
@* Transform raw data to cross-section demeaned ones and delete missing values @
@* *************************************************************************** @

dmdat = zeros(maxt*N,k+1);
di = 1;
do until di > k+1;
   dtemp = dat[.,di];
   { dmdat[.,di] } = DEMEAN(dtemp,maxt,n,missing);
   di = di + 1;
endo;

@* ************************************************************ @
@* Deleting missing values of the raw data if any, and counting
*  the number of time periods for each group i = 1,..., N. @
@* ************************************************************ @

{ Rdata, Ddata, T } = DELMISS(dat,dmdat,maxt,n,missing); @ {T(1),...,T(N)} @
sumt = sumc(t);                                   @ sum of T(1), ..., T(N) @

cls;
format /rd 3,0;
"** Data loaded and transformed sucessfully. ** ";"";
IF sub_op eq 2;
"** You excluded" n_ex "group(s) from the total of" maxn "groups. **";"";
ENDIF;

@* ******************************************** @
@* OPTION for presence of the relative variable @
@* ******************************************** @

@ ******** @
  BACKREL:
@ ******** @

"You may have included a relative variable (for example, income relative
to U.S. income) in the X regressors.

See also Comment 5 in the data loading step on how to deal with this issue.

Type 1 if the regressors X do not include a relative variable.
Type 2 if the regressors X include such a relative variable.";"";

rel_op = con(1,1);"";

if     rel_op eq 1;  NN = n;
elseif rel_op eq 2;  NN = n - 1;
else; "WARNING! Try again.";
      "Press any key to continue."; wait; cls;
      goto BACKREL;
endif;

@* *************************************************** @
@* OPTION for the full or small version of the program @
@* *************************************************** @

@ ******* @
  BACKSM:
@ ******* @

cls;
"Due to computer memory or workspace limit, you may not be able to run the full
version of the program. In this case you may choose to run the small version of
the program and still obtain the key estimation results without the
accompanying standard errors.

Type 1 if you want to run the full version of the program.
Type 2 if you want to run the small version of the program.";"";

sm_op = con(1,1);

if     sm_op eq 1;
elseif sm_op eq 2;
else; "WARNING! Try again.";
      "Press any key to continue."; wait; cls;
      goto BACKSM;
endif;

@  *************************** @
@* Saving the initial settings @
@  *************************** @

save sub_op t sumt Rdata Ddata n gn_op1 gn_op2 grname;
save rel_op nn sm_op;

if    sub_op eq 1;
else;
      if gn_op1 eq 1;
         save n_ex n_pos n_subgr;
      else;
         save n_ex n_pos n_subgr na_subgr;
      endif;
endif;

/*
***********
*  RESTART:
*  ********
*  If you use the same data set, and initial settings
*  are already created, then restart here.
* ***************************************************
 */
@ ****** @
  SKIP1:
@ ****** @

if beg_op eq 2;
   load sub_op t sumt Rdata Ddata n gn_op1 gn_op2 grname;
   load rel_op nn sm_op;
   if    sub_op eq 2;
         if gn_op1 eq 1;
            load n_ex n_pos n_subgr;
         else;
            load n_ex n_pos n_subgr na_subgr;
         endif;
   endif;
endif;

@* ************************************************************** @
@* OPTION for a choice of the raw or cross-section de-meaned data @
@* ************************************************************** @

@ ******* @
  BACK6:
@ ******* @

"Type 1 if you want to use the raw data set.
Type 2 if you want to use the (cross-section) demeaned data to deal with
the common factor problem.";

dem_op = con(1,1);"";

IF     dem_op eq 1;    Data = Rdata;  @ Raw data @
ELSEIF dem_op eq 2;    Data = Ddata;  @ Demeaned data @
ELSE; "WARNING! Try again.";
      "Press any key to continue."; wait; cls;
      goto BACK6;
ENDIF;

Y = data[.,1];             @ dependent variable @
X = data[.,2:cols(data)];  @ (stochastic) independent variables @

/*
*******************************************************************************
*  Step 4. Selection of the lag orders used in dynamic regressions
*  ---------------------------------------------------------------
*  There are two options of selecting the lag orders.
*
*  (1) Fix the lag orders by yourself.
*  (2) Choose the lag orders by model selection criteria such as
*     AIC, SBC and HQ.
*******************************************************************************
*/
@ ******* @
  LBACK0:
@ ******* @
output off;cls;

if rel_op eq 2; lag_op = 1; goto SKIPLAG;
endif;

"There are two options for selecting the lag orders of the Y and X variables.

Type 1 if you want to choose the lag orders yourself.
Type 2 if you want to choose the lag orders by model selection criteria
such as AIC (Akaike), SBC (Schwarz) and HQ (Hannan and Quinn).";"";

lag_op = con(1,1);"";
if     lag_op eq 1;
elseif lag_op eq 2;
else; "WARNING! Try again.";
      "Press any key to continue."; wait; cls;
      goto LBACK0;
endif;

@ ******** @
  SKIPLAG:
@ ******** @

{PLAG, QLAG, ms_op, maxlag} = LAGSEL(lag_op,Y,X,Z,n,k,T,rel_op);

/*
********************************************************************************
*  Step 5. ML estimation option:
*  -----------------------------
*  There are two options for the pooled ML estimation in Panels.
*
*  (1) Restrict all the long-run parameters to be the same across groups.
*  (2) Restrict only the subset of the long-run parameters to be the same.
********************************************************************************
*/

if k eq 1; ml_op = 1; goto SKIP2;
endif;

output off; cls;
@ ******** @
  MLBACK0:
@ ******** @

"There are two options for pooled maximum likelihood estimation in panels.

Type 1 if you want to restrict all the long-run parameters to be the same
across groups.
Type 2 if you want to restrict only a subset of the long-run parameters
to be the same across groups.";

ML_op = con(1,1);"";
IF     ml_op eq 1;
ELSEIF ml_op eq 2;
       {X1,X2,X1name,X2name,k1,k2,Q1LAG,Q2LAG}
       = MLEOP(X,Xname,k,sumt,QLAG,rel_op);
ELSE; "WARNING! Try again.";
      "Press any key to continue."; wait; cls;
      goto MLBACK0;
ENDIF;

@ ****** @
  SKIP2:
@ ****** @

/*
********************************************************************************
*  Step 6. OLS estimation for each group
*  -------------------------------------
*  Here we compute the followings using the lag orders selected by Step 5:
*
*  ML Estimation option 1:
*  (All long-run parameters restricted to be ths same across groups)
*  -----------------------------------------------------------------

*  OLSRES: OLS estimates of coefficients, standard errors and t-rarios.
*          Results are stored in the order of theta (long-run on X),
*          phi (error correction), beta (short-run on X), dynamic parameters
*          on dy's and dX's and the coeeficients on Z variables for each group.
*
*          Number of parameters in each group is
*          m = k + 1 + k + numsi(i),
*          so the dimension of OLSRES is  mNx3
*
*  OLSDRES: Rbar-suqare, log-likelihood, and diagnostic results for the OLS
*           estimation results. The dimension is 11xN, and the row order is
*           RSQ, RBARSQ, LL, AIC, SBC, HQ, chiSC, chiFF, chiNO, chiHE and sige.
*
*  Theta: Estimates of the LR parameters which are restricted to be the
*         same across groups. But, here they are obtained from each OLS
*         regression, so the dimension is kxN.
*
*  URSLL: Sum of unrestricted log-likelihood for N groups.
*
*  Numsi: Number of other paramters (e.g. on the first differences of y
*         and x, and their lagged values, and Z variables) for each group.
*         This may be different and depends on the selected lag values.
*
*  CFname: Name of coefficients crorresponding to OLSRES.
*
*  ML Estimation option 2:
*  (Subset of long-run parameters restricted to be ths same across groups)
*  -----------------------------------------------------------------------
*  OLSRES: OLS results for the coefficients, standard errors and t-rarios
*          Results are stored in the order of theta1, theta2, phi, beta1,
*          beta2, other parameters including Z variables for each group.
*
*          Number of parameters in each group is
*          m = k1 + k2 + 1 + k1 + k2 + numsi(i),
*          so the dimension of OLSRES is  mNx3.
*
*  OLSDRES: R-suqare, log-likelihood, and diagnostic results for the OLS
*           estimation results. The dimension is 11xN, and the row-order is
*           RSQ, RBARSQ, LL, AIC, SBC, HQ, chiSC, chiFF, chiNO, chiHE and sige.
*
*  Theta1: Estimates of the LR parameters which are restricted to be the
*          same across groups. But, here they are obtained from each OLS
*          regression, so the dimension is k1xN.
*
*  Theta2: Estimates of the LR parameters which are allowed to be different
*          across groups. The dimension is k2xN.
*
*  URSLL, Numsi and CFname are similar as before.
********************************************************************************
*/

IF     ml_op eq 1;
       {OLSRES, OLSDRES, THETA, URSLL, NUMSI, CFname}
       = OLS1(Y,X,Z,PLAG,QLAG,k,NN,T,Yname,Xname,Zname);

ELSEIF ml_op eq 2;
       {OLSRES, OLSDRES, THETA1, THETA2, URSLL, NUMSI, CFname}
       = OLS2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,k1,k2,NN,T,Yname,X1name,X2name,Zname);
ENDIF;

/*
********************************************************************************
*  Step 7. Computation of the (initial) estimates of the LR parameters
*  -------------------------------------------------------------------
*  Here we obtain the initial estimates of the long-run parameters
*  using one of the mean group estimator (MGE), the (static) fixed effects
*  estimator (SFE) or the dynamic fixed effects estimator (DFE).
*
*  Of course, you may give your own initial estimates of the LR parameters.
*
*  These initial estimates are used for the pooled ML estimatiom of
*  by yourself.
********************************************************************************
*/

@* Mean Group estimator of the long-run parameters @

IF     ml_op eq 1;
       MGERES = zeros(k,3);        @ estimate, se and t-ratio @
       MGEdev = zeros(NN,k);
       ii = 1;
       do until ii > k;
         temp = theta[ii,.]';
         { MGERES[ii,.], MGEdev[.,ii] } = AVG(temp);
         ii = ii + 1;
       endo;
       MGECOV = mgedev'mgedev./(NN*(NN-1));  @ Covariance for MGE theta @

ELSEIF ml_op eq 2;
       MGERES1 = zeros(k1,3);        @ estimate, se and t-ratio @
       MGEdev1 = zeros(NN,k1);
       ii = 1;
       do until ii > k1;
         temp = theta1[ii,.]';
         { MGERES1[ii,.], MGEdev1[.,ii] } = AVG(temp);
         ii = ii + 1;
       endo;
       MGECOV = mgedev1'mgedev1./(NN*(NN-1));  @ Covariance for MGE theta1 @
ENDIF;

@* The static fixed effects estimates of the long-run parameters @

IF     ml_op eq 1;
ELSEIF ml_op eq 2;
       X = X1~X2; k = k1 + k2;
       Xname = X1name|X2name;
ENDIF;
{SFERES, SFEIND, RSFERES, SFEDRES, sfeobs, sfet} = SFE(Y,X,Z,n,T,k);

@* The dynamic fixed effects estimates of the long-run parameters @

If lag_op eq 1;
   numot = numsi[1] - 1;  @ number of other parameters excluding intercept @
   IF      ml_op eq 1;
           {DFERES, DFEIND, RDFERES, DFEDRES, dfeobs, dfet}
           = DFE1(Y,X,Z,n,T,k,numot,plag,qlag);

    ELSEIF ml_op eq 2;
          {DFERES, DFEIND, RDFERES, DFEDRES, dfeobs, dfet1, dfet2}
          = DFE2(Y,X1,X2,Z,n,T,k1,k2,numot,plag,q1lag,q2lag);
    ENDIF;
Endif;

@* ************************************************************* @
@* Selection of the initial estimates of the long-run parameters @
@* ************************************************************* @

@ ********* @
  INIBACK1:
@ ********* @

output off; cls;

"** Now select the procedure for obtaining initial estimate(s) of the
    long-run parameter(s) necessary for the pooled maximum likelihood
    estimation of the long-run parameters. **

Type 1 if you want to choose the mean group estimates.
Type 2 if you want to choose the static fixed effects OLS estimates.
Type 3 if you want to choose the dynamic fixed effects OLS estimates.
Type 4 if you want to choose your own initial values.";"";

if lag_op eq 2;
   "When using the model selction criteria to choose the lag orders,
Option 3 is not available.";
endif;

if rel_op eq 2;
   if ml_op eq 2;
      "When using a relative variable and choosing the option that a
subset of the long-run parameters is restricted to be the same across
groups, Option 1 is not available.";"";
   endif;
endif;

ini_op = con(1,1);

if     ini_op eq 1;
elseif ini_op eq 2;
elseif ini_op eq 3;
elseif ini_op eq 4;
else; "WARNING! Try again.";
      "Press any key to continue."; wait; goto INIBACK1;
endif;

if lag_op eq 2;
   if ini_op eq 3;
      "WARNING! The dynamic fixed effects estimates are not available! Try
again."
      "Press any key to continue."; wait; goto INIBACK1;
   endif;
endif;"";

if rel_op eq 2;
   if ml_op eq 2;
      if ini_op eq 1;
         "WARNING! The mean group estimates are not available! Try again.";
         "Press any key to continue."; wait; goto INIBACK1;
      endif;
   endif;
endif;

IF     ml_op eq 1;
   if     ini_op eq 1;   ini_LR = mgeres[.,1];
   elseif ini_op eq 2;   ini_LR = sfet;
   elseif ini_op eq 3;   ini_LR = dfet;
   elseif ini_op eq 4;   format /rd 3,0;
      "Type your own" k "initial estimates.";
       ini_LR = con(k,1);"";
   endif;

ELSEIF ml_op eq 2;
   if     ini_op eq 1;  ini_LR1 = mgeres1[.,1]; ini_LR2 = theta2;
   elseif ini_op eq 2;  ini_LR1 = sfet[1:k1];   ini_LR2 = sfet[k1+1:k1+k2];
   elseif ini_op eq 3;  ini_LR1 = dfet1;        ini_LR2 = dfet2;
   elseif ini_op eq 4;   format /rd 3,0;
      "Type your own" k1 "initial estimates of Theta1.";
       ini_LR1 = con(k1,1);"";
      "Type your own" k2 "initial estimates of Theta2.";
       ini_LR2 = con(k2,1);"";
   endif;

ENDIF;

/*
************************************************************
*  Step 8. Pooled ML estimation of the long-run parameters
************************************************************
*/

@ ********  @
  PMGBACK:
@ ********  @

cls; output off;

"** There are two different numerical algorithms. **

Type 1 if you want to use the Back-Substitution method which uses only
the first derivative of the log-likelihood function.

Type 2 if you want to use the Newton-Raphson method which uses both
the first and the second derivative of the log-likelihood function.";"";

if     sm_op eq 1;  nu_op = con(1,1);
elseif sm_op eq 2;  nu_op = 1;
   "If you have chosen the small version of the program, the Newton-Raphson
method is not available, so the Back-Substitution method is selected
automatically by the program.";
   "Press any key to continue."; wait; cls;
endif;

if     nu_op eq 1; goto BSA;
elseif nu_op eq 2; goto NRA;
ELSE; "WARNING! No option is available. Try again.";
      "Press any key to continue."; wait; cls;
      goto PMGBACK;
endif;

/*
*******************************************************************************
*  Step 8A. Pooled ML estimation of the LR parameters by BS algorithm
*  -------------------------------------------------------------------
*  After obtaining the converged pooled ML estimates of long-run parameters,
*  we then obtain other parameters including the coefficients on Z.
*
*  Then, we also compute diagnostic statistics and so on for each group
*  using the lag orders selected by Step 5 in the individual OLS regressions.
*******************************************************************************
*/

@ ****  @
  BSA:
@ ****  @

cflag = 0;   @ Not converged if 0, and converged if 1. @

IF     ml_op eq 1; @ ALL LR PARA THE SAME @

   oldtheta = ini_LR;   @ initial estimates of Theta @

   @ Tolerance level used in convergence criteria @

   tol1 = 10^(-4);     @ check if Like(i) - Like(i-1) < tol1 @
   tol2 = 10^(-4)*k;   @ check if |Theta1(i)- Theta1(i-1)| < tol2 @

   iter = 1;
   do until iter > maxiter;

      {THETA,RSLL} = PMLEBS1(Y,X,Z,PLAG,QLAG,k,N,T,OLDTHETA);

      if iter eq 1; OLDRSLL = RSLL; oldtheta = theta;
                    goto repeatA1;
      endif;

      @ check if the convergence occurred @

      convlike = RSLL - oldRSLL;
      convt = abs(theta - oldtheta);
      convtr = sumc(convt);
      convcr = convlike|convtr;

      cls; screen on; output off; format /rd 12,8;
      "Number of Iterations:";; iter;
      "Restricted Sum of Log-Likelihood Functions: ";;RSLL;
      "RSLL - RSLL(-1):";;convlike;
      "ABS(theta - oldtheta):";;convt'; theta;

      if    convlike < tol1;
         if    convtr < tol2;  cflag = 1;
               break;  @ convergence occurred @
         endif;
      endif;

      OLDRSLL = RSLL; oldtheta = theta;

@ ********* @
  REPEATA1:
@ ********* @

      iter = iter + 1;"";
endo;

ELSEIF ml_op eq 2;   @ SUBSET OF LR PARA ARE THE SAME @

   @ Initial estimates of Theta1 and Theta2 @

   oldthta1 = ini_LR1; oldthta2 = ini_LR2;

   @ Tolerance level used in checking the convergence criteria @

   tol1 = 10^(-4);     @ check if LL(i) - LL(i-1) < tol1 @
   tol2 = 10^(-4)*k1;  @ check if |Theta1(i)- Theta1(i-1)| < tol2 @

   iter = 1;
   do until iter > maxiter;

      {THETA1,THETA2,RSLL} =
    PMLEBS2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,k1,k2,N,T,OLDTHTA1,OLDTHTA2,ini_op,iter);

      if iter eq 1; OLDRSLL = RSLL; oldthta1 = theta1; oldthta2 = theta2;
                    goto repeatA2;
      endif;

   @ check if the convergence occurred @

      convlike = RSLL - oldRSLL;
      convt1 = abs(theta1 - oldthta1);
      convtr = sumc(convt1);
      convcr = convlike|convtr;

      cls; screen on; output off; format /rd 12,8;
      "Number of Iterations:";; iter;
      "Restricted Sum of Log-Likelihood Functions: ";;RSLL;
      "RSLL - RSLL(-1):";;convlike;
      "ABS(theta1 - oldtheta1):";;convt1'; theta1;

      if    convlike  < tol1;
         if    convtr < tol2; cflag = 1;
               break;  @ convergence occurred @
         endif;
      endif;

      OLDRSLL = RSLL;   oldthta1 = theta1;   oldthta2 = theta2;

@ ********* @
  REPEATA2:
@ ********* @

      iter = iter + 1;"";
   endo;

ENDIF;

GOTO NUMEND;

/*
**********************************************************************
*  Step 8B. Pooled ML estimation of the LR parameters by NR algorithm
*  ------------------------------------------------------------------
*  See also 8A
**********************************************************************
 */

@ ****  @
  NRA:
@ ****  @

cflag = 0;   @ If 0, not converged, if 1, converged. @

IF     ml_op eq 1;

   oldtheta = ini_LR;      @ theta(j-1) @

   { oldphi } = INIPHI1(Y,X,Z,PLAG,QLAG,k,N,T,OLDTHETA); @ phi(j-1) @

   @ Tolerance level used in convergence criteria @

   tol1 = 10^(-4);     @ check if Like(i) - Like(i-1) < tol1 @
   tol2 = 10^(-4)*k;   @ check if |Theta1(i)- Theta1(i-1)| < tol2 @

   iter = 1;
   do until iter > maxiter;

      {THETA,PHI,RSLL}
      = PMLENR1(Y,X,Z,PLAG,QLAG,k,N,T,OLDTHETA,OLDPHI);

      if iter eq 1; OLDRSLL = RSLL; oldtheta = theta; oldphi = phi;
                    goto repeatB1;
      endif;

   @ check if the convergence occurred @

      convlike = RSLL - oldRSLL;
      convt = abs(theta - oldtheta);
      convtp = abs(phi - oldphi);
      convtr = sumc(convt);
      convtrp = sumc(convtp);
      convg = convtr + convtrp;
      convcr = convlike|convtr;

      cls; screen on; output off; format /rd 12,8;
      "Number of Iterations:";; iter;
      "Restricted Sum of Log-Likelihood Functions: ";;RSLL;
      "RSLL - RSLL(-1):";;convlike;
      "ABS(theta - oldtheta):";;convt'; theta;

      if    convlike < tol1;
         if    convtr < tol2;  cflag = 1;
               break;  @ convergence occurred @
         endif;
      endif;

      OLDRSLL = RSLL; oldtheta = theta; oldphi = phi;

@ ********* @
  REPEATB1:
@ ********* @

      iter = iter + 1;"";
endo;

ELSEIF ml_op eq 2;

   @ initial estimates of Theta1 and Theta2 @

   oldthta1 = ini_LR1; oldthta2 = ini_LR2;

   { oldphi } = INIPHI2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,k1,k2,N,T,
                        OLDTHTA1,OLDTHTA2,ini_op);

   @ Tolerance level used in the convergence criteria @

   tol1 = 10^(-4);     @ check if LL(i) - LL(i-1) < tol1 @
   tol2 = 10^(-4)*k1;  @ check if |Theta1(i)- Theta1(i-1)| < tol2 @

   iter = 1;
   do until iter > maxiter;

      {THETA1,THETA2,PHI,RSLL} = PMLENR2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,
                     k1,k2,N,T,OLDTHTA1,OLDTHTA2,OLDPHI,ini_op,iter);

      if iter eq 1; OLDRSLL = RSLL;
                    oldthta1 = theta1; oldthta2 = theta2; oldphi = phi;
                    goto repeatB2;
      endif;

   @ check if the convergence occurred @

      convlike = RSLL - oldRSLL;
      convt1 = abs(theta1 - oldthta1);
      convtr = sumc(convt1);
      convcr = convlike|convtr;

      cls; screen on; output off; format /rd 12,8;
      "Number of Iterations:";; iter;
      "Restricted Sum of Log Likelihood Functions: ";;RSLL;
      "RSLL - RSLL(-1):";;convlike;
      "ABS(theta1 - oldtheta1):";;convt1'; theta1;

      if    convlike  < tol1;
         if    convtr < tol2; cflag = 1;
               break;  @ convergence occurred @
         endif;
      endif;

      OLDRSLL = RSLL; oldthta1 = theta1; oldthta2 = theta2; oldphi = phi;

@ ********* @
  REPEATB2:
@ ********* @

      iter = iter + 1;"";
   endo;

ENDIF;

@ ******* @
  NUMEND:
@ ******* @

/*
****************************************************************************
*  Step 9. The pooled ML estimation of all parameters
*  --------------------------------------------------
*  Here we obtain the pooled ML estimation results of all parameters
*  with their standard errors obtained from the information matrix.
*
*  ML estimation option 1
*  ----------------------
*  PMLESRES: Pooled ML estimation results for the coefficients, standard
*            errors and t-rarios. Results are stored in the order of theta,
*            phi(1), ..., phi(N), beta(1), ..., beta(N),
*            si(1), ..., si(N). So, the dimension is
*            (k + n + k*n +[numsi(1)+...+numsi(N)]) by 3.
*
*  PMLEDRES: R-suqare, LL, and diagnostic results for the PML estimation
*            results. The dimension is 11xN, and the order is RSQ, RBARSQ,
*            LL, AIC, SBC, HQ, chiSC, chiFF, chiNO, chiHE and sige.
*
*  RSLL: Sum of restricted log-likelihood for N groups.
*
*  ML estimation option 2
*  ----------------------
*  PMLESRES: Pooled ML estimation results for the coefficients, standard
*            errors and t-rarios. Results are stored in the order of theta1,
*            theta2(1), ..., theta2(N), phi(1), ..., phi(N),
*            beta1(1), ..., beta1(N), beta2(1), ..., beta2(N),
*            si(1), ..., si(N). So, the dimension is
*            (k1 + k2*n + n + k1*n + k2*n + [numsi(1)+...+numsi(N)])x3.
*
*  PMLEDRES: R-suqare, LL, and diagnostic results for the PML estimation
*            results. The dimension is 11xN, and the order is RSQ, RBARSQ,
*            LL, AIC, SBC, HQ, chiSC, chiFF, chiNO, chiHE and sige.
*
*  RSLL: Sum of restricted log-likelihood for N groups.
*******************************************************************************
*/

If sm_op eq 2; goto SKIPSM0;  @ SMALL VERSION @
endif;

IF     ml_op eq 1;
       {PMLERES, PMLEDRES, RSLL, PMLECOV}
       = PMLESR1(Y,X,Z,PLAG,QLAG,k,NN,T,NUMSI,THETA);
ELSEIF ml_op eq 2;
       {PMLERES, PMLEDRES, RSLL, PMLECOV}
       = PMLESR2(Y,X1,X2,Z,PLAG,Q1LAG,Q2LAG,k1,k2,NN,T,NUMSI,THETA1,THETA2);
ENDIF;

/*
************************************
*  Step 10. Printing the results
************************************
*/

@ ******** @
  SKIPSM0:
@ ******** @

output on;

@* ************************************** @
@* Printing the options you have selected @
@* ************************************** @

format /rd 3,0;

@ ********************** @
@ Printing date and time @
@ ********************** @

xdate = date;
xdateout = xdate[1:3];
xtime = time;
xtimeout = xtime[1:2];
"Created at " xtime[1] ":" xtime[2] "on " xdate[3] "-" xdate[2] "- " xdate[1];
"";"";


format /rd 3,0;

"****************************************************************************";
"The following information and options have been used for the PMG estimation:";
"****************************************************************************";
"";

if ml_op eq 1;
   "There are" k "stochastic regressors (X) and" cols(Z) "deterministic
regressors (Z).";
elseif ml_op eq 1;
   "There are" k1+k2 "stochastic regressors (X) and" cols(Z) "deterministic
regressors (Z).";
endif;

"";
"The maximum number of time periods and groups are:" Maxt~Maxn;

"";
"The maximum number of iterations for the pooled maximum likelihood estimation
is set to: " maxiter;

"";
"The value for the missing observations is set to:"; missing;

"";
if     gn_op1 eq 1;
       "Group names are given in numeric order.";
elseif gn_op1 eq 2;
    if     gn_op2 eq 1;
           "Group names are given as numeric codes.";
    elseif gn_op2 eq 2;
           "Group names are given as alphanumerics.";
    endif;
endif;

"";
IF     sub_op eq 1; "The whole group data is selected.";
ELSEIF sub_op eq 2;
   "A subgroup of the data is selected.";
   "You excluded" n_ex "group(s) from the total of" maxn "groups.";
ENDIF;

if rel_op eq 2;"";
   "You have included a relative variable as one of the regressors (X).";
endif;

if sm_op eq 2;"";
   "You have selected the small version of the program. Only the summary
estimation results without standard errors will be provided.

If you wish to have more detailed results and there is sufficient workspace
on your computer, use the full version of the program.";
endif;

"";
If dem_op eq 1;
   "The data you have used are the raw ones. If there is a common
factor problem, use demeaned data.";

ELSEIf dem_op eq 2;
   "The data you have used are the cross-section demeaned ones.";
ENDIF;

"";
if     lag_op eq 1; "The fixed lag specification has been selected.";
elseif lag_op eq 2;
   if     ms_op eq 1;
      "AIC (Akaike) has been used to select the lag orders for each group.";
   elseif ms_op eq 2;
      "SBC (Schwarz) has been used to select the lag orders for each group.";
   elseif ms_op eq 3;
 "HQ (Hannan & Quinn) has been used to select the lag orders for each group.";
   endif;
   "Notice here that dynamic fixed effects OLS estimates are not available.";
endif;

"";
IF     ml_op eq 1;
   "All the long-run parameters have been restricted to be the same
across groups.";
ELSEIF ml_op eq 2;
   "Only the long-run parameters on the following subset of the X regressors
have been restricted to be the same across groups:";
   call printfm(X1name,0,"s"~8~8);"";
ENDIF;

"";
if ini_op eq 1;

"The mean group estimates have been used as initial estimate(s) of the long-run
parameter(s) for the pooled maximum likelihood estimation.";

elseif ini_op eq 2;

"The static fixed effects OLS estimates have been used as initial estimate(s)
of the long-run parameter(s) for the pooled maximum likelihood estimation.";

elseif ini_op eq 3;

"Your own values have been used as initial estimate(s) of the long-run
parameter(s) for the pooled maximum likelihood estimation.";

endif;

if rel_op eq 2;
   if ml_op eq 2;"";

      "When using a relative variable and choosing the option that
a subset of the long-run parameters is restricted to be the same across
groups, the mean group estimates cannot be used as the initial estimates.";

   endif;
endif;

"";
if     nu_op eq 1;
   "You have chosen to use the Back-Substitution method which uses only
the first derivative of the log-likelihood function.";
elseif nu_op eq 2;
   "You have chosen to use the Newton-Raphson method which uses both
the first and the second derivative of the log-likelihood function.";
endif;

if sm_op eq 2;"";
   "When you choose the small version of the program, the Newton-Raphson
method which uses both the first and the second derivative of
the log-likelihood function is not available.";
endif;

"";
"Detailed estimation results follow.";"";

"";

If sm_op eq 2; goto SKIPSM;  @ SMALL VERSION @
endif;

@* ********************************************** @
@* Print the summary results for the full version @
@* ********************************************** @

if     sub_op eq 1;  n_pos = 0;
endif;

if     ml_op eq 1;  n_LR = k;
                    Alag = plag~qlag;
                    PMGTHETA = THETA;
                    if lag_op eq 1;
                       COMPRES = Xname~MGERES~SFERES[1:k,.]~DFERES[1:k,.];
                    else;
                       COMPRES = Xname~MGERES~SFERES[1:k,.];
                    endif;
elseif ml_op eq 2;  n_LR = k1;
                    Alag = plag~q1lag~q2lag;
                    PMGTHETA = THETA1;
                    if lag_op eq 1;
                       COMPRES = X1name~MGERES1~SFERES[1:k1,.]~DFERES[1:k1,.];
                    else;
                       COMPRES = X1name~MGERES1~SFERES[1:k1,.];
                    endif;
endif;

page = 1;
format /rd 3,0;
"PAGE "  page;"";

{ Tobs }= PRSUM(cflag,iter,sub_op,n_pos,dem_op,ml_op,ini_op,RSLL,URSLL,n_LR,
          n,T,ms_op,Alag,COMPRES,maxlag,nu_op,PMGtheta,sm_op);

@ Print the ML and OLS estimation results for each group @

if sub_op eq 1; n_subgr = 0; na_subgr = 0;
endif;
if gn_op1 eq 1; grname = 0; na_subgr = 0;
endif;

IF     ml_op eq 1;
       {page,PMLESUM} = PRGROUP1(OLSRES,PMLERES,CFNAME,OLSDRES,PMLEDRES,maxlag,
    k,nn,T,Tobs,numsi,sub_op,n_subgr,grname,na_subgr,gn_op1,gn_op2,yname,page);

ELSEIF ml_op eq 2;
       {page,PMLESUM} = PRGROUP2(OLSRES,PMLERES,CFNAME,OLSDRES,PMLEDRES,maxlag,
 k1,k2,nn,T,Tobs,numsi,sub_op,n_subgr,grname,na_subgr,gn_op1,gn_op2,yname,page);

ENDIF;

@ Print the fixed effects estimation results @

""; page = page + 1; format /rd 3,0;
"PAGE "  page;"";

if ml_op eq 1;

   PRFIX(Yname,Xname,Zname,SFERES,RSFERES,SFEDRES,sfeobs);

   if lag_op eq 1;
      pp = plag[1];      kk = k;
      PRFIX1(CFNAME,Yname,DFERES,RDFERES,DFEDRES,dfeobs,kk,numot,pp);
   endif;

elseif ml_op eq 2;

   XXname = X1name|X2name;
   PRFIX(Yname,XXname,Zname,SFERES,RSFERES,SFEDRES,sfeobs);

   kk = k1 + k2;
   if lag_op eq 1;
      pp = plag[1];
      PRFIX2(CFNAME,Yname,DFERES,RDFERES,DFEDRES,dfeobs,kk,numot,pp,k1,k2);
   endif;
endif;

@ Print the average of all the OLS and MLE parameters @

""; page = page + 1; format /rd 3,0;
"PAGE "  page;"";

if ml_op eq 1;  kk = k;
   PRAVG1(ms_op,OLSRES,PMLERES,nn,numsi,kk,plag,qlag,
          CFname,maxlag,yname,xname,zname,MGECOV,PMLECOV);

elseif ml_op eq 2;
   PRAVG2(ms_op,OLSRES,PMLERES,nn,numsi,k1,k2,plag,q1lag,q2lag,CFname,maxlag,
          yname,x1name,x2name,zname,MGECOV,PMLECOV);
endif;

@ Print the (summary) OLS estimation results of the LR parameters
  with RBARSQ, LL, ... @

""; page = page + 1; format /rd 3,0;"PAGE "  page;"";

if ml_op eq 1;   XXname = Xname;
else;            XXname = X1name|X2name;
endif;

"                     OLS VERSION";"";

PRSUMLR(Yname,XXname,ms_op,OLSRES,OLSDRES,nn,numsi,
        gn_op1,gn_op2,kk,grname,sub_op,n_subgr,na_subgr,page);

If sm_op eq 1;  @ Full Version @

""; page = page + 1; format /rd 3,0;"PAGE "  page;"";
"                     PMLE VERSION";"";

PRSUMLR(Yname,XXname,ms_op,PMLESUM,PMLEDRES,nn,numsi,
        gn_op1,gn_op2,kk,grname,sub_op,n_subgr,na_subgr,page);

endif;

goto PRend;

@ ******* @
  SKIPSM:
@ ******* @

output on;

@ Print the summary results for the small version @

if     sub_op eq 1;  n_pos = 0;
endif;

if     ml_op eq 1;  n_LR = k;
                    Alag = plag~qlag;
                    PMGTHETA = THETA;
                    if lag_op eq 1;
                       COMPRES = Xname~MGERES~SFERES[1:k,.]~DFERES[1:k,.];
                    else;
                       COMPRES = Xname~MGERES~SFERES[1:k,.];
                    endif;
elseif ml_op eq 2;  n_LR = k1;
                    Alag = plag~q1lag~q2lag;
                    PMGTHETA = THETA1;
                    if lag_op eq 1;
                       COMPRES = X1name~MGERES1~SFERES[1:k1,.]~DFERES[1:k1,.];
                    else;
                       COMPRES = X1name~MGERES1~SFERES[1:k1,.];
                    endif;
endif;

page = 1;
format /rd 3,0;
"PAGE "  page;"";
{ Tobs }= PRSUM(cflag,iter,sub_op,n_pos,dem_op,ml_op,ini_op,RSLL,URSLL,n_LR,
          n,T,ms_op,Alag,COMPRES,maxlag,nu_op,PMGtheta,sm_op);

""; page = page + 1; format /rd 3,0;
"PAGE "  page;"";

@ Print the fixed effects estimation results @

if ml_op eq 1;

   PRFIX(Yname,Xname,Zname,SFERES,RSFERES,SFEDRES,sfeobs);

   if lag_op eq 1;
      pp = plag[1];      kk = k;
      PRFIX1(CFNAME,Yname,DFERES,RDFERES,DFEDRES,dfeobs,kk,numot,pp);
   endif;

elseif ml_op eq 2;

   XXname = X1name|X2name;
   PRFIX(Yname,XXname,Zname,SFERES,RSFERES,SFEDRES,sfeobs);

   kk = k1 + k2;
   if lag_op eq 1;
      pp = plag[1];
      PRFIX2(CFNAME,Yname,DFERES,RDFERES,DFEDRES,dfeobs,kk,numot,pp,k1,k2);
   endif;
endif;

if sub_op eq 1; n_subgr = 0; na_subgr = 0;
endif;
if gn_op1 eq 1; grname = 0; na_subgr = 0;
endif;

@ Print the MG Estimates @

""; page = page + 1; format /rd 3,0;
"PAGE "  page;"";

if ml_op eq 1;  kk = k;
   PRAVG1SM(ms_op,OLSRES,nn,numsi,kk,plag,qlag,
          CFname,maxlag,yname,xname,zname);

elseif ml_op eq 2;
   PRAVG2SM(ms_op,OLSRES,nn,numsi,k1,k2,plag,q1lag,q2lag,CFname,maxlag,
          yname,x1name,x2name,zname);
endif;

@ Print the (summary) individual OLS estimation results of LR parameters
  and ECM coefficients with diagnostics such as RBARSQ, LL, ... @

""; page = page + 1; format /rd 3,0;
"PAGE "  page;"";

if ml_op eq 1;   XXname = Xname; kk = k;
else;            XXname = X1name|X2name; kk = k1 + k2;
endif;

PRSUMLR(Yname,XXname,ms_op,OLSRES,OLSDRES,nn,numsi,
        gn_op1,gn_op2,kk,grname,sub_op,n_subgr,na_subgr,page);

@ ****** @
  PRend:
@ ****** @

END;





