@ * ********************************************* @
@ * JASA: Application 1: Consumption function     @
@ * ********************************************* @
@ * Loading the raw data and data transformations @
@ * ********************************************* @

new;

@* Specify the name of the output file, if necessary. @

OUTPUT FILE = DATA.OUT RESET;

@* Set the value for the missing observations - the default value, 8934567,
** is the one used in Microfit. This can be changed by the user. @

MISSING = 8934567;

@* Specify the maximum number of time periods (MAXT) and groups (MAXN). @

MAXT = 34;  @ maximum number of time periods @
MAXN = 24;  @ maximum number of groups @

temp = maxt*maxn;

@* Loading the (raw) data, each variable on N different groups in separate
** files. There are no restrictions on the number of variables that can be
** read in at this stage. @

@* Load the raw data. @

load cp[temp,1] = cp.dat;
   @ private final consumption expenditure (current) @

load rcp[temp,1] = rcp.dat;
   @ real private final consumption expenditure at 1989 prices @

load ndi[temp,1] = ndi.dat;
   @ national disposable income (current) @

load pdi[temp,1] = pdi.dat;
   @ private disposable labor income (current) @

load popul[temp,1] = pop.dat;  @ population data @

@* Data transformation I @

lpc    = ln(rcp./popul);  @ natural log of per capita real consumption @
prcon  = cp./rcp;         @ consumption price deflator @
prcon1 = lag(prcon);      @ lagged price deflator @
lndi = ln(ndi./(prcon.*popul));
   @ natural logarithm of per capita real disposable income @
lpdi = ln(pdi./(prcon.*popul));
   @ natural logarithm of per capita private real disposable income @

@* Checking for missing values @

lpc    = chkmiss(lpc,rcp,popul,temp,missing);
prcon  = chkmiss(prcon,cp,rcp,temp,missing);
prcon1 = chkmiss1(prcon1,maxt,maxn,missing);
lndi   = chkmiss2(lndi,ndi,popul,prcon,temp,missing);
lpdi   = chkmiss2(lpdi,pdi,popul,prcon,temp,missing);

@* Data transformation II @

dp = ln(prcon) - ln(prcon1);  @ inflation (deflator) @
dp = chkmiss(dp,prcon,prcon1,temp,missing);

@* Saving the transformed data (ASCII) to be used in the application @

format /rd 15,6;
OUTPUT FILE = lpc.dat RESET;  lpc;
OUTPUT FILE = lndi.dat RESET; lndi;
OUTPUT FILE = lpdi.dat RESET; lpdi;
OUTPUT FILE = dp.dat RESET;   dp;

end;

/* ************************************************* */
    proc (1) = chkmiss(chk,chk1,chk2,temp,missing);
/* ************************************************* */
local ii,chk1,chk2;

ii = 1;
do until ii > temp;

   if chk1[ii] eq missing; chk[ii] = missing;
   endif;
   if chk2[ii] eq missing; chk[ii] = missing;
   endif;

   ii = ii + 1;
endo;


RETP(CHK);ENDP;

/* ********************************************* */
   proc (1) = chkmiss1(chk,maxt,maxn,missing);
/* *********************************************
*  Notice here that initial observations for each group are not defined.
*  So you MUST put the missing values as these initial values.  */

local ii,jj;

ii = 1;
do until ii > maxn;

   jj = (ii-1)*maxt+1;
   chk[jj] = missing;

   ii = ii + 1;
endo;

RETP(chk);ENDP;


/* ****************************************************** */
    proc (1) = chkmiss2(chk,chk1,chk2,chk3,temp,missing);
/* ****************************************************** */
local ii;


ii = 1;
do until ii > temp;

   if chk1[ii] eq missing; chk[ii] = missing;
   endif;
   if chk2[ii] eq missing; chk[ii] = missing;
   endif;
   if chk3[ii] eq missing; chk[ii] = missing;
   endif;

   ii = ii + 1;
endo;


RETP(CHK);ENDP;


