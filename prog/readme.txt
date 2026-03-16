Pooled Mean Group Estimation of Dynamic Heterogeneous Panels
------------------------------------------------------------

These files allow you to replicate the estimates in the empirical section of
Pesaran, Shin and Smith, Pooled Mean Group Estimation of Dynamic Hetero-
geneous Panels, JASA (forthcoming).

There are two examples:
1) the consumption function in 24 OECD countries;
2) energy demand in 10 Asian developing countries.

Section 1 of the GAUSS program files can also be edited to allow you to
carry out Pooled Mean Group estimation on your own data sets.

The files:

GAUSS program files:

   JASA1.PRG for the consumption function example
   JASADAT1.PRG prepares and transforms raw data for example 1
   JASA2.PRG for the energy demand in Asia example

Subroutines which must be placed in the same directory as the relevant
program file(s):

   SUB0.PRC: common subprocedures
   SUB1.PRC: subprocedures necessary for the PMG estimation when all the
             long-run parameters are restricted to be the same across groups
   SUB2.PRC: subprocedures necessary when a subset of the long-run parameters
             is restricted to be the same across groups
   FIX.PRC:  subprocedures necessary for fixed effects estimation
   PR.PRC:   subprocedures for printing the full version results
   PRSM.PRC: subprocedures for printing the results for a small version

Data files:

Example 1:

Raw data files

        CP.DAT  Nominal private consumption
        RCP.DAT Real private consumption
        NDI.DAT National disposable income
        PDI.DAT Private non-property income
        POP.DAT Population

Transformed data files

        DP.DAT  Inflation
        LNDI.DAT Log real per-capita national disposable income
        LPC.DAT  Log real per-capita private consumption
        LPDI.DAT Log real percapita private non-property income

Example 2:

        EDA.DAT  log energy demand and log prices

Group names for each example:

        J1NAME.DAT
        J2NAME.DAT

Output files:

        JASA1A.OUT
        JASA1B.OUT
        JASA2.OUT

Summary of instructions:

The programs are annotated with instructions. Below follows a brief
summary on how to edit the program.

When editing Section 1 of the program one needs to specify:

a name for the output file;
the maximum number of iterations for the pooled ML estimation;
a value for missing observations;
the maximum number of time periods MAXT;
the maximum number of groups MAXN;
the data to load;
the list of group names to load;
the dependent variable - Y;
the distributed lag regressors - X;
the deterministic regressors - Z.

The data must be in the form of ASCII files containing no alphabetics
with a column for each variable with temp=MAXT*MAXN elements arranged
with all the time series observations on the first group, followed by
all the time series observations on the second group, etc. with missing
values explicitly included.

If the files contain single variables called var1 and var2, the command
takes the form:

load var1[temp,1]=var1.dat;
load var2[temp,1]=var2.dat;

If the file contains two variables the command takes the form:

load var[temp,2]=var.dat;
var1=var[.,1];
var2=var[.,2];

In the program this is referred to as Data Option 1. Data Option 2 is
not available in the program at the moment.

The program will provide deterministic regressors for intercept, trend,
and monthly or quarterly seasonals. The intercept must be specified as
the last element in Z.

The program will also allow you to use data relative to some reference
variable. Such a series must be specified as the last column of X.

You can also provide a file with a list of group names which will be
used in the output.

Having edited the program to provide the initial information and loaded
the data it can be run in GAUSS. You will be presented with various
choices.

The program will ask you whether:

- you have loaded names or numeric codes for the groups
- you want to use all the groups or a subset
- the regressors include a relative variable
- you want to run the full version or the small version of the program
- you want to express the data as deviations from the cross-section means
  for each period, to remove common time effects?
- you want to set the lag orders yourself or have them chosen by a
  model selection criterion
- you want all the long-run parameters to be restricted to be the same
  across groups or only a subset of them
- you want initial values for the long-run parameters to be set by
  yourself or set equal to
        (1) the mean of the OLS estimates,
        (2) the static fixed effect estimates, or
        (3) the dynamic fixed effect estimates?

  (*Option (3) cannot be used when lag lengths are chosen by a model
  selection criterion)

  (*Option (1) cannot be used when using a relative variable and choosing
  the option that a subset of the long-run parameters is restricted to be
  the same across groups)

- you want to use the Back Substitution Algorithm (BSA) or the Newton-
  Raphson Algorithm (NRA)? If you chose the small version of the program,
  you must use BSA and only the summary estimation results without
  standard errors will be provided.

The full version of the program then gives you an output file which contains

- Summary on the information and the options used in the PMG estimation
- Summary on the alternative estimates of the long-run paramters
- Pooled Mean Group and OLS estimates for each group
- Pooled Mean Group and Mean Group estimates
- Static and Dynamic Fixed Effects estimates
- Various diagnostic statistics
- OLS-based group-specific summary results for the long-run coefficients
- PMG-based group-specific summary results for the long-run coefficients

The small version of the program gives an output file which contains

- Summary on the information and the options used in the PMG estimation
- Summary on the alternative estimates of the long-run paramters
- Mean Group estimates
- Static and Dynamic Fixed Effects estimates
- OLS-based group-specific summary results for the long-run coefficients

f
