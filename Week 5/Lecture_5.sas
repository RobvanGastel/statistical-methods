LIBNAME SAS "D:\Documents\PhD\Teaching\Applied statistics\2018-2019\SAS code";

data CU;
	input CU@@;
	datalines;
102 104 101 97 95 106 103 98 96 99
;
run;

PROC UNIVARIATE DATA=CU;
	VAR CU;
	QQPLOT CU/NORMAL(MU=est SIGMA=est);
RUN;


DATA RCT;
 SET SAS.Cellsaver;
RUN;

DATA RCT1;
SET RCT;
where TIME=3;
run;

DATA RCT2;
SET RCT;
where TIME=3 AND ID<100;
run;

proc sort data=rct2;
by trt;
run;

ods select histogram;
PROC UNIVARIATE DATA=RCT2;
	BY TRT;
	VAR RESP;
	HISTOGRAM RESP /NORMAL;
RUN;

ods select probplot;
PROC UNIVARIATE DATA=RCT2;
	BY TRT;
	VAR RESP;
	PROBPLOT RESP /NORMAL(MU=est SIGMA=est);
RUN;

PROC MIXED DATA=RCT2 METHOD=TYPE3 CL;
	CLASS TRT;
	MODEL RESP = TRT /SOLUTION OUTP=PRED CL;
RUN;

ods select probplot;
PROC UNIVARIATE DATA=PRED;
	VAR RESID;
	PROBPLOT RESID /NORMAL(MU=est SIGMA=est);
RUN;

PROC UNIVARIATE DATA=PRED;
	VAR RESID;
	HISTOGRAM RESID /NORMAL;
RUN;

PROC UNIVARIATE DATA=PRED NORMAL;
	VAR RESID;
RUN;

proc iml;


DATA approx;
*N=712;
*G1=0.43510643;
*G2=0.51594979;
N=96;
G1=0.88861759;
G2=2.63603486;
b1=(N-2)*G1/(sqrt(N*(N-1)));
b2=G2*((N-2)*(N-3))/((N+1)*(N-1))+3*(N-1)/(N+1);
*JB=N*(b1**2/6+(G2-3)**2/24);
Cn=(3*(N**2+27*N-70)*(N+1)*(N+3))/((N-2)*(N+5)*(N+7)*(N+9));
Wn2=-1+SQRT(2*(Cn-1));
Alphan=SQRT(2/(Wn2-1));
Dn=1/sqrt(log(sqrt(Wn2)));
Bn=sqrt((N+1)*(N+3)/(6*(N-2)))*b1;
Ts=Dn*log(Bn/Alphan+sqrt(1+(Bn/Alphan)**2));
Mun=3*(N-1)/(N+1);
Sigman=sqrt((24*N*(N-2)*(N-3))/((N+3)*(N+5)*(N+1)**2));
Gamma1n=((6*(N**2-5*N+2))/((N+7)*(N+9)))*sqrt(6*(N+3)*(N+5)/(N*(N-2)*(N-3)));
An=6+(8/(Gamma1n))*(2/Gamma1n+sqrt(1+4/(Gamma1n**2)));
Un=(b2-Mun)/Sigman;
Tk=sqrt(9*An/2)*((9*An-2)/(9*An)-((1-2/An)/(1+Un*sqrt(2/(An-4))))**(1/3));
K2=Tk**2+Ts**2;
Ps=2*min(cdf('Normal',Ts,0,1),1-cdf('Normal',Ts,0,1));
Pk=2*min(cdf('Normal',Tk,0,1),1-cdf('Normal',Tk,0,1));
PK2=1-cdf('chisq',K2,2);
*PJB=1-cdf('chisq',JB,2);
run;

proc print data=approx;
var Ts Tk K2 Ps Pk PK2 ;
run;

PROC GLM DATA=PRED;
	CLASS TRT;
	MODEL RESID = TRT;
	MEANS TRT/ HOVTEST=Bartlett;
RUN;

*Check approx;
data A;
do i = 1 to 96;
   X = rand("Normal");     
   output;
end;
run;

proc univariate data=A normal;
var X;
run;

DATA approx;
*N=712;
*G1=0.43510643;
*G2=0.51594979;
N=96;
G1=-0.0910325;
G2=-0.0075493;
b1=(N-2)*G1/(sqrt(N*(N-1)));
b2=G2*((N-2)*(N-3))/((N+1)*(N-1))+3*(N-1)/(N+1);
*JB=N*(b1**2/6+(G2-3)**2/24);
Cn=(3*(N**2+27*N-70)*(N+1)*(N+3))/((N-2)*(N+5)*(N+7)*(N+9));
Wn2=-1+SQRT(2*(Cn-1));
Alphan=SQRT(2/(Wn2-1));
Dn=1/sqrt(log(sqrt(Wn2)));
Bn=sqrt((N+1)*(N+3)/(6*(N-2)))*b1;
Ts=Dn*log(Bn/Alphan+sqrt(1+(Bn/Alphan)**2));
Mun=3*(N-1)/(N+1);
Sigman=sqrt((24*N*(N-2)*(N-3))/((N+3)*(N+5)*(N+1)**2));
Gamma1n=((6*(N**2-5*N+2))/((N+7)*(N+9)))*sqrt(6*(N+3)*(N+5)/(N*(N-2)*(N-3)));
An=6+(8/(Gamma1n))*(2/Gamma1n+sqrt(1+4/(Gamma1n**2)));
Un=(b2-Mun)/Sigman;
Tk=sqrt(9*An/2)*(1-2/(9*An)-((1-2/An)/(1+Un*sqrt(2/(An-4))))**(1/3));
K2=Tk**2+Ts**2;
Ps=2*min(cdf('Normal',Ts,0,1),1-cdf('Normal',Ts,0,1));
Pk=2*min(cdf('Normal',Tk,0,1),1-cdf('Normal',Tk,0,1));
PK2=1-cdf('chisq',K2,2);
*PJB=1-cdf('chisq',JB,2);
run;

proc print data=approx;
var Ts Tk K2 Ps Pk PK2;
run;

*Week 5-2;

data L10;
	input OBS Y @@;
	datalines;
	1 	104.3 	6 	99.0 	11 	112.5 	16 	114.0 
	2 	132.4 	7 	109.4 	12 	98.8 	17 	98.9 
	3 	112.4 	8 	101.9 	13 	97.0 	18 	112.1 
	4 	100.7 	9 	100.5 	14 	114.8 	19 	100.6 
	5 	105.3 	10 	110.5 	15 	110.7 	20 	119.3 
	;
run;

*Doornbos test;
proc means data=L10 mean std n;
	var Y;
	output out=L10_sumstat mean=mean median=median std=std n=n;
run;

data L10;
set L10;
if _n_=1 then set L10_sumstat;
drop _TYPE_ _FREQ_;
run;


data DOORNBOS;
	SET L10;
	U = (Y-MEAN)/STD;
	W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
	DOORNBOS_Y=ABS(W);
	CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
	P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),1-CDF('T', DOORNBOS_Y, N-2))*N,1);
RUN;

PROC SORT DATA =Doornbos;
	BY descending DOORNBOS_Y;
RUN;

proc print data=Doornbos;
run;

*Grubbs test;
DATA Grubbs;
	SET L10;
	U = (Y-MEAN)/STD;
	Grubbs=ABS(U);
	*Adapt C_exact manually for different N;
	C_onesided_exact= 2.56;
	C_twosided_exact= 2.71;	
	t = quantile("t", 0.05 / (2*N), N-2);
	C_twosided_approx = (n-1) * sqrt(t**2 / (n * (t**2 + n - 2)));
	u_inv = u*sqrt((n-2)) / sqrt(n-1-u**2);
	p_twosided_approx = min(2*n*MIN(1-cdf("t", u_inv, n-2),cdf("t", u_inv, n-2)), 1);
RUN;

PROC SORT DATA =Grubbs;
	BY descending Grubbs;
RUN;


proc print data=Grubbs;
run;

*Hampel's rule;
DATA Hampel;
	SET L10;
	D = ABS(Y-MEDIAN);
RUN;

proc means data=Hampel;
	var D;
	output out=Hampel_med median=medianD;
run;

data Hampel;
set Hampel;
if _n_=1 then set Hampel_med;
drop _TYPE_ _FREQ_;
run;

DATA Hampel;
	SET Hampel;
	Z = ABS(Y-MEDIAN)/MEDIAND;
	H = (Z>3.5);
RUN;


Proc sort data=Hampel;
by descending H;
run;

PROC PRINT DATA=Hampel;
RUN;

*Tukey's method;
PROC BOXPLOT DATA=L10;
  PLOT Y*MEAN/BOXSTYLE=SCHEMATIC;
RUN;

DATA TUKEY;
SET L10;
IQR=p75-p25;
LOWERT = p25 - 1.5*IQR;
UPPERT = p75 + 1.5*IQR;
RUN;

DATA TUKEY;
SET TUKEY;
T=(Y>UPPERT OR Y<LOWERT);
RUN;

PROC SORT data=TUKEY;
by descending T;
run;

PROC PRINT DATA=TUKEY;
RUN;

*Normality test with and without outlier;
proc univariate data=L10 NORMAL;
var Y;
run;

data L10_NEW;
SET L10;
WHERE OBS NE 2;
run;

proc univariate data=L10_NEW NORMAL;
var Y;
run;
