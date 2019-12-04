LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
run;
/* PROC PRINT data=WEEK1 NOOBS; */

/* Question 1.1 */
PROC FREQ data=WEEK1;
	tables TRT;
Run;

/* Question 1.2 */
/* a) Compute the mean and variance of AGEM. */
PROC MEANS DATA=WEEK1 mean var std;
	VAR AGEM;
RUN;

/* prediction interval */
PROC IML;
use WEEK1;
read all var{agem}; close WEEK1;
alpha=0.05;
Ybar=mean(agem); 
s=var(agem);
n=nrow(agem); 
qT=quantile('t',alpha/2,n-1); 
UPL=Ybar-qT*sqrt((n+1)*s/n); 
LPL=Ybar+qT*sqrt((n+1)*s/n); 
A=Ybar||LPL||UPL;

create DATA from A[colname={'mean' 'LPL' 'UPL'}]; 
append from A;
close DATA;
quit;

/* Another way of PI */
proc reg data=WEEK1;
	model AGEM= / cli alpha=0.05;
run;

/* b) How many mothers were 40 years old or 
older when they became pregnant? */
/* Answer: 3 */
PROC FREQ data=WEEK1;
	tables AGEM;
	where AGEM >= 40;
Run;

/* Is AGEM really normal distributed? */
/* It seems like it is. */
ods select HISTOGRAM, QQPLOT, PPPLOT;
PROC UNIVARIATE data=WEEK1;
   VAR AGEM;
   HISTOGRAM AGEM/NORMAL;
   QQPLOT    AGEM/NORMAL;
   PPPLOT 	 AGEM/NORMAL;
RUN;

/* c) Compute a 95% confidence interval for 
the variance of AGEM.*/
ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1 cibasic(alpha=0.05);
   var AGEM;
run;

proc univariate data=WEEK1 cibasic;
	var AGEM;
	histogram AGEM/ normal;
run;


/* Question 1.3 */
DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
run;

/* a) */
/* What is the interquartile range of the birth weight? */
ods select QUANTILES;
PROC UNIVARIATE data=WEEK1 cipctldf;
  var BW;
run;

/* p25 p50 p75 quantiles in MEANS */
PROC MEANS data=WEEK1 qrange;
	var BW;
run;

/* b) */
/* TODO */
DATA IQRRANGE;
	SET WEEK1;
	IF BW < 3365 + 870 THEN count= 1; 
	IF BW > 3365 - 870 THEN count= 1; 
	ELSE count=0;
run;

PROC FREQ data=IQRRANGE;
	tables count;
Run;

/* c) */
/* λ ∈ {−2,−1/2,0,1/2,2} */
DATA WEEK1BOXCOX; 
	SET WEEK1;
	AGEMMINUS2 = (-1/2)*(AGEM**-2 -1); 
	AGEMMINUS1 = (-1)*(AGEM**-1 -1); 
	AGEMMINUS12 = (-2)*(AGEM**-(0.5)-1); 
	AGEM0 = log(AGEM);
	AGEMPLUS12 = (2)*(AGEM**(1/2) -1); 
	AGEMPLUS2 = (0.5)*(AGEM**(2) -1);
run;

ods select histogram;
PROC UNIVARIATE data=WEEK1BOXCOX; 
   	histogram AGEM /normal;
	histogram AGEMMINUS2 /normal; 
	histogram AGEMMINUS1 /normal;
   	histogram AGEMMINUS12 /normal;
   	histogram AGEM0 /normal;
  	histogram AGEMPLUS12 /normal;
   	histogram AGEMPLUS2 /normal;
run;
/* Becomes normal at λ = 2 */

/* d) Prediction Interval */
/* TODO  */
proc reg data=WEEK1;
     model BW= / cli alpha=0.05;
run;

/* e) */
/* Use means, or show data and manually look up */
PROC MEANS data=WEEK1 max;
	var bw;
run;

/* Min and Max of cols in general */
proc iml;
use WEEK1;
read all var _NUM_ into X[c=varNames];
close WEEK1;
 
minC = X[><, ];    /* row vector contains min of columns */
maxC = X[<>, ];    /* row vector contains max of columns */
print (minC//maxC)[r={"Min" "Max"} c=varNames];
run;

/* f) */
/* boys (SEX = 1) or girls (SEX = 0)  */
data week1_girl;
	set WEEK1;
	where SEX = 0;
run;

data week1_boy;
	set WEEK1;
	where SEX = 1;
run;

proc freq data=week1;
	tables sex;
run;

/* Girls */
proc means data=week1_girl mean var n skew kurt;
	var BW;
run;

proc means data=week1_boy mean var n skew kurt;
	var BW;
run;

/* g) */
/* TODO */

/* Question 1.4 */
DATA WEEK1_CUSTOM; 
   INPUT value; 
   DATALINES; 
	25.0 
	27.4 
	17.1 
	22.1 
	20.8 
	21.3 
	22.5 
	29.2 
	27.9 
	25.7 
	24.7 
	18.8
 ;
RUN; 

/* a) */
proc means data=WEEK1_CUSTOM mean var skew kurt;
run;

/* Extra histogram plot */
ods select HISTOGRAM;
PROC UNIVARIATE data=WEEK1_CUSTOM;
   HISTOGRAM value/NORMAL;
   QQPLOT    value/NORMAL;
   PPPLOT 	 value/NORMAL;
RUN;

/* b) Confidence Interval */
ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1_CUSTOM cibasic(alpha=0.05);
   var value;
run;

/* c) Prediction Interval */
proc means data=WEEK1_CUSTOM mean std n; 
	var value;
    output out=custom_sumstat;
run;

proc transpose data=custom_sumstat out=value_PI (DROP= _TYPE_ _FREQ_ _NAME_ _LABEL_);
	by _type_ _freq_;
	id _stat_;
run;

data value_PI;
	set value_PI;
	T = QUANTILE("T", 1 - 0.01/2, N-1); 
	LPL = MEAN - T * std*sqrt((N+1)/ N); 
	UPL = MEAN + T * std*sqrt((N+1)/ N);
run;

proc print DATA=value_PI; 
	var LPL UPL;
run;

/* Question 1.5 */
/* 95% confidence interval for log(μ) equal to (−0.137, 1.128) */
/* a) Use the confidence interval of log(μ) to derive a 95% 
      confidence interval for μ. */
/* e^log(μ) = μ */
/* (0.871970226, 3.08947137) */

/* b) Test H0: μ=3 vs. H1: μ !=3 at a significance level of α=5%. */
/* The confidence interval contains 3 so we can't reject H0. */

/* Question 1.6 */
/* a) Use the transformation LOG(44 - GA) to compute a 95% prediction 
      interval for a single new observation of the gestational age. */
DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
run;

DATA WEEK1;
	set WEEK1;
	AG = log(44 - GA);
run;

proc means data=WEEK1 mean std n; 
	var AG;
    output out=gal_sumstat;
run;

proc transpose data=gal_sumstat out=gal_PI (DROP= _TYPE_ _FREQ_ _NAME_ _LABEL_);
	by _type_ _freq_;
	id _stat_;
run;

data gal_PI;
	set gal_PI;
	T = QUANTILE("T", 1 - 0.05/2, N-1); 
	LPL = MEAN - T * std*sqrt((N+1)/ N); 
	UPL = MEAN + T * std*sqrt((N+1)/ N);
run;

proc print DATA=gal_PI; 
	var LPL UPL;
run;
/* PI: (0.60957, 2.38299) */

/* b) Create a 95% confidence interval for the mean of the transformed 
      variable LOG(44 - GA). Can you use this confidence interval to create a 
      95% confidence interval for the mean of GA? Explain your reasoning. */
DATA WEEK1_Q6;
	set WEEK1;
	GAL = log(44 - GA);
run;

ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1_Q6 cibasic(alpha=0.05);
   var GAL;
run;
/* Can't convert Log(44 - GA) to work for GA */

/* c)  */
ods select HISTOGRAM;
PROC UNIVARIATE data=WEEK1;
   HISTOGRAM GA/NORMAL;
   QQPLOT    GA/NORMAL;
   PPPLOT 	 GA/NORMAL;
RUN;

/* d-g) */
/* TODO */



/* Question 1.7 */
%macro samples(dataset=,ns=,n=);
	proc surveyselect data=&dataset NOPRINT method=urs n=&n out=FINAL;
run;

data FINAL;
	set FINAL; 
	sampleno = 1;
run;

%do sn = 2 %to &ns;
	proc surveyselect data=&dataset NOPRINT 
	method=urs n=&n out=SAMPLEI;
run;

data SAMPLEI;
	set SAMPLEI;
	sampleno = &sn;
run;

data FINAL;
	set Final SAMPLEI;
run;
%end;

proc dataset library=work NOPRINT;
	delete SAMPLEI;
run;
%mend;

%samples(dataset=WEEK1, ns=1000, n=10);

/* a-c) */
/* TODO */

/* Question 1.8 */
/* a-h) */
/* TODO */

/* Example slides */
DATA WEEK1;
	set SASDATA.IVF;
	where PER=4; AGEMB = (AGEM<30); 
	keep AGEMB AGEM;
run;

PROC FREQ data=WEEK1;
 tables AGEMB /binomial(wald wilson exact level=2) alpha=0.05;
run;



