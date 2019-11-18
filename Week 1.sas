LIBNAME SASDATA "/folders/myfolders/applied-stats/data";

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
/* a) */
PROC MEANS DATA=WEEK1 mean var std;
	VAR AGEM;
RUN;

/* confidence intervals */
ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1 cibasic(alpha=0.05);
   var AGEM;
run;
/* TODO: 95% prediction interval for a single new observation of AGEM. */

/* b) How many mothers were 40 years old or older when they became pregnant? */
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

/* Question 1.3 */
/* a) */
DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
run;

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
proc means data=WEEK1 mean std n; 
	var BW;
    output out=agem_sumstat;
run;

proc transpose data=agem_sumstat out=bw_PI (DROP= _TYPE_ _FREQ_ _NAME_ _LABEL_);
	by _type_ _freq_;
	id _stat_;
run;

data bw_PI;
	set bw_PI;
	T = QUANTILE("T", 1 - 0.05/2, N-1); 
	LPL = MEAN - T * std*sqrt((N+1)/ N); 
	UPL = MEAN + T * std*sqrt((N+1)/ N);
run;

proc print DATA=bw_PI; 
	var LPL UPL;
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

/* proc sort data=WEEK1 out=WEEK1_SORT; */
/* 	by descending BW; */
/* run; */

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
/* a) Use the confidence interval of log(μ) to derive a 95% confidence interval for μ. */


/* b) Test H0: μ=3 vs. H1: μ !=3 at a significance level of α=5%. */


/* Example slides */
DATA WEEK1;
	set SASDATA.IVF;
	where PER=4; AGEMB = (AGEM<30); 
	keep AGEMB AGEM;
run;

PROC FREQ data=WEEK1;
 tables AGEMB /binomial(wald wilson exact level=2) alpha=0.05;
run;



