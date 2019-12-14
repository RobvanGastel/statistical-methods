LIBNAME SASDATA "/folders/myfolders/statistical-methods/DATA";

DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
RUN;
/* PROC PRINT DATA=WEEK1 NOOBS; */

/* Question 1.1 */
PROC FREQ data=WEEK1;
	tables TRT;
RUN;

/* Question 1.2 */
/* a) Compute the mean and variance of AGEM. */
PROC MEANS data=WEEK1 mean var std;
	var AGEM;
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
QUIT;

/* Another way of PI */
PROC REG data=WEEK1;
	model AGEM= / cli alpha=0.05;
RUN;

/* b) How many mothers were 40 years old or 
older when they became pregnant? */
/* Answer: 3 */
PROC FREQ data=WEEK1;
	tables AGEM;
	where AGEM >= 40;
RUN;

/* Is AGEM really normal distributed? */
/* It seems like it is. */
ods select HISTOGRAM, QQPLOT, PPPLOT;
PROC UNIVARIATE data=WEEK1;
   var AGEM;
   histogram value/normal;
   qqplot    value/normal;
   ppplot 	 value/normal;
RUN;

/* c) Compute a 95% confidence interval for 
the variance of AGEM.*/
ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1 cibasic(alpha=0.05);
   var AGEM;
RUN;

PROC univariate data=WEEK1 cibasic;
	var AGEM;
	histogram AGEM/ normal;
RUN;


/* Question 1.3 */
DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
RUN;

/* a) */
/* What is the interquartile range of the birth weight? */
ods select QUANTILES;
PROC UNIVARIATE data=WEEK1 cipctldf;
  var BW;
RUN;

/* p25 p50 p75 quantiles in MEANS */
PROC MEANS data=WEEK1 qrange;
	var BW;
RUN;

/* b) */
/* TODO */
DATA IQRRANGE;
	set WEEK1;
	if BW < 3365 + 870 then count= 1; 
	if BW > 3365 - 870 then count= 1; 
	else count=0;
RUN;

PROC FREQ data=IQRRANGE;
	tables count;
RUN;

/* c) */
/* λ ∈ {−2,−1/2,0,1/2,2} */
DATA WEEK1BOXCOX; 
	set WEEK1;
	AGEMMINUS2 = (-1/2)*(AGEM**-2 -1); 
	AGEMMINUS1 = (-1)*(AGEM**-1 -1); 
	AGEMMINUS12 = (-2)*(AGEM**-(0.5)-1); 
	AGEM0 = log(AGEM);
	AGEMPLUS12 = (2)*(AGEM**(1/2) -1); 
	AGEMPLUS2 = (0.5)*(AGEM**(2) -1);
RUN;

ods select histogram;
PROC UNIVARIATE data=WEEK1BOXCOX; 
   	histogram AGEM /normal;
	histogram AGEMMINUS2 /normal; 
	histogram AGEMMINUS1 /normal;
   	histogram AGEMMINUS12 /normal;
   	histogram AGEM0 /normal;
  	histogram AGEMPLUS12 /normal;
   	histogram AGEMPLUS2 /normal;
RUN;
/* Becomes normal at λ = 2 */

/* d) Prediction Interval */
/* TODO  */
PROC REG data = WEEK1;
     model BW = / cli alpha=0.05;
RUN;

/* e) */
/* Use MEANS, or show DATA and manually look up */
PROC MEANS data = WEEK1 max;
	var bw;
RUN;

/* Min and Max of cols in general */
PROC IML;
use WEEK1;
read all var _NUM_ into X[c=varNames];
close WEEK1;
 
minC = X[><, ];    /* row vector contains min of columns */
maxC = X[<>, ];    /* row vector contains max of columns */
print (minC//maxC)[r={"Min" "Max"} c=varNames];
RUN;

/* f) */
/* boys (SEX = 1) or girls (SEX = 0)  */
DATA WEEK1_GIRL;
	set WEEK1;
	where SEX = 0;
RUN;

DATA WEEK1_BOY;
	set WEEK1;
	where SEX = 1;
RUN;

PROC FREQ data=week1;
	tables sex;
RUN;

/* Girls */
PROC MEANS data=WEEK1_GIRL mean var n skew kurt;
	var BW;
RUN;

PROC MEANS data=WEEK1_BOY mean var n skew kurt;
	var BW;
RUN;

/* g) */
/* TODO */

/* Question 1.4 */
DATA WEEK1_CUSTOM; 
   input value; 
   datalines; 
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
PROC MEANS data=WEEK1_CUSTOM mean var skew kurt;
RUN;

/* Extra histogram plot */
ods select HISTOGRAM;
PROC UNIVARIATE data=WEEK1_CUSTOM;
   histogram value/normal;
   qqplot    value/normal;
   ppplot 	 value/normal;
RUN;

/* b) Confidence Interval */
ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1_CUSTOM cibasic(alpha=0.05);
   var value;
RUN;

/* c) Prediction Interval */
PROC MEANS data=WEEK1_CUSTOM mean std n; 
	var value;
    output out=custom_sumstat;
RUN;

PROC transpose data=custom_sumstat out=value_PI (DROP= _TYPE_ _FREQ_ _NAME_ _LABEL_);
	by _type_ _freq_;
	id _stat_;
RUN;

DATA value_PI;
	set value_PI;
	T = quantile("T", 1 - 0.01/2, N-1); 
	LPL = MEAN - T * std*sqrt((N+1)/ N); 
	UPL = MEAN + T * std*sqrt((N+1)/ N);
RUN;

PROC PRINT data=value_PI; 
	var LPL UPL;
RUN;

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
RUN;

DATA WEEK1;
	set WEEK1;
	AG = log(44 - GA);
RUN;

PROC MEANS data=WEEK1 mean std n; 
	var AG;
    output out=gal_sumstat;
RUN;

PROC transpose data=gal_sumstat out=gal_PI (drop= _TYPE_ _FREQ_ _NAME_ _LABEL_);
	by _type_ _freq_;
	id _stat_;
RUN;

DATA gal_PI;
	set gal_PI;
	T = quantile("T", 1 - 0.05/2, N-1); 
	LPL = MEAN - T * std*sqrt((N+1)/ N); 
	UPL = MEAN + T * std*sqrt((N+1)/ N);
RUN;

PROC PRINT data=gal_PI; 
	var LPL UPL;
RUN;
/* PI: (0.60957, 2.38299) */

/* b) Create a 95% confidence interval for the mean of the transformed 
      variable LOG(44 - GA). Can you use this confidence interval to create a 
      95% confidence interval for the mean of GA? Explain your reasoning. */
DATA WEEK1_Q6;
	set WEEK1;
	GAL = log(44 - GA);
RUN;

ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1_Q6 cibasic(alpha=0.05);
   var GAL;
RUN;
/* Can't convert Log(44 - GA) to work for GA */

/* c)  */
/* Confidence interval median */
ods select HISTOGRAM;
PROC UNIVARIATE data=WEEK1;
   histogram GA/normal;
   qqplot    GA/normal;
   ppplot 	 GA/normal;
RUN;
/* So not normally distributed */

ods select Quantiles;
PROC UNIVARIATE data=WEEK1 cipctldf cipctlnormal;
	var GA;
RUN;

/* d) */
DATA WEEK1_Q6;
	set WEEK1;
	if GA <= 38 then PRE_38 = 1;
	else PRE_38 = 0;
RUN;

/* e) */
PROC SQL;
	SELECT avg(FIS)
	FROM WEEK1_Q6
	WHERE PRE_38 = 1;
QUIT;

PROC SQL;
	SELECT avg(FIS)
	FROM WEEK1_Q6
	WHERE PRE_38 = 0;
QUIT;

PROC IML;	
	use WEEK1_Q6 where(PRETERM=1);
		read all var{FIS};
	close WEEK1_Q6;
	
	avg = mean(FIS);
	
	A = avg;
	create DATA from A[colname={'percentage'}];
		append from A;
	close DATA;
RUN;

/* f) */
PROC FREQ data=WEEK1_Q6;
	where PRE_38 = 1;
	tables FIS /binomial(wald wilson exact level=2) alpha=0.1;
RUN;

/* g) */
/* For exact intervals we have to know the distribution */


/* Question 1.7 */
%macro samples(DATAset=,ns=,n=);
	PROC surveyselect DATA=&DATAset NOPRINT method=urs n=&n out=FINAL;
RUN;

DATA FINAL;
	set FINAL; 
	sampleno = 1;
RUN;

%do sn = 2 %to &ns;
	PROC surveyselect DATA=&DATAset NOPRINT 
	method=urs n=&n out=SAMPLEI;
RUN;

DATA SAMPLEI;
	set SAMPLEI;
	sampleno = &sn;
RUN;

DATA FINAL;
	set Final SAMPLEI;
RUN;
%end;

PROC DATASETS library=work NOPRINT;
	delete SAMPLEI;
RUN;
%mend;

%samples(DATAset=WEEK1, ns=1000, n=10);

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
RUN;

PROC FREQ data=WEEK1;
 tables AGEMB /binomial(wald wilson exact level=2) alpha=0.05;
RUN;



