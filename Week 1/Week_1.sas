LIBNAME SASDATA "/folders/myfolders/statistical-methods/DATA";

DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
RUN;

/* Question 1.1 */
/* Frequence between control group */
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
	read all var{agem}; 
	close WEEK1;
	
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
	model AGEM= /cli alpha=0.05;
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

/* c) 
/* Confidence Interval for the variance of AGEM. */
ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1 cibasic(alpha=0.05);
   var AGEM;
RUN;

PROC UNIVARIATE data=WEEK1 cibasic;
	var AGEM;
	histogram AGEM/normal;
RUN;


/* Question 1.3 */
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
/* Birthweight differs atmost 1 IQR from median */
PROC SQL;
	SELECT COUNT(BW) / (SELECT count(*) FROM WEEK1)
	FROM WEEK1 
	WHERE BW BETWEEN 3365 - 870 AND 3365 + 870; 
QUIT;	/* median +- interquartile range */

/* c) */
/* Box-Cox transformation */
/* λ ∈ {−2,−1/2,0,1/2,2} */
DATA WEEK1BOXCOX; 
	set WEEK1;
	BWMINUS2 = (-1/2)*(BW**-2 -1); 
	BWMINUS1 = (-1)*(BW**-1 -1); 
	BWMINUS12 = (-2)*(BW**-(0.5)-1); 
	BW0 = log(BW);
	BWPLUS12 = (2)*(BW**(1/2) -1); 
	BWPLUS2 = (0.5)*(BW**(2) -1);
RUN;

ods select histogram;
PROC UNIVARIATE data=WEEK1BOXCOX; 
   	histogram BW /normal;
	histogram BWMINUS2 /normal; 
	histogram BWMINUS1 /normal;
   	histogram BWMINUS12 /normal;
   	histogram BW0 /normal;
  	histogram BWPLUS12 /normal;
   	histogram BWPLUS2 /normal;
RUN;
/* Becomes normal at λ = 2 */

/* d) */
/* Prediction Interval */
PROC IML;
	use WEEK1BOXCOX;
	read all var{BWPLUS2};
	close WEEK1BOXCOX;
	
	alpha=0.05;
	Ybar=mean(BWPLUS2);
	s=var(BWPLUS2);
	n=nrow(BWPLUS2);
	qT=quantile('t', alpha/2, n-1);
	
	UCL=Ybar - qT*sqrt(s/n);
	LCL=Ybar + qT*sqrt(s/n);
	UPL=Ybar - qT*sqrt((n+1) * s/n);
	LPL=Ybar + qT*sqrt((n+1) * s/n);
	
	/* Reverse the transform */
	Ybar = sqrt(2*Ybar + 1);
	LCL = sqrt(2*LCL + 1);
	UCL = sqrt(2*UCL + 1);
	LPL = sqrt(2*LPL + 1);
	UPL = sqrt(2*UPL + 1);
		
	A=Ybar||UPL||LPL||UCL||LCL;
	
	create DATA from
	A[colname={'mean', 'UPL', 'LPL', 'UCL', 'LCL'}];
	append from A;
	close DATA;
QUIT;

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

PROC FREQ data=WEEK1;
	tables sex;
RUN;

/* Compute statistics of subsets */
PROC MEANS data=WEEK1_GIRL mean var n skew kurt;
	var BW;
RUN;

PROC MEANS data=WEEK1_BOY mean var n skew kurt;
	var BW;
RUN;

/* g) */
/* The prediction interval for boys is wider as */
/* the variance is bigger for boys. */
PROC IML;
	use WEEK1BOXCOX;
	read all var{BWPLUS2} into GIRL where(SEX=:0);
	read all var{BWPLUS2} into BOY where(SEX=:1);
	close WEEK1BOXCOX;
	
	alpha = 0.05;
	Ybar_b = mean(BOY);
	Ybar = mean(BOY);
	s = var(BOY);
	n = nrow(BOY);
	qT = quantile('t', alpha/2, n-1);
	
	UPL=Ybar - qT*sqrt((n+1) * s/n);
	LPL=Ybar + qT*sqrt((n+1) * s/n);
	
	LPL_b = sqrt(2*LPL + 1);
	UPL_b = sqrt(2*UPL + 1);
	
	alpha = 0.05;
	Ybar_g = mean(GIRL);
	Ybar = mean(GIRL);
	s = var(GIRL);
	n = nrow(GIRL);
	qT = quantile('t', alpha/2, n-1);
	
	UPL=Ybar - qT*sqrt((n+1) * s/n);
	LPL=Ybar + qT*sqrt((n+1) * s/n);
	
	LPL_g = sqrt(2*LPL + 1);
	UPL_g = sqrt(2*UPL + 1);
	
	A = Ybar_b||LPL_b||UPL_b||Ybar_g||LPL_g||UPL_g;
	create DATA from A[colname={'mean Boy' 'LPL Boy' 'UPL Boy' 'mean Girl' 'LPL Girl' 'UPL Girl'}];
		append from A;
	close DATA;
QUIT;

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

/* Plots */
ods select HISTOGRAM;
PROC UNIVARIATE data=WEEK1_CUSTOM;
   histogram value/normal;
   qqplot    value/normal;
   ppplot 	 value/normal;
RUN;

/* b) 
/* Confidence Interval */
ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1_CUSTOM cibasic(alpha=0.05);
   var value;
RUN;

/* c) 
/* Prediction Interval */
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
/* 95% CI for log(μ) is equal to (−0.137, 1.128) */
/* a) 
/* Use the CI of log(μ) to derive a 95% CI for μ.  */
/* We can apply e^log(μ) = μ so, (0.871970226, 3.08947137) */

/* b) 
/* Test H0: μ=3 vs. H1: μ !=3 at a significance level  */
/* of α=5%. */
/* The CI contains 3 so we can't reject H0. */

/* Question 1.6 */
/* a) 
/* Use the transformation LOG(44 - GA) to compute a 95%  */
/* PI for a single new observation of GA */
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

DATA gal_PI;
	set gal_PI;
	LPL = 44 - exp(LPL);
	UPL = 44 - exp(UPL);
RUN;

PROC PRINT data=gal_PI;
	var UPL LPL;
RUN;
/* PI: (33.1628, 42.1604) */

/* b)  */
/* Can't convert Log(44 - GA) to work for GA */
DATA WEEK1_Q6;
	set WEEK1;
	GAL = log(44 - GA);
RUN;

ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1_Q6 cibasic(alpha=0.05);
   var GAL;
RUN;

/* c)  */
/* Confidence interval median */
ods select histogram;
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
	use WEEK1_Q6 where(PRE_38=1);
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
/* Its done the same as F but use Clopper-Pearson (exact) */

/* Question 1.7 */
%macro samples(DATAset=,ns=,n=);
	PROC SURVEYSELECT data=&dataset noprint method=urs n=&n out=FINAL;
RUN;

DATA FINAL;
	set FINAL; 
	sampleno = 1;
RUN;

%do sn = 2 %to &ns;
	PROC SURVEYSELECT data=&dataset noprint
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

PROC DATASETS library=work noprint;
	delete SAMPLEI;
RUN;
%mend;

/* a) */
%samples(dataset=WEEK1, ns=1000, n=10);

PROC PRINT DATA=FINAL;
RUN;

/* b) */
PROC MEANS data=FINAL nway;
	class sampleno;
	var AGEM;
RUN;

/* Or by SQL */
PROC SQL;
	CREATE TABLE AVG_FINAL AS
		SELECT avg(AGEM) as mean, var(AGEM) as variance
		FROM FINAL
		GROUP BY SAMPLENO;
RUN;

/* c) */
/* We can see a normal distribution for the mean */
/* and a Chi-Square distribution for the variance */
ods select hisotgram;
PROC UNIVARIATE data=AVG_FINAL;
   histogram mean/normal;
   histogram variance/normal;
RUN;

/* Question 1.8 */
/* a) */
PROC TTEST h0=3200 sides=2 data=WEEK1;
	var BW;
RUN;
/* We reject H0 with p-value = 0.285 and */
/* test statistic = 2.20 */

/* b) */
/* the p-value = 0.0143 */
PROC TTEST h0=3200 sides=U data=WEEK1;
	var BW;
RUN;

/* c) */
/* Does the CLT apply? */
/* Seems like it does */
%samples(dataset=WEEK1, ns=100, n=253);

PROC MEANS data=FINAL mean noprint; 
	var BW;
	by sampleno;
	output out=MEANSBW mean=BW_MEAN; 
RUN;

PROC UNIVARIATE data=MEANSBW; 
	hist BW_MEAN /normal;
run; 
%mend;

/* d) */
/* 3 identical copies will not change the effect size */

/* e) */
/* It equals the p-value of question a */

/* f) */
ods graphics off;
ods exclude all;
PROC TTEST data=FINAL h0=3200 sides=u alpha=0.05;
	var BW;
	by SAMPLENO;
	ods output ttests=POWER(keep=PROBT);
RUN;
ods exclude none; 
ods graphics on;

/* g) */
/* Power is about 66%, Reject H0 */
DATA RESULTS; 
   set POWER; 
   RejectH0 = (Probt <= 0.05); 
RUN; 
 
/* Compute proportion (# reject H0)/ # samples */
/* and CI */
PROC FREQ data=RESULTS; 
   tables RejectH0 / nocum binomial(level='1');
RUN;

/* h) */
/* Power is expected to go down with smaller  */
/* sample size */


