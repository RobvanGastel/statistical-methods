LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Summary data on the heart disease mortality rate */
/* (per 100000 people) per ethnicity group in 15 different */
/* large cities between 2014 and 2016. */
DATA BCHC;
	SET SASDATA.BCHC;
	ID = _n_; 
RUN;

/* Place: The city the data originates from. */
/* Year: The calendar year in which the data has been collected. */
/* Race_Ethnicity: Ethical background; Black, White, Hispanic,  */
/* 				   Asian/PI. */
/* Mortality: Number of people that died as a consequence  */
/*            of a heart disease per 100000. */
PROC PRINT DATA=BCHC;
RUN;

/* Both Asian/PI and Hispanic have 8 missing values for mortality */
PROC MEANS data=BCHC n mean std var median min max nmiss;
	class Race_Ethnicity;
RUN;

/* The data-set is balanced with all 4 groups having nearly 40 observatio */
PROC FREQ DATA=BCHC;
	TABLES Race_Ethnicity;
RUN;

/* Check for normality assumptions per group */
PROC SORT DATA=BCHC;
	by Race_Ethnicity;
RUN;
ods select histogram;

PROC UNIVARIATE data=BCHC normal;
	by Race_Ethnicity;
	var mortality;
	histogram mortality/normal;
RUN;

PROC UNIVARIATE data=BCHC normal;
	histogram mortality/normal;
RUN;


/* Therfore, we use the boxcox transform to choose */
/* a lambda which makes it closest to a normal distribution */
DATA BCHC_transformed;
	SET BCHC;
	Mortality13 = (3)*(Mortality**(1/3) -1);
RUN;

proc univariate data=BCHC_transformed normaltest;
	var Mortality13 Mortality;
	histogram Mortality13 /normal;
	histogram Mortality /normal;
RUN;


/* Hampel's Rule ' */
%Hampel(dataset=BCHC_transformed, var=Mortality13, id=ID);
PROC SQL;
	SELECT * FROM BCHC_transformed
	WHERE ID IN (138, 142, 144, 140);
RUN;
/* Extreme observations San Antonio */

%Hampel(dataset=BCHC_transformed, var=Mortality13, id=ID);


/* MCAR */
/* Assumptions of MCAR is the data is normally distributed, */
/* so we use the transformation for lambda 1/3. */

/* In fact, simulation studies suggest that mean imputation */
/* is possibly the worst missing data handling method  */
/* available (because the variability decreases). */
/* Consequently, in no situation is mean  */
/* imputation defensible, and you should absolutely avoid  */
/* this approach. (Enders, C. K. (2010). Applied missing  */
/* data analysis. New York: Guilford.) */

/* Create imputations */
PROC MI data=BCHC_transformed out=BCHC_transformed_mv nimpute=1 seed=1 minimum=0;
	class Race_Ethnicity Place Year;
	fCS logistic(Year);
	FCS discrim(Race_Ethnicity Place);
	var Mortality13 Year Race_Ethnicity Place;
RUN;


/* Dataset with imputation values */
DATA BCHC_IM;
	SET BCHC_transformed_mv;
	final_mortality = ((Mortality13 / 3) + 1) ** 3;
RUN;


/* Question 1: Is the impact of race ethnicity on mortality */
/* rate significant? If so, can this impact be quantified? */
/* Add columns with numeric categories for race and place */
DATA BCHC_IM;
	set BCHC_IM;
	if Race_Ethnicity = "Asian/PI" then R = 0;
	else if Race_Ethnicity = "Black" then R = 1;
	else if Race_Ethnicity = "Hispanic" then R = 2;
	else if Race_Ethnicity = "White" then R = 3;
	else put 'ERROR: Race not found ' Race_Ethnicity;
RUN;


PROC GLM data=BCHC_IM;
	class Race_Ethnicity;
	model final_mortality = Race_Ethnicity;
	means Race_Ethnicity / hovtest=levene hovtest=bf;
RUN;

/* Wilcoxon rank-sum */
ods output WilcoxonScores=WMW_01 (keep= Class N SumOfScores);
PROC NPAR1WAY data=BCHC_IM correct=NO;
	where R = 0 or R = 1;
	class R;
	var final_mortality;
	title 'Asian/PI & Black';
	exact wilcoxon /mc;
	ods select WilcoxonTest;
RUN;
PROC IML;
	use WMW_01;
		read all var{Class N SumOfScores};
	close WMW_01;
	
	G = Num(Class);
	U = SumOfScores - N#(N+1)/2;
	P = U / prod(N);
	
	A = G||N||U||P;
	create MWU from A [colname={'Group' 'N' 'U' 'P'}];
		append from A;
	close MWU;
	
	PROC PRINT data=MWU;
		title 'Asian/PI & Black';
	RUN;
QUIT;

ods output WilcoxonScores=WMW_02 (keep= Class N SumOfScores);
PROC NPAR1WAY data=BCHC_IM correct=NO;
	where R = 0 or R = 2;
	class R;
	var final_mortality;
	title 'Asian/PI & Hispanic';
	exact wilcoxon /mc;
	ods select WilcoxonTest;
RUN;
PROC IML;
	use WMW_02;
		read all var{Class N SumOfScores};
	close WMW_02;
	
	G = Num(Class);
	U = SumOfScores - N#(N+1)/2;
	P = U / prod(N);
	
	A = G||N||U||P;
	create MWU from A [colname={'Group' 'N' 'U' 'P'}];
		append from A;
	close MWU;
	
	PROC PRINT data=MWU;
		title 'Asian/PI & Hispanic';
	RUN;
QUIT;

ods output WilcoxonScores=WMW_03 (keep= Class N SumOfScores);
PROC NPAR1WAY data=BCHC_IM correct=NO;
	where R = 0 or R = 3;
	class R;
	var final_mortality;
	title 'Asian/PI & White';
	exact wilcoxon /mc;
	ods select WilcoxonTest;
RUN;
PROC IML;
	use WMW_03;
		read all var{Class N SumOfScores};
	close WMW_03;
	
	G = Num(Class);
	U = SumOfScores - N#(N+1)/2;
	P = U / prod(N);
	
	A = G||N||U||P;
	create MWU from A [colname={'Group' 'N' 'U' 'P'}];
		append from A;
	close MWU;
	
	PROC PRINT data=MWU;
		title 'Asian/PI & White';
	RUN;
QUIT;

ods output WilcoxonScores=WMW_12 (keep= Class N SumOfScores);
PROC NPAR1WAY data=BCHC_IM correct=NO;
	where R = 1 or R = 2;
	class R;
	var final_mortality;
	title 'Black & Hispanic';
	exact wilcoxon /mc;
	ods select WilcoxonTest;
RUN;
PROC IML;
	use WMW_12;
		read all var{Class N SumOfScores};
	close WMW_12;
	
	G = Num(Class);
	U = SumOfScores - N#(N+1)/2;
	P = U / prod(N);
	
	A = G||N||U||P;
	create MWU from A [colname={'Group' 'N' 'U' 'P'}];
		append from A;
	close MWU;
	
	PROC PRINT data=MWU;
		title 'Black & Hispanic';
	RUN;
QUIT;

ods output WilcoxonScores=WMW_13 (keep= Class N SumOfScores);
PROC NPAR1WAY data=BCHC_IM correct=NO;
	where R = 1 or R = 3;
	class R;
	var final_mortality;
	title 'Black & White';
	exact wilcoxon /mc;
	ods select WilcoxonTest;
RUN;
PROC IML;
	use WMW_13;
		read all var{Class N SumOfScores};
	close WMW_13;
	
	G = Num(Class);
	U = SumOfScores - N#(N+1)/2;
	P = U / prod(N);
	
	A = G||N||U||P;
	create MWU from A [colname={'Group' 'N' 'U' 'P'}];
		append from A;
	close MWU;
	
	PROC PRINT data=MWU;
		title 'Black & White';
	RUN;
QUIT;

ods output WilcoxonScores=WMW_23 (keep= Class N SumOfScores);
PROC NPAR1WAY data=BCHC_IM correct=NO;
	where R = 2 or R = 3;
	class R;
	var final_mortality;
	title 'Hispanic & White';
	exact wilcoxon /mc;
	ods select WilcoxonTest;
RUN;
PROC IML;
	use WMW_23;
		read all var{Class N SumOfScores};
	close WMW_23;
	
	G = NUM(Class);
	U = SumOfScores - N#(N+1)/2;
	P = U / prod(N);
	
	A = G||N||U||P;
	create MWU from A [colname={'Group' 'N' 'U' 'P'}];
		append from A;
	close MWU;
	
	PROC PRINT data=MWU;
		title 'Hispanic & White';
	RUN;
QUIT;

/* Normaltest */
PROC UNIVARIATE data=BCHC_IM cibasic normaltest;
	where Race_Ethnicity = 'Asian/PI';
	var final_mortality;
	histogram final_mortality / normal;
 	qqplot final_mortality /normal(mu=est sigma=est);
	title 'Asian/PI';
RUN;
PROC UNIVARIATE data=BCHC_IM cibasic normaltest;
	where Race_Ethnicity = 'Black';
	var final_mortality;
	histogram final_mortality/ normal;
 	qqplot final_mortality /normal(mu=est sigma=est);
	title 'Black';
RUN;
PROC UNIVARIATE data=BCHC_IM cibasic normaltest;
	where Race_Ethnicity = 'Hispanic';
	var final_mortality;
	histogram final_mortality/ normal;
 	qqplot final_mortality /normal(mu=est sigma=est);
	title 'Hispanic';
RUN;
PROC UNIVARIATE data=BCHC_IM cibasic normaltest;
	where Race_Ethnicity = 'White';
	var final_mortality;
	histogram final_mortality/ normal;
 	qqplot final_mortality /normal(mu=est sigma=est);
	title 'White';
RUN;

/* mean tests */
PROC TTEST data=BCHC_IM;
	where R = 0 or R = 1;
	class R;
	var final_mortality;
	title 'Asian/PI & Black';
RUN;
PROC TTEST data=BCHC_IM;
	where R = 0 or R = 2;
	class R;
	var final_mortality;
	title 'Asian/PI & Hispanic';
RUN;
PROC TTEST data=BCHC_IM;
	where R = 0 or R = 3;
	class R;
	var final_mortality;
	title 'Asian/PI & White';
RUN;
PROC TTEST data=BCHC_IM;
	where R = 1 or R = 2;
	class R;
	var final_mortality;
	title 'Black & Hispanic';
RUN;
PROC TTEST data=BCHC_IM;
	where R = 1 or R = 3;
	class R;
	var final_mortality;
	title 'Black & White';
RUN;
PROC TTEST data=BCHC_IM;
	where R = 2 or R = 3;
	class R;
	var final_mortality;
	title 'Hispanic & White';
RUN;


PROC MIXED data=BCHC_IM method=TYPE3 cl;
	class Race_Ethnicity;
	model Mortality = Race_Ethnicity /solution cl;
	lsmeans Race_Ethnicity adjust=tukey;
RUN;

PROC MEANS data=BCHC_IM;
	var mortality final_mortality;
RUN;



/* Question 2: quantify how much of the variability between cities  */
/* can be explained by the race, and communicate this in a clear and */
/* practically interpretable way. */

/* Question 2a */
PROC MIXED data=BCHC_IM method=TYPE3 cl;	
	class Place;
	model final_mortality = /solution cl;
	random Place;
RUN;

/* Question 2b */
ods output "Solution for Fixed Effects" = Fixed_Anova;
PROC MIXED data=BCHC_IM method=TYPE3 cl;	
	class Place Race_Ethnicity;
	model final_mortality = Place /solution cl  outpm=RM outp=RC;
	random Race_Ethnicity;
RUN;


PROC IML;
	/* Estimate COV of random effect */
	S_g = 2704.55; 
	/* Estimate COV of residuals */
	S_res = 826.80;
	ICC = S_g / (S_g + S_res);
	
	A = ICC;
	create ICC from A [colname={'ICC'}]; 
	append from A;
	close ICC;
RUN;

/* • Normality of the residuals */
PROC FREQ data=RC; /* No ties */
	tables Resid;
RUN;

PROC UNIVARIATE data=RC normaltest;
	var resid;
	probplot resid /normal(mu=est sigma=est);
	histogram resid/normal;
RUN;

/* • Homogeneity of residual variance across groups */
/* We don’t have to test for homogeneity of  */
/* variance among the groups defined by random effects */
PROC GLM data=RC;
	class Place;
	model resid = Place;
	means Place/ hovtest=BF;
RUN;


/* • Normality of the random effects */


/* Question 2c */
PROC MIXED data=BCHC_IM method=TYPE3 cl;	
	class Place Race_Ethnicity Year;
	model final_mortality = Place /solution cl;
	random Race_Ethnicity Year;
RUN;
/* ICC = 0 */


/* Question 3: Which places are relevant for further analysis? */
PROC SQL;
	SELECT place, stdErr/estimate as cv 
	FROM Fixed_Anova
	ORDER BY cv desc;
RUN;


/* Hampel's Rule */
%MACRO Hampel(dataset, var, id);
	ods select none;
	PROC MEANS data=&dataset median;
		var &var;
     	output out=median median=median; 
    RUN;
    
    DATA hampel;
    	set &dataset(keep=&id &var); 
    	if _n_=1 then set median; 
    	abs_dev = abs(&var - median);
    RUN;
    
    PROC MEANS data=hampel median;
		var abs_dev;
		output out=abs_dev_median median=abs_dev_median;
	RUN;
	
	DATA hampel;
		set hampel;
		if _n_=1 then set abs_dev_median;
		abs_norm_val = abs_dev / abs_dev_median;
		if abs_norm_val <= 3.5 then delete;
	RUN;
	
	PROC SORT data=hampel;
		by descending abs_norm_val;
	RUN;
	ods select all;
    
    title1 "Outliers found using hampel’s rule.";
    title2 "Variable: &var";
	PROC PRINT data=hampel noobs label;
		var &id &var abs_norm_val;
      	label abs_norm_val="absolute normalized value (z_k)";
	RUN;
	title; 
%MEND;
