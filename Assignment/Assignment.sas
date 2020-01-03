LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";


DATA BCHC;
	set SASDATA.BCHC;
RUN;


/* Exploration data */
PROC MEANS data=BCHC n mean std var median min max;
	class Race_Ethnicity;
RUN;

DATA BCHCT;
	set BCHC;
	if Race_Ethnicity = "White" then RE = 0;
	if Race_Ethnicity = "Black" then RE = 1;
	if Race_Ethnicity = "Asian/PI" then RE = 2;
	if Race_Ethnicity = "Hispanic" then RE = 3;
RUN;


PROC SQL;
	SELECT count(*) FROM BCHCT
	WHERE mortality is null;
RUN;
/* Amount of empty mortality rows is 16 */

/* Distributions of mortality per Ethnicity group */
/* TODO: Could extend with Skewness and Kurtosis test */
PROC UNIVARIATE data=BCHCT normal;
	where RE = 0;
	var mortality;
	histogram mortality/normal;
RUN;

PROC UNIVARIATE data=BCHCT normal;
	where RE = 1;
	var mortality;
	histogram mortality/normal;
RUN;

PROC UNIVARIATE data=BCHCT normal;
	where RE = 2;
	var mortality;
	histogram mortality/normal;
RUN;

PROC UNIVARIATE data=BCHCT normal;
	where RE = 3;
	var mortality;
	histogram mortality/normal;
RUN;

PROC UNIVARIATE data=BCHCT normal;
	histogram mortality/normal;
RUN;

/* Box Cox transform */
/* λ ∈ {−2,−1/2,0,1/2,2} */
DATA BCHC_BC; 
	set BCHCT;
	MTMINUS2 = (-1/2)*(mortality**-2 -1); 
	MTMINUS1 = (-1)*(mortality**-1 -1); 
	MTMINUS12 = (-2)*(mortality**-(0.5)-1); 
	MT0 = log(mortality);
	MTPLUS12 = (2)*(mortality**(1/2) -1); 
	MTPLUS2 = (0.5)*(mortality**(2) -1);
RUN;

PROC UNIVARIATE data=BCHC_BC normaltest;
	where RE = 0;
	var Mortality MTMINUS2 MTMINUS1 MTMINUS12
		MT0 MTPLUS12 MTPLUS2;
   	histogram Mortality /normal;
	histogram MTMINUS2 /normal; 
	histogram MTMINUS1 /normal;
   	histogram MTMINUS12 /normal;
   	histogram MT0 /normal;
  	histogram MTPLUS12 /normal;
   	histogram MTPLUS2 /normal;
RUN;
/* Significant MTMINUS2, MTMINUS1, MTMINUS12 */
/* TODO: Justify a test and pick one */

PROC UNIVARIATE data=BCHC_BC normaltest;
	where RE = 1;
	var Mortality MTMINUS2 MTMINUS1 MTMINUS12
		MT0 MTPLUS12 MTPLUS2;
   	histogram Mortality /normal;
	histogram MTMINUS2 /normal; 
	histogram MTMINUS1 /normal;
   	histogram MTMINUS12 /normal;
   	histogram MT0 /normal;
  	histogram MTPLUS12 /normal;
   	histogram MTPLUS2 /normal;
RUN;
/* Significant MTMINUS12, MT0, MTPLUS12 (Not fully MTMINUS1) */
/* TODO: Justify a test and pick one */


PROC UNIVARIATE data=BCHC_BC normaltest;
	where RE = 2;
	var Mortality MTMINUS2 MTMINUS1 MTMINUS12
		MT0 MTPLUS12 MTPLUS2;
   	histogram Mortality /normal;
	histogram MTMINUS2 /normal; 
	histogram MTMINUS1 /normal;
   	histogram MTMINUS12 /normal;
   	histogram MT0 /normal;
  	histogram MTPLUS12 /normal;
   	histogram MTPLUS2 /normal;
RUN;
/* Significant MTMINUS2, MTMINUS1, (Not fully MTMINUS12)  */
/* TODO: Justify a test and pick one */

PROC UNIVARIATE data=BCHC_BC normaltest;
	where RE = 3;
	var Mortality MTMINUS2 MTMINUS1 MTMINUS12
		MT0 MTPLUS12 MTPLUS2;
   	histogram Mortality /normal;
	histogram MTMINUS2 /normal; 
	histogram MTMINUS1 /normal;
   	histogram MTMINUS12 /normal;
   	histogram MT0 /normal;
  	histogram MTPLUS12 /normal;
   	histogram MTPLUS2 /normal;
RUN;
/* Not significant for any transformation highest yield is */
/* MTMINUS12 */
/* TODO: Justify a test and pick one */

/* MTMINUS12 looks OK */



/* Multiple groups (ANOVA) on not transformed data */
PROC MIXED data=BCHCT method=TYPE3 cl;
	class RE;
	model Mortality = /solution cl outpm=RM outp=RC;
	random RE /solution;
RUN;

/* ANOVA on transformed data */
PROC MIXED data=BCHC_BC method=TYPE3 cl;
	class RE;
	model MTMINUS12 = /solution cl outpm=RM outp=RC;
	random RE /solution;
RUN;

/* TODO: Two-Way ANOVA  */


/* CHECK ANOVA ASSUMPTIONS: */
/* • Normality of the residuals */
/* • Homogeneity of residual variance across groups */
/* • Normality of the random effects */




/* Removing outliers */
/* TODO: Justify choice of test for outliers */
DATA BCHCT_n;
	set BCHCT;
	ID = _n_;
RUN;

%Grubbs_test(dataset=BCHCT_n, var=mortality, id=ID);
%Doornbos_test(dataset=BCHCT_n, var=mortality, id=ID);
%Tukey_method(dataset=BCHCT_n, var=mortality, id=ID);






/* Measure dependence */
PROC CORR data=BCHCT spearman kendall pearson 
		fisher(biasadj=no)
		plots=scatter(ellipse=none);
	var RE mortality;
RUN;
/* Rejecting all H0 */
/* There is some form of dependence? */

/* MACRO's */
%MACRO Grubbs_test(dataset, var, id, alpha=0.05);
    ods select none;
	PROC MEANS data=&dataset mean var n;
		var &var;
	    output out=out mean=mean var=var n=n; 
	RUN;
	
	DATA outliers;
		set &dataset(keep=&id &var);
		if _n_=1 then set out;
		/* statistic */
		u = abs((&var - mean) / sqrt(var));
		/* critical value */
		t = quantile("t", &alpha / (2*n), n-2);
		c = (n-1) * sqrt(t**2 / (n * (t**2 + n - 2))); 
		/* check if this is an outlier */
		if(u > c) then outlier = "yes"; 
		else outlier = "no"; /* p-value */
		u_inv = u*sqrt((n-2)*n) / sqrt(1 - (u**2-(n-2))*n); 
		p_value = min(2*n*(1-cdf("t", u_inv, n-2)), 1); 
		keep &id p_value u c outlier;
	RUN;
	
	PROC SORT data=outliers; 
		by descending u;
	RUN;
	ods select all;
	              
	title1 "Grubbs test";
	title2 "Variable: &var";
	PROC PRINT data=outliers(obs=1) label noobs;
		var &id p_value u c outlier; 
		label &id="id" p_value="p-value"
		u="statistic (|u|)" c="critical value" 
		outlier="is outlier?";
    RUN;
    title;
%MEND;

/* Tukey's method */
%MACRO Tukey_method(dataset, var, id);
	ods select none;
	PROC MEANS data=&dataset median p25 p75;
		var &var;
		output out=quartiles p25=p25 p75=p75;
	RUN;
	
	DATA tukey;
		set &dataset(keep=&id &var);
		if _n_=1 then set quartiles;
		iqr = p75-p25;
		lower = p25 - 1.5*IQR;
		upper = p75 + 1.5*IQR;
		if &var >= lower and &var <= upper then delete;
	RUN;
	ods select all;
	
	title1 "Outliers found using tukey’s method."; 
	title2 "Variable: &var";
	PROC PRINT data=tukey noobs label;
		var &id &var lower upper iqr; 
		label lower="Lower limit" upper="Upper limit"
		iqr="Interquartile range (IQR)";
	RUN; 
	title;
%MEND;

/* Doornbos test */
%MACRO Doornbos_test(dataset, var, id, alpha=0.05);
	ods select none;
	PROC MEANS data=&dataset mean var n;
		var &var;
		output out=out mean=mean var=var n=n; 
	RUN;
	
	DATA outliers;
		set &dataset(keep=&id &var);
		if _n_=1 then set out;
		/* leave-one-out variance */ 
		var_loo = ((n - 1) / (n - 2)) * var - (n / ((n - 1)*(n - 2)))*(&var - mean)**2;
		/* statistics */
		w = abs((&var - mean) / sqrt(var_loo * (n - 1) / n)); 
		/* critical value */
		c = quantile("t", 1 - &alpha / (2*n), n-2);
		/* check if this is an outlier */
		if(w > c) then outlier = "yes"; 
		else outlier = "no"; /* p-value */
		p_value = min(2*n*(1-cdf("t", w, n-2)), 1);
		keep &id p_value w c outlier; 
	RUN;
	
	PROC SORT data=outliers; 
		by descending w;
	RUN;
	ods select all;
	
	title1 "Doornbos test";
	title2 "Variable: &var";
	PROC PRINT data=outliers(obs=1) label noobs;
		var &id p_value w c outlier;
		label &id="id" p_value="p-value" 
		w="statistic (|W|)" c="critical value" 
		outlier="is outlier?";
	RUN;
	title; 
%MEND;

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
