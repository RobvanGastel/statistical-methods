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
	ID = _n_; 
	if Race_Ethnicity = "White" then RE = 0;
	if Race_Ethnicity = "Black" then RE = 1;
	if Race_Ethnicity = "Asian/PI" then RE = 2;
	if Race_Ethnicity = "Hispanic" then RE = 3;
RUN;

/* Distributions of mortality per Ethnicity group */
/* TODO: Justify choise of test */

/* Per Group */
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

/* Total sample */
PROC UNIVARIATE data=BCHCT normal;
	histogram mortality/normal;
RUN;

PROC MEANS data=BCHCT skewness kurtosis n;
	var mortality;
RUN;

%Skewness_Kurtosis_test(skewness=1.1839264, kurtosis=19256820, n=148);
/* Not significant */


/* Dealing with NaN values */
PROC SQL;
	SELECT count(*) FROM BCHCT
	WHERE mortality is null;
RUN;
/* Amount of empty mortality rows is 16 */

/* In fact, simulation studies suggest that mean imputation */
/* is possibly the worst missing data handling method  */
/* available. Consequently, in no situation is mean  */
/* imputation defensible, and you should absolutely avoid  */
/* this approach. (Enders, C. K. (2010). Applied missing  */
/* data analysis. New York: Guilford.) */

/* TODO: Took MCAR approach and justify it */

/* Multiple Imputation approach */
/* Overview before filling in missing values */
PROC MEANS DATA=BCHCT mean var std median min max;
	class Race_Ethnicity;
RUN;

/* TODO: Maybe do this by group with where */
/* MCAR */
PROC SORT data=BCHCT;
	by Race_Ethnicity;
RUN;

/* We have 16 empty fields */
PROC MI data=BCHCT nimpute=5 out=BCHCT_MI seed=42 minimum=0;
	by Race_Ethnicity;
	var mortality year;
RUN;
/* Use as many variables for var as possible accoridng */
/* to the paper */
/* https://support.sas.com/rnd/app/stat/papers/multipleimputation.pdf */
/* Page 7 */

PROC PRINT data=BCHCT_MI;
RUN;

/* Overview after filling in missing values for mortality */
/* Increased the mean for asain/PI by +11, for hispanic by +3 */
PROC MEANS DATA=BCHCT_MI mean var std median min max;
	class Race_Ethnicity;
RUN;

/* Replace dataset */
DATA BCHCT;
	set BCHCT_MI;
RUN;

/* Removing outliers */
/* TODO: Justify choice of test for outliers */

/* From lecture: */
/* Perform your analyses also on the complete data  */
/* and report on the influence of excluding observations. */
/* • Sensitivity analysis. */
/* Discuss extreme observations with the data collector  */
/* when possible. */


PROC MEANS data=BCHCT;
RUN;

DATA BCHCT_mean;
	set BCHCT;
	mean = 163.9783784;
RUN;


/* All tests asumme normality which is why we need */
/* BoxCox transform for lambda = 0, total set is normal */
/* Both test for only 1 outliers and Grubbs is used */
/* More in practice. */
%Grubbs_test(dataset=BCHCT, var=mortality, id=ID);

%Doornbos_test(dataset=BCHCT, var=mortality, id=ID);

%Tukey_method(dataset=BCHCT, var=mortality, id=ID);

/* For Tukey method */
PROC BOXPLOT data=BCHCT_mean;
	plot Mortality*mean/boxstyle=schematic;
RUN;

%Hampel(dataset=BCHCT, var=mortality, id=ID);

/* Some of the outliers */
PROC SQL;
	SELECT * FROM BCHCT_n
	WHERE ID IN (138, 139, 140, 142, 143, 144, 34);
RUN;
/* Tukey and Hampel give the same outliers except for 34 */
/* 34 is an extreme value of Detriot. */
/* The other values \approx 138-144 \139 are all from */
/* San Antonio which indicate these are not outliers */

/* (BOX COX TRANSFORM BELOW SHOULD RUN FIRST) */
/* As we can transform data for our outlier test and know */
/* The mortality is normal after transform lambda = 0 */
/* We apply the hampel's rule as this one doesn't give us */
/* The place san antonio with the extreme values */

%Hampel(dataset=BCHC_BC, var=MT0, id=ID);
PROC SQL;
	SELECT * FROM BCHCT_n
	WHERE ID IN (83, 163);
RUN;
/* TODO: If we exclude these values we need a good */
/* justification and sensitivity analysis. */

/* Only finds NaN values */
%Tukey_method(dataset=BCHC_BC, var=MT0, id=ID);

/* Can't find outliers */
%Grubbs_test(dataset=BCHC_BC, var=MT0, id=ID);
%Doornbos_test(dataset=BCHC_BC, var=MT0, id=ID);

/* Transforming the data for normality */
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

/* ANOVA on transformed data MINUS12 */
PROC MIXED data=BCHC_BC method=TYPE3 cl;
	class RE;
	model MTMINUS12 = /solution cl outpm=RM outp=RC;
	random RE /solution;
RUN;

/* TODO: Two-Way ANOVA  */
PROC MIXED data=BCHCT method=TYPE3 cl;
	class RE;
	model Mortality = /solution cl DDFM=SAT outpm=RM outp=RC;
	random RE /solution;
	LSMEANS RE /CL;
RUN;


/* CHECK ANOVA ASSUMPTIONS: */
/* • Normality of the residuals */
/* • Homogeneity of residual variance across groups */
/* • Normality of the random effects */








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

%MACRO Skewness_Kurtosis_test(skewness=, kurtosis=, n=);
	DATA approx;
		N=&n;
		G1=&skewness;
		G2=&kurtosis;
		b1=(N-2)*G1/(sqrt(N*(N-1)));
		b2=G2*((N-2)*(N-3))/((N+1)*(N-1))+3*(N-1)/(N+1);

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
	RUN;
	
	PROC PRINT data=approx noobs;
	RUN;
%MEND;
