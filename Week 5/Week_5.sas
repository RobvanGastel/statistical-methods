LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 5.1 */
/* you assumed that the age of the mother (AGEM) */
/* is normally distributed. In this question you */
/* will verify these assumptions. */

DATA WEEK5_Q1;
	set SASDATA.IVF;
RUN;

/* a) */
/* TODO: Doesn't generate same results as answer */
PROC UNIVARIATE data=WEEK5_Q1 normaltest;
	var AGEM;
	histogram AGEM/normal;
	probplot AGEM/normal(mu=est sigma=est);
RUN;

/* Shapiro-Wilk */
/* Test statistic = 0.992907 and p-value = 0.0011 */
/* Conclusion: */
/* Is sensitive to non symmetric departures from normality. */
/* This distribution seems relatively symmetri. */

/* b) */
/* Kolmogorov-Smirnov */
/* Test statistic = 0.040049 and p-value < 0.0100 */
/* Conclusion: */
/* As it compares the vertical distance (supremum) between */
/* the distributions we can see that there is a relatively */
/* large distance between F_n and F. */

/* c) */
/* Cramer-von Mises */
/* Test statistic = 0.178199 and p-value = 0.0099 */
/* Conclusion: */
/* Puts more weight on the tails as the tails deviate less */
/* it gets a higher score then the other tests as the tails */
/* Seem relatively normal. */

/* d) */
/* Anderson-Darling */
/* Test statistic = 1.271699 and p-value < 0.0050 */
/* Conclusion: */
/* As it looks at the integral of the distance between F_n */
/* and F. This has the smallest p-value of all the tests. */

/* e) */
/* TODO */

/* f) and g) */
%MACRO Skewness_Kurtosis_test(skewness=, kurtosis=, n=);
	/* if n >= 100, b_1 \approx Y^hat_1 and */
	/* b_2 \approx Y^hat_2 + 3 */
	DATA approx;
		N=&n;
		G1=&skewness;
		G2=&kurtosis;
		b1=(N-2)*G1/(sqrt(N*(N-1))); /* Sample Skewness */
		b2=G2*((N-2)*(N-3))/((N+1)*(N-1))+3*(N-1)/(N+1); /* Sample Kurtosis */

		Cn=(3*(N**2+27*N-70)*(N+1)*(N+3))/((N-2)*(N+5)*(N+7)*(N+9));
		Wn2=-1+sqrt(2*(Cn-1));
		
		Alphan=sqrt(2/(Wn2-1));
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

/* Skewness = -0.1773114, Kurtosis = -0.1983322 */
/* and n = 759 */

%Skewness_Kurtosis_test(skewness=-0.1773114, kurtosis=-0.198322, n=759);

/* h) */
/* TODO */

/* Question 5.2 */
/* BW is approximately normal at lambda = 2 */
DATA WEEK5_Q2;
	set SASDATA.IVF;
RUN;

/* a) */
/* By use of the grubbs test we determine there are */
/* no outliers. */

/* For n =< 100 look up the critical value in the table */
/* After sorting if the first value > approx or absolute */
/* value, we reject H0 and determine there is a outlier. */
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

%Grubbs_test(dataset=WEEK5_Q2, var=BW, id=ID);

PROC UNIVARIATE data=WEEK5_Q2 normaltest;
	var BW;
	histogram BW/normal;
	probplot BW/normal(mu=est sigma=est);
RUN;

/* Shows alot of ties (freq > 1) */
PROC FREQ data=WEEK5_Q2;
	tables BW;
RUN;

/* Anderson-Darling */
/* Test statistic = 1.125654 and p-value = 0.0063 */

/* b) */
PROC MEANS data=WEEK5_Q2 mean std n;
	var BW;
	output out=ds_out mean=mean median=median std=std n=n;
RUN;
	
DATA WEEK5_Q2;
	set WEEK5_Q2;
	if _n_= 1 then set ds_out;
	drop _TYPE_ _FREQ_;
RUN;
	
/* Box plots */
PROC BOXPLOT data=WEEK5_Q2;
	plot BW*MEAN/boxstyle=schematic;
RUN;

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
		label lower="Lower limit" upper="Upper limit" iqr="Interquartile range (IQR)";
	RUN; 
	title;
%MEND;

/* ID = 98 and ID = 294 are outliers by the Tukey's method */
%Tukey_method(dataset=WEEK5_Q2, var=BW, id=ID);

/* Dataset without outliers */
DATA WEEK5_Q2_b;
	set SASDATA.IVF;
	if ID = 98 then delete;
	if ID = 294 then delete;
RUN;
	
PROC UNIVARIATE data=WEEK5_Q2_b normaltest;
	var BW;
	histogram BW/normal;
	probplot BW/normal(mu=est sigma=est);
RUN;
/* The Test statistics and p-values remain mostly the same */

/* c) */
/* Box-Cox transform with lambda = 2. */
DATA WEEK5_Q2;
	set SASDATA.IVF;
	BWPLUS2 = (0.5)*(BW**(2) - 1);
RUN;

/* H0 couldn't be rejected, no outliers */
%Grubbs_test(dataset=WEEK5_Q2, var=BWPLUS2, id=ID);

/* H0 couldn't be rejected, no outliers */
%Tukey_method(dataset=WEEK5_Q2, var=BWPLUS2, id=ID);

PROC UNIVARIATE data=WEEK5_Q2 normaltest;
	var BWPLUS2;
	histogram BWPLUS2/normal;
	probplot BWPLUS2/normal(mu=est sigma=est);
RUN;
/* Because there are some outliers we use Anderson-Darling test */

/* Question 5.3 */
DATA WEEK5_Q3;
	set SASDATA.IVF;
	where PER = 4;
	GAL = log(44 - GA);
RUN;

/* a) */
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

/* Perform outlier test before transform */
%Doornbos_test(dataset=WEEK5_Q3, var=GA, id=ID);

/* Outliers: id = 174 */
DATA WEEK5_Q3;
	set WEEK5_Q3;
	if ID = 174 then delete;
RUN;

PROC UNIVARIATE data=WEEK5_Q3 normaltest;
	var GA;
	histogram GA/normal;
	probplot GA/normal(mu=est sigma=est);
RUN;
/* SHOULD'VE BEEN ANDERSON-DARLING */
/* As we have a slight right skew we use the shapiro-Wilk */
/* test, with T statistic = 0.9912354 and p-value = 0.1378 */
/* ALSO DO TEST FOR TIES */

/* b) */
DATA WEEK5_Q3;
	set SASDATA.IVF;
	where PER = 4;
	GAL = log(44 - GA);
RUN;

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

/* Outliers found */
%Hampel(dataset=WEEK5_Q3, var=GA, id=ID);

/* Remove the outliers from the dataset */
PROC SQL;
	CREATE TABLE WEEK5_Q3_b AS
		SELECT * FROM WEEK5_Q3
		WHERE ID NOT IN (SELECT ID FROM hampel);
RUN;

PROC UNIVARIATE data=WEEK5_Q3 normaltest;
	var GA;
	histogram GA/normal;
	probplot GA/normal(mu=est sigma=est);
RUN;
/* Anderson-Darling test */
/* T statistic = 0.519412 and p-value = 0.1940 */

/* c) */
DATA WEEK5_Q3;
	set SASDATA.IVF;
	where PER = 4;
	GAL = log(44 - GA);
RUN;

%Doornbos_test(dataset=WEEK5_Q3, var=GAL, id=ID);
%Hampel(dataset=WEEK5_Q3, var=GAL, id=ID);
/* Hampel tests finds some outliers */

PROC UNIVARIATE data=WEEK5_Q3 normaltest;
	var GA;
	histogram GA/normal;
	probplot GA/normal(mu=est sigma=est);
RUN;

/* Question 5.4 */
DATA COAG; 
	input Patient C K@@;
	datalines;
	1 120 132 8 145 133 15 117 123
	2 114 116 9 120 123 16 125 108
	3 129 135 10 129 116 17 136 131 
	4 128 115 11 126 127 18 151 119
	5 155 134 12 136 140 19 130 129 
	6 105 56 13 135 140 20 136 124
	7 114 114 14 125 114 21 113 112
	;
RUN;

/* K numeric value 1, C is numeric value 2 */
PROC SQL;
	CREATE TABLE COAG_T AS
		SELECT c.Patient, c.K as Value, 1 as Type,
		(g.K + g.C)/2 as mean
		FROM COAG c, COAG g
		WHERE c.Patient = g.Patient
		UNION
		SELECT c.Patient, c.C, 2,
		(g.K + g.C)/2 as mean
		FROM COAG c, COAG g
		WHERE c.Patient = g.Patient;
RUN;

/* a) */
PROC MIXED data=COAG_T method=TYPE3 cl;
	class patient type;
	model value = /solution cl outpm=RM outp=RC;
	random patient;
RUN;

/* Shows alot of ties (freq > 1) */
PROC FREQ data=RC;
	tables Resid;
RUN;

PROC UNIVARIATE data=RC normaltest;
	var Resid;
	histogram Resid/normal;
	probplot Resid/normal(mu=est sigma=est);
RUN;

/* Anderson-Darling */
/* Test Statistic = 1.158151 and p-value <0.0050 */

/* b) */
/* TODO */

/* c) */
/* TODO */

/* Question 5.5 */

/* a) */
/* TODO */
/* b) */
/* TODO */

/* Question 5.6 */
DATA WEEK5_Q6; 
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
PROC UNIVARIATE data=WEEK5_Q6 normaltest;
	histogram value/normal;
	probplot value/normal(mu=est sigma=est);
RUN;

/* b) */
/* Add ID */
DATA WEEK5_Q6;
	set WEEK5_Q6;
	ID = _n_;
RUN;

%Doornbos_test(dataset=WEEK5_Q6, var=value, id=ID);
%Grubbs_test(dataset=WEEK5_Q6, var=value, id=ID);
%Hampel(dataset=WEEK5_Q6, var=value, id=ID);
%Tukey_method(dataset=WEEK5_Q6, var=value, id=ID);

/* Question 5.8 */
/* a)-d) */
/* Are all implemented */


/* Question 5.9 */
DATA WEEK5_Q9;
	do i = 1 to 100000;
	  X = rand('normal', 0, 1);
	  output;
	end;
RUN;

/* a) */

/* We can see its approximately standard normally */
/* distributed. */
PROC MEANS data=WEEK5_Q9;
RUN;

%Tukey_method(dataset=WEEK5_Q9, var=X, id=i);

PROC MEANS data=Tukey n;
RUN;
/* For n = 1.000.000 we observed 711 outliers */
/* which is 0,711% of the data was classifed as */
/* outliers. */
/* Which is close to the number we expected which is */
/* 4,30%. */

/* b) */
DATA WEEK5_Q9;
	do i = 1 to 100000;
	  X = rand('exponential', 1);
	  output;
	end;
RUN;

%Tukey_method(dataset=WEEK5_Q9, var=X, id=i);

PROC MEANS data=Tukey n;
RUN;
/* 4.909% of the data is classified as outlier and the */
/* actual value we expected is 6,61%. */
