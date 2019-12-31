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
	probplot AGEM/normal;
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
	BWPLUS2 = (0.5)*(BW**(2) - 1);
RUN;

/* a) */
/* By use of the grubbs test we determine there are */
/* no outliers. */

/* For n =< 100 look up the critical value in the table */
/* After sorting if the first value > approx or absolute */
/* value, we reject H0 and determine there is a outlier. */
%MACRO Grubbs_test(ds=, var=);
	PROC MEANS data=&ds mean std n;
		var &var;
		output out=ds_out mean=mean median=median std=std n=n;
	RUN;
	
	DATA &ds;
		set &ds;
		if _n_=1 then set ds_out;
		drop _TYPE_ _FREQ_;
	RUN;
	
	/* By sorting the value and  */
	DATA Grubbs;
		set &ds;
		U = (&var - mean)/std;
		Grubbs=abs(U);
		
		/* For n =< 100 look up the critical value in */
		/* the table */
		/* 	C_onesided_exact= 2.56; */
		/* 	C_twosided_exact= 2.71;	 */
		
		t = quantile("t", 0.05 / (2*N), N-2);
		u_inv = u*sqrt((n-2)) / sqrt(n-1-u**2);
		
		C_twosided_approx = (n-1) * sqrt(t**2 / (n * (t**2 + n - 2)));
		p_twosided_approx = min(2*n*min(1-cdf("t", u_inv, n-2),cdf("t", u_inv, n-2)), 1);
	RUN;
	
	PROC SORT data=Grubbs;
		by descending Grubbs;
	RUN;
	
	PROC PRINT data=Grubbs;
	RUN;
%MEND Grubbs_test; 

%Grubbs_test(ds=WEEK5_Q2, var=BWPLUS2);

%Grubbs_test(ds=WEEK5_Q2, var=BW);

PROC UNIVARIATE data=WEEK5_Q2 normaltest;
	var BWPLUS2;
	histogram BWPLUS2/normal;
	probplot BWPLUS2/normal;
RUN;

/* Anderson-Darling */
/* Test statistic = 1.125654 and p-value = 0.0063 */

/* b) */
PROC MEANS data=WEEK5_Q2 mean std n;
	var BW;
	output out=ds_out mean=mean_bw median=median_bw std=std_bw n=n;
RUN;
	
DATA WEEK5_Q2;
	set WEEK5_Q2;
	if _n_=1 then set ds_out;
	drop _TYPE_ _FREQ_;
RUN;
	
/* Box plots */
PROC BOXPLOT data=WEEK5_Q2;
	plot BW*MEAN_BW/boxstyle=schematic;
RUN;

PROC BOXPLOT data=WEEK5_Q2;
	plot BWPLUS2*MEAN/boxstyle=schematic;
RUN;

/* Tukey's method */

/* (Reset Dataset) */
DATA WEEK5_Q2;
	set SASDATA.IVF;
	BWPLUS2 = (0.5)*(BW**(2) - 1);
RUN;

PROC MEANS data=WEEK5_Q2 mean std n;
	var BW;
	output out=ds_out mean=mean median=median std=std n=n p25=p25 p75=p75;
RUN;

DATA TUKEY;
	set WEEK5_Q2;
	IQR=p75-p25;
	LOWERT = p25 - 1.5*IQR;
	UPPERT = p75 + 1.5*IQR;
RUN;












