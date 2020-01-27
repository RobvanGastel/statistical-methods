LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 4.1 */
DATA WEEK4_Q1;
	input supplement$ iq@@; 
	datalines;
suppB 104.3 
suppB 99.0 
suppB 132.4 
suppB 109.4
suppB 112.4 
suppB 101.9 
suppB 100.7 
suppB 100.5
suppB 105.3 
suppB 110.5
suppB 112.5 
suppB 114.0
suppB 98.8
suppB 98.9
suppB 97.0
suppB 112.1
suppB 114.8
suppB 100.6
suppB 110.7
suppB 119.3
suppA 108.1
suppA 97.0 
suppA 106.7
suppA 96.7 
suppA 105.4 
suppA 106.3 
suppA 99.7
suppA 109.5 
suppA 99.8 
suppA 102.7 
suppA 106.3
suppA 97.1
suppA 105.6
suppA 110.0
suppA 108.4
suppA 106.3
suppA 93.7
suppA 107.7
suppA 97.7
suppA 107.3
	;
RUN;

/* a) */
/* Fixed model, */
/* 𝒀_ij = 𝝁 + 𝜶_i + 𝒆_ij */

/* b) */
/* H0:𝛼_1=𝛼_2=⋯=𝛼_m=0 vs H1:∃𝑖,𝛼_𝑖≠0 */
PROC MIXED data=WEEK4_Q1 method=TYPE3 cl; 
	class supplement;
	model iq = supplement /solution cl;
	lsmeans supplement /diff cl;
RUN;
/* For 0 the estimate is -4.1550 and p-value = 0.0747 */
/* Thus, we can't reject H0. */

/* c) */
/* Create the CI for mu */
PROC UNIVARIATE data=WEEK4_Q1 cibasic;
RUN;

/* TODO: This estimation is off by +-0.06 for mu */
*mu;
PROC UNIVARIATE data=WEEK4_Q1 cibasic;
	var iq;
RUN;

*a_1;
PROC IML;
	n = 40;
	alpha=0.05;
	y_bar = 2.0825;
	s = 5.77953;
	qt = quantile("t", alpha/2, n-1);
	
	UCL = y_bar + qt * (s/(sqrt(n)));
	LCL = y_bar - qt * (s/(sqrt(n)));
	
	PRINT(UCL||LCL);
RUN;

*a_1;
PROC IML;
	n = 40;
	alpha=0.05;
	y_bar = -2.0775;
	s = 5.77923;
	qt = quantile("t", alpha/2, n-1);
	
	UCL = y_bar + qt * (s/(sqrt(n)));
	LCL = y_bar - qt * (s/(sqrt(n)));
	
	PRINT(UCL||LCL);
RUN;

/* d) */
PROC TTEST data=WEEK4_Q1;
	class supplement;
	var iq; 
RUN;
/* Test statistic = -1.83 and p-value = 0.0747 */
/* So this is the same as the ANOVA. */

/* e) */
/* If value > 1 there are ties */
PROC FREQ data=WEEK4_Q1;
	table iq;
RUN;

/* Use WRS to calculate Kruskal-Wallis test */
PROC NPAR1WAY data=WEEK4_Q1 wilcoxon correct=NO;
	class supplement;
	var iq;
RUN;
/* Test statistic = 2.5916 and p-value = 0.1074 */
/* Thus we can't reject H0 */

/* f) */
/* As normality is not violated the ANOVA is a lot */
/* more powerful than Kruskal-Wallis. */

/* g) */
DATA WEEK4_Q1_g;
	input supplement$ iq@@; 
	datalines;
suppB 104.3 
suppB 99.0 
suppB 132.4 
suppB 109.4
suppB 112.4 
suppB 101.9 
suppB 100.7 
suppB 100.5
suppB 105.3 
suppB 110.5 
suppA 97.0 
suppA 106.7
suppA 96.7 
suppA 105.4 
suppA 106.3 
suppA 99.7
suppA 109.5 
suppA 99.8 
suppA 102.7 
suppA 106.3
suppB 112.5 
suppB 114.0
suppB 98.8
suppB 98.9
suppB 97.0
suppB 112.1
suppB 114.8
suppB 100.6
suppB 110.7
suppB 119.3
suppA 108.1
suppA 97.1
suppA 105.6
suppA 110.0
suppA 108.4
suppA 106.3
suppA 93.7
suppA 107.7
suppA 97.7
suppA 107.3
suppC 103.3
suppC 104.0
suppC 117.5
suppC 119.0
suppC 135.4
suppC 113.4
suppC 103.8
suppC 103.9
suppC 115.4
suppC 106.9
suppC 102.0
suppC 117.1
suppC 105.7
suppC 105.5
suppC 119.8
suppC 105.6
suppC 110.3
suppC 115.5
suppC 115.7
suppC 124.3
	;
RUN;

PROC UNIVARIATE data=WEEK4_Q1_g;
RUN;
/* Gives mean 107.8533 */

PROC MIXED data=WEEK4_Q1_g method=TYPE3 cl; 
	class supplement;
	model iq = supplement /solution cl; 
	lsmeans supplement /diff cl;
RUN;

/* h) */
/* The paired T-test is inline with results */
/* Take s_pooled as H_0 of the ANOVA expects the sum */
/* of a_i to be equal to 0. */
PROC TTEST data=WEEK4_Q1_g;
	where supplement in ("suppB", "suppC");
	class supplement;
	var iq; 
RUN;

PROC TTEST data=WEEK4_Q1_g;
	where supplement in ("suppB", "suppA");;
	class supplement;
	var iq; 
RUN;

PROC TTEST data=WEEK4_Q1_g;
	where supplement in ("suppA", "suppC");;
	class supplement;
	var iq; 
RUN;
/* This results in a different standard error */

/* Question 4.4 */
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
ods output "Covariance Parameter Estimates" = Cov;
PROC MIXED data=COAG_T method=TYPE3 cl;
	class patient type;
	model value = /solution cl;
	random patient;
RUN;

/* b) */
PROC PRINT data=Cov;
RUN;

PROC IML;
	/* Estimate COV of random effect */
	S_g = 122.21; 
	/* Estimate COV of residuals */
	S_res = 123.81;
	ICC = S_g / (S_g + S_res);
	
	A = ICC;
	create ICC from A [colname={'ICC'}]; 
	append from A;
	close ICC;
RUN;

/* c) */
PROC SORT data=COAG_T;
	by type patient;
RUN;

ods output SolutionR=COAG_Solution;
PROC MIXED data=COAG_T method=TYPE3 cl;
	class patient type;
	model value = /solution cl outp=RC outpm=RM;
	random patient /solution;
RUN;
/* Patient 9, EBLUP = -2.0071 */

/* d) */
/* Output table */
PROC PRINT data=RM;
PROC PRINT data=RC;

DATA RM;
	set RM;
	ID = _n_;
RUN;

DATA RC;
	set RC;
	ID = _n_;
RUN;

DATA RM;
	set RM; 
	RESIDM=RESID; 
	PredM=Pred; 
	drop RESID Pred; 
RUN;

PROC SORT data=RM;
	by Patient;
RUN;


DATA PREDMERGE; 
	merge RM RC; 
	by ID TYPE; 
RUN;

DATA PREDMERGE; 
	set PREDMERGE;
	EBLUP=PRED-PREDM;
RUN;

%MACRO Runs_test(data=, var=, alpha=);
	PROC IML;
		use &data;
		read all var {&var};
		close &data;
		
		X=&var;
		n=nROW(X);
		MED=median(X);
		
		XC=X;
		do i=1 to n by 1;
			if (XC[i] >= MED) then XC[i]=1;
			else XC[i]=0;
		end;
		
		n1C=sum(XC);
		n2C=n-n1C;
		
		RC=1;
		do i=2 to n by 1;
			if(XC[i] ^= XC[i-1]) then RC=RC+1;
		end;
		
		MUC=1+(2*n1C*n2C)/(n1C+n2C);
		VARC=2*n1C*n2C*(2*n1C*n2C-n1C-n2C)/((n1C+n2C-1)*(n1C+n2C)**2);
		
		SC=(RC-MUC)/sqrt(VARC);
		TC=quantile('NORMAL',&alpha/2);
		TCU=quantile('NORMAL',1-&alpha/2);
		PC=(1-cdf('NORMAL',abs(SC)))*2;
		
		XUC=repeat(0,n-1,1);
		TIES=0;
		do i=1 to (n-1) by 1;
			if (X[i+1] > X[i]) then XUC[i]=1;
			if (X[i+1] < X[i]) then XUC[i]=0;
			if (X[i+1] = X[i]) then XUC[i]=XUC[i-1];
			if (X[i+1] = X[i]) then TIES=TIES+1;
		end;
		
		RUC=1;
		do i=2 to (n-1) by 1;
			if(XUC[i] ^= XUC[i-1]) then RUC=RUC+1;
		end;
		
		MUUC=(2*(n-TIES)-1)/3;
		VARUC=(16*(n-TIES)-29)/90;
		
		SUC=(RUC-MUUC)/sqrt(VARUC);
		TUC=quantile('NORMAL',&alpha/2);
		TUCU=quantile('NORMAL',1-&alpha/2);
		PUC=(1-cdf('NORMAL',abs(SUC)))*2;
		
		A = RC||MUC||sqrt(VARC)||PC||SC||TC||TCU||n;
		create CondRuns from A [colname={'runs','mean runs','std runs','p-value','Normalized statistic','Critical values L','Critical values U','n'}];
		append from A; 
		
		B = RUC||MUUC||sqrt(VARUC)||PUC||SUC||TUC||TUCU||n-TIES;
		create UncondRuns from B [colname={'runs','mean runs','std runs','p-value','Normalized statistic','Critical values L','Critical values U','n'}]; 
		append from B;
		
		C = ties;
		create Colties from C [colname={'Ties'}];
		append from C;
		close UncondRuns CondRuns Colties;
				
		PROC PRINT data=CondRuns noobs;
			title "Median based (conditional) runs test";
		RUN;
		
		PROC PRINT data=UncondRuns noobs;
			title "(unconditional) runst test for serial randomness";
		RUN;
		
		PROC PRINT data=Colties noobs;
			title "Number of ties for both";
		RUN;
	quit;
%MEND;

/* Look at the conditional runs test */
/* As we're not interested in serial correlation */
%Runs_test(data=PREDMERGE, 
		 var=EBLUP,
		 alpha=0.05);

/* e) */
/* Through the randomness tests we check the  */
/* independence of the residuals.  */

/* f) */
/* TODO */
PROC SORT data=COAG_SOLUTION;
	by Patient;
RUN;

%Runs_test(data=COAG_Solution, 
		 var=Estimate,
		 alpha=0.05);


PROC SORT data=predmerge;
	by EBLUP;
RUN;		 

%Runs_test(data=PREDMERGE, 
		 var=RESIDM,
		 alpha=0.05);


/* Question 4.6 */
DATA RCT;
	set SASDATA.RCT;
	where TIME = 1;
RUN;

/* a) */
/* Yij = mu + alpha_i + e_ij */

/* b) */
PROC MIXED data=RCT method=TYPE3 cl;
	class ID;
	model RESP = /solution cl outpm=RM outp=RC;
	random ID /solution;
RUN;
/* EBLUPS Patient 1 = -0.01764 and Patient 2 = -0.5176 */

PROC MIXED data=rct method=type3;
	class ID;
	model RESP = / outp=RC;
	random ID /solution;
RUN;

*c;
PROC SORT data=RC;
	by center;
RUN;

%Runs_test(data=RC, var=Resid, alpha=0.05);

/* Question 4.7 */
DATA WEEK4_Q7;
	input input@@; 
	datalines;
	104.3 132.4 112.4 100.7 105.3 99.0 109.4 101.9 100.5 110.5
	112.5 98.8 97.0 114.8 110.7 114.0 98.9 112.1 100.6 119.3
	;
RUN;

/* a) */
%runs_test(data=WEEK4_Q7, var=input, alpha=0.05);
/* conditional test, p-value = 0.64590 */
/* Thus can't reject, H0 */

/* b) */
%runs_test(data=WEEK4_Q7, var=input, alpha=0.05);
/* conditional test, p-value = 0.26603 */
/* Thus can't reject, H0 */

/* c) & d) */
/* TODO: */
PROC RANK data=WEEK4_Q7 out=WEEK4_Q7_R;
      var input;
      ranks r_input;
RUN;

%runs_test(data=WEEK4_Q7_R, var=r_input, alpha=0.05);
/* Values stayed the same */

/* Question 4.8 */
/* TODO: Implement MACRO */
/* %macro sercor(dataset,var1,var2,lag);  */
/* --- YOUR CODE --- */
/* %mend; */

/* Question 4.9 */
PROC IMPORT datafile="/folders/myfolders/statistical-methods/data/Groningen.csv"
	out = WEEK4_Q9
	dbms = CSV;
RUN;

PROC PRINT data=WEEK4_Q9;

/* a) */
PROC REG data=WEEK4_Q9 plots=none;
	model X = /dwProb;
RUN;
/* 1st Order Autocorrelation = 0.991 */

/* b) */
DATA Residuals;
	set WEEK4_Q9;
	yhat = 0.027*(X**2) - 0.275*X;
	resid = y - yhat;
run;

/* It appears we reject the null hypothesis. */
%runs_test(data=Residuals, var=resid, alpha=0.05);


