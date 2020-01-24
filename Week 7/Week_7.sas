LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 7.1 */
DATA WEEK7_Q1;
	set SASDATA.RCT;
	where ID <=75;
RUN;

/* a) */
PROC MIXED data=WEEK7_Q1 method=TYPE3 cl;
	class TRT ID;
	model RESP=TRT /solution cl ddfm=satterthwaite;
	random ID;
RUN;

/* b) */
PROC MIXED data=WEEK7_Q1 method=TYPE3 cl;
	class TRT;
	model RESP=TRT /solution cl ddfm=satterthwaite;
RUN;

/* TODO: Compare */
/* c) */
PROC MIXED data=WEEK7_Q1 method=TYPE3 cl;
	class TRT TIME ID;
	model RESP=TRT /solution cl ddfm=satterthwaite;
	random ID TIME;
RUN;

/* d) */
/* TODO: Answer is close but still off. */
PROC FREQ data=WEEK7_Q1;
	tables TRT*RESP*TIME /cmh2 scores=rank noprint;
RUN;

/* e) */
/* TODO */

/* Question 7.2 */
DATA IVF;
	set SASDATA.IVF;
RUN;

/* You can assume that the infant motor profile (IMP), */
/* that is observed at three different dates, can be */
/* modeled with a mixed effects ANOVA model. We will */
/* consider the fixed effects treatment (TRT) and period */
/* (PER), their interaction, and the random effect ID. */

/* a) */
/* fixed effects TRT, PER */
/* random effect ID */
/* effect trying to model is IMP */
/* Y_ijk = mu + a_i + b_j + c_k + e_ijk */

/* b) */
/* The traditional hypothesis for a fixed effects vars, */
/* H0: 𝜶𝟏 + ⋯ + 𝜶𝒎 = 𝟎  vs. H1: 𝜶𝟏 + ⋯ + 𝜶𝒎 != 𝟎 */
/* TODO: p-value is off by 0.0001 */
ods output SolutionR = solR;
PROC MIXED data=IVF method=TYPE3 cl;
	class TRT PER ID;
	model IMP = TRT PER /solution cl ddfm=SAT outpm=RM outp=RC;
	random ID(TRT) /solution;
RUN;

/* c) */
/* TODO */

/* d) */
/* TODO: Also this estimate is slightly off */
PROC IML;
	/* Estimate COV of random effect */
	S_g = 1.5477; 
	/* Estimate COV of residuals */
	S_res = 7.1145;
	ICC = S_g / (S_g + S_res);
	
	A = ICC;
	create ICC from A [colname={'ICC'}]; 
	append from A;
	close ICC;
RUN;

DATA _;
	set solR;
	where ID = 9 or ID = 11 or ID = 33;
RUN;

/* e) */
PROC UNIVARIATE data=RC normaltest;
	var resid;
	probplot resid/normal(mu=est sigma=est);
	histogram resid/normal;
RUN;

PROC GLM data=RC;
	class TRT PER;
	model resid = TRT*PER;
	means TRT*PER/ hovtest=Bartlett;
RUN;

/* f) */
PROC MIXED data=IVF method=TYPE3 cl;
	class TRT PER ID;
	model IMP = TRT PER /solution cl ddfm=SAT outpm=RM outp=RC;
	random ID(TRT) /solution;
	lsmeans TRT/diff=control adjust=tukey cl;
RUN;

/* g) */
/* TODO */

/* Question 7.3 */
DATA RCT;
	set SASDATA.RCT;
RUN;

/* a) */
/* TODO */

/* b) */
/* TODO */
PROC MIXED data=RCT method=TYPE3 cl;
	class TRT TIME CENTER ID;
	model RESP = TRT TIME CENTER TRT*TIME TRT*CENTER TIME*CENTER TRT*TIME*CENTER /solution cl ddfm=SAT outpm=RM outp=RC;
	random ID(TRT*CENTER) /solution;
	lsmeans TRT/diff=control adjust=tukey cl;
RUN;

/* c) */
























