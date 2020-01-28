LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 7.1 */
DATA WEEK7_Q1;
	set SASDATA.RCT;
	where ID <=75;
RUN;

/* a) */
PROC MIXED data=WEEK7_Q1 method=TYPE3 cl;
	class TRT ID;
	model RESP=TRT /solution cl ddfm=SAT;
	random ID;
RUN;

/* b) */
PROC MIXED data=WEEK7_Q1 method=TYPE3 cl;
	class TRT;
	model RESP=TRT /solution cl ddfm=SAT;
RUN;

/* c) */
PROC MIXED data=WEEK7_Q1 method=TYPE3 cl;
	class TRT TIME ID;
	model RESP=TRT /solution cl ddfm=SAT;
	random ID TIME;
RUN;

/* d) */
PROC FREQ data=WEEK7_Q1;
	tables TRT*ID*TIME*RESP /cmh2 scores=rank noprint;
RUN;

/* e) */
/* Friedman as the assumptions of the ANOVA model are not */
/* valid. */

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
/* We reject the H0 as p-value = 0.0031 */
ods output solutionR = solR;
PROC MIXED data=IVF method=TYPE3;
	class TRT PER ID;
	model IMP = TRT PER TRT*PER /solution outp=RC outpm=RM;
	random ID(TRT) /solution;
RUN;

/* c) */
/* We can't say age has an influence on the perforamence */
/* of the treatment. As the p-value = 0.1205 */

/* d) */
PROC IML;
	/* Estimate COV of random effect */
	S_g = 1.5449; 
	/* Estimate COV of residuals */
	S_res = 7.0802;
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
	model resid = TRT*PER*TRT*PER;
	means TRT*PER*TRT*PER/ hovtest=Bartlett;
RUN;

/* f) */
PROC MIXED data=IVF method=TYPE3 cl;
	class TRT PER ID;
	model IMP = TRT PER TRT*PER /solution ddfm=SAT;
	random ID(TRT) /solution;
	lsmeans TRT/diff=control('0') adjust=dunnett cl;
RUN;

/* g) */
PROC MIXED data=IVF method=TYPE3 cl;
	class TRT PER ID;
	model IMP = TRT PER TRT*PER /solution ddfm=SAT;
	random ID(TRT) /solution;
	lsmeans TRT/diff=control('0') adjust=dunnett cl;
RUN;

/* Question 7.3 */
DATA RCT;
	set SASDATA.RCT;
RUN;

/* a) */
/* - */

/* b) */
ods output SolutionR = SolR;
PROC MIXED data=RCT method=TYPE3 cl;
	class TRT TIME CENTER ID;
	model RESP = TRT TIME CENTER TRT*TIME TRT*CENTER TIME*CENTER TRT*TIME*CENTER /solution cl ddfm=SAT outpm=RM outp=RC;
	random ID(TRT*CENTER) /solution;
RUN;

/* c) */
/* - */

/* d) */
PROC IML;
	/* Estimate COV of random effect */
	S_g = 0.3299; 
	/* Estimate COV of residuals */
	S_res = 0.3882;
	ICC = S_g / (S_g + S_res);
	
	A = ICC;
	create ICC from A [colname={'ICC'}]; 
	append from A;
	close ICC;
RUN;

DATA _;
	set SolR;
	where ID > -1 and ID < 4;
RUN;

/* e) */
ods output SolutionR = SolR;
PROC MIXED data=RCT method=TYPE3 cl;
	class TRT TIME CENTER ID;
	model RESP = TRT TIME CENTER TRT*TIME TRT*CENTER TIME*CENTER TRT*TIME*CENTER /solution cl ddfm=SAT outpm=RM outp=RC;
	random ID(TRT*CENTER);
	lsmeans TRT /diff=control('0') adjust=dunnett cl;
RUN;

/* f) */
/* TODO */

/* g) */
data RCT_DATASET_NEW;
	set RCT;
	where CENTER =1;
	if TIME < 4 then SENSOR = 1; else SENSOR = 2;
	keep ID RESP SENSOR;
run;

PROC MIXED data=RCT_DATASET_NEW method=TYPE3 cl;
	class ID SENSOR;
	model RESP = ID /solution ddfm=SAT outpm=RM outp=RC;
	random SENSOR(ID);
RUN;

/* h) */
/* TODO */

/* i) */
/* TODO */





