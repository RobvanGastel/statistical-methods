LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

DATA WEEK2;
	set SASDATA.IVF;
	where PER=4;
	drop IMP PER AGE;
RUN;

/* Question 2.1 */
/* Assume mother's age (AGEM) is normally distributed */

/* a) */
/* H0: σ_2(FIS=0) = σ_2(FIS=1) */

/* T-Test and F-Test */
PROC TTEST data=WEEK2; 
	class FIS;
	var AGEM; 
RUN;
/* test statistic = 1.03 and p-value = 0.8490 */

/* Bartlett's test */
PROC GLM data=WEEK2; 
	class FIS;
	model AGEM = FIS;
	means FIS / hovtest=BARTLETT;
RUN;
/* test statistic = 0.0288 and p-value = 0.8653 */

/* Levene's test */
PROC GLM data=WEEK2; 
	class FIS;
	model AGEM = FIS;
	means FIS / hovtest=LEVENE;
RUN;
/* test statistic = 0.0200 and p-value = 0.8813 */
/* All tests have approximately the same p-values. */
/* And we don't reject the H0 */

/* b) */
/* H0: μ_2(F, FIS=0) = μ_2(F, FIS=1)*/

/* T-Test */
PROC TTEST data=WEEK2 alpha=0.05 test=diff ci=equal; 
	class FIS;
	var AGEM;
RUN;
/* Because of HOV tests we assume variance is equal */
/* Thus we used pooled T-Test. */
/* Test Statistic = −0.83 and p-value = 0.4065. */

/* Question 2.2 */
DATA WEEK2_Q2;
	set WEEK2;
	where ID <= 100 and TRT = 0 or TRT = 1;
RUN;

PROC FREQ data=WEEK2_Q2;
	tables TRT;
RUN;

/* or SQL */
PROC SQL;
	CREATE TABLE WEEK2_Q2 AS
	SELECT TRT, BW
	FROM WEEK2
	WHERE ID <= 100
	AND TRT = 0 OR TRT = 2;
RUN;

/* a) */
/* Compare treatment (TRT) on birth weight (BW) */
/* H0: m_1 = m_2 vs H1: m_1 != m_2 */
/* TODO */
ods output WilcoxonScores=WRS (keep= Class N SumOfScores); 
PROC NPAR1WAY data=WEEK2_Q2 correct=no;
   class TRT;
   var BW;
   exact wilcoxon/mc;
RUN;


/* b) */
/* TODO */
PROC IML;
	use WRS;
	read all var{N SumOfScores};
	close WRS;
		
	G={0, 1};
	U=SumOfScores-N#(N+1)/2;
	P=U/prod(N);
		
	A=G||N||U||P;
	
	create MWU from A [colname={'Group' 'N' 'U' 'P'}]; 
	append from A;       
	close MWU;
QUIT;

/* c) */

/* d) */

/* e) */

/* f) */


/* Question 2.3 */
/* a) */
%MACRO mann_whitney_u(dataset, class, var);
	ods select none;
	PROC NPAR1WAY data=&dataset;
		var &var;
		class &class;
		exact wilcoxon/mc;
		ods output WilcoxonScores=OUT_SCORES(
			rename=(SumOfScores=S));
		ods output WilcoxonTest=OUT_TEST(
			rename=(cValue1=P_VALUE) where =
			( Name1 =" P2_WIL "));
	RUN;
	ods select all;
	
	DATA OUT_SCORES;
		set OUT_SCORES; 
		CLASS_ID = _N_ - 1;
	RUN;
	
	PROC TRANSPOSE data=OUT_SCORES
			out=OUT_N(drop=_NAME_) prefix=N;
		id CLASS_ID;
		var N;
	RUN;

	PROC TRANSPOSE data=OUT_SCORES
			out=OUT_S(drop=_NAME_ _LABEL_) PREFIX=S;
		id CLASS_ID;
		var S;
	RUN;
	
	DATA RESULT;
		merge OUT_N OUT_S OUT_TEST(keep=P_VALUE); 
		U0 = S0 - N0 * (N0+1)/2;
		U1 = S1 - N1 * (N1+1)/2;
		P0 = U0 / (N0*N1);
		P1 = U1 / (N0*N1);
	RUN;

	title "Mann Whitney U test";
	PROC PRINT data=OUT_SCORES label noobs;
		var CLASS_ID CLASS; 
		label CLASS_ID="class"
		CLASS="group identifier";
	RUN;
	title;
	
	PROC PRINT data=RESULT;
		var P_VALUE U0 U1 P0 P1;
		label P_VALUE="p-value"
		U0="statistic (U0)" 
		U1="statistic (U1)" 
		P0="P(class0 > class1)" 
		P1="P(class0 <= class1)";
	RUN;
%MEND;

%mann_whitney_u(WEEK2_Q2, TRT, BW);

/* b) */
/* TODO */

/* Question 2.4 */
DATA WEEK2_Q4;
	SET WEEK2;
	keep GA BW;
RUN;

/* a) */
DATA BW_HEAVY;
	set WEEK2_Q4;
	if BW > 4000 then heavy=1;
	else heavy=0;
RUN;

/* b) */
/* H0: μ(BW>4000) = μ(BW<4000) */
PROC TTEST data=BW_HEAVY;
	class heavy;
	var GA;
RUN;
/* F-Test */
/* Test statistic = 5.30 and p-value = .0001 */
/* Thus we can not assume equal variance, */
/* Test statistic = -9.07 and p-value = .0001 */
/* Therefore we reject H0 */

/* c) */
DATA GA_LATE;
	SET WEEK2_Q4;
	if GA > 41 then late=1;
	else late=0;
RUN;

/* d) */
/* H0: μ(GA>41) = μ(GA<41) */
PROC TTEST data=GA_LATE;
	class late;
	var BW;
RUN;
/* F-test */
/* Test statistic = 1.77 and p-value = 0.0407 */
/* Thus we can not assume equal variance, */
/* Test statistic = -6,55 and p-value = 0.0001 */
/* Therefore we reject the H0 */

/* e) */
/* As we're using the mean with a relatively large */
/* sample size n=253. Based on the CLT the results */
/* Are probably reliable. */

/* Question 2.5 */
/* It's a reasonable approximation as the count */
/* in each cell > 5. */
/* H0: p_1 = p_2 */
DATA WEEK2_Q5; 
   input treatment $ low $ high $ level@@; 
   datalines; 
	1 77 23 2
	2 81 19 1
 ;
RUN;

PROC FREQ data=WEEK2_Q5;
	tables TREATMENT*LEVEL;
RUN;


Data CU;
Input BATCH$ MEASUREMENT@@; datalines;
1 102 1 104 1 102 1 97 1 99
1 101 1 103 1 98 1 96 1 97
2 99 2 97 2 99 2 100 2 99
2 96 2 99 2 98 2 97 2 98
;
run;

Data CU;
	Set CU;
	If MEASUREMENT<97 OR MEASUREMENT>103 then OSPECLIM=1; Else OSPECLIM=0;
run;

PROC PRINT data=CU;
RUN;

Proc freq data=CU;
Table BATCH*OSPECLIM / chisq; Exact chisq;
RUN;



