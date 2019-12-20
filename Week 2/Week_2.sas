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
PROC SQL;
	CREATE TABLE WEEK2_Q2 AS
		SELECT TRT, BW
		FROM WEEK2
		WHERE ID <= 100;
RUN;

/* a) */
/* Compare treatment (TRT) on birth weight (BW) */
/* H0: m_1 = m_2 vs H1: m_1 != m_2 */
/* 1 := (TRT = 1) etc. */

PROC SORT data=WEEK2_Q2; 
	by TRT;
	PROC MEANS mean median std n skewness kurtosis; 
		var BW; 
		by TRT; 
RUN;

DATA WEEK2_Q2_a;
	set WEEK2_Q2; 
	if TRT < 2;
	ods output WilcoxonScores=WRS_01 (keep= Class N SumOfScores);
	PROC NPAR1WAY wilcoxon; 
		class TRT;
		var BW;
	title 'TRT=0 vs TRT=1'; 
RUN;

DATA WEEK2_Q2_a; 
	set WEEK2_Q2; 
	if TRT > 0;
	ods output WilcoxonScores=WRS_12 (keep= Class N SumOfScores);
	PROC NPAR1WAY wilcoxon; 
		class TRT;
		var BW;
	title 'TRT=0 vs TRT=2'; 
RUN;

DATA WEEK2_Q2_a; 
	set WEEK2_Q2; 
	if TRT NE 1; /* NE is the same as != */
	ods output WilcoxonScores=WRS_02 (keep= Class N SumOfScores);
	PROC NPAR1WAY wilcoxon; 
		class TRT;
		var BW;
		exact wilcoxon/mc;
	title 'TRT=1 vs TRT=2'; 
RUN;

/* Table: */
/* TRT a vs b | Statistic (S) | P-value */
/* 01 			256				0.1868  */
/* 02 			519.5			0.0164  */
/* 12 			261.5			0.4891  */

/* TRT=2 has a different median birth weight */
/* as the p-value < 0.05 */

/* b) */
/* P-value doesn't change from WRS */
PROC IML;
	use WRS_01;
	read all var{Class N SumOfScores};
	close WRS;
		
	G=Num(Class); /* Num convert to numeric values */
	U=SumOfScores-N#(N+1)/2;
	P=U/prod(N);
		
	A=G||N||U||P;
	
	create MWU from A [colname={'Group' 'N' 'U' 'P'}]; 
	append from A;       
	close MWU;
QUIT;

PROC IML;
	use WRS_02;
	read all var{Class N SumOfScores};
	close WRS;
		
	G=Num(Class); /* Num convert to numeric values */
	U=SumOfScores-N#(N+1)/2;
	P=U/prod(N);
		
	A=G||N||U||P;
	
	create MWU from A [colname={'Group' 'N' 'U' 'P'}]; 
	append from A;       
	close MWU;
QUIT;

PROC IML;
	use WRS_12;
	read all var{Class N SumOfScores};
	close WRS;
		
	G=Num(Class); /* Num convert to numeric values */
	U=SumOfScores-N#(N+1)/2;
	P=U/prod(N);
		
	A=G||N||U||P;
	
	create MWU from A [colname={'Group' 'N' 'U' 'P'}]; 
	append from A;       
	close MWU;
QUIT;


/* Table: */
/* TRT a vs b | Statistic (U) | P-value */
/* 01 			277				0.1868  */
/* 02 			170.5			0.4891  */
/* 12 			538.5			0.0164  */


/* c) */
/* Looking at the probability of PROC IML results, */
/* 0 vs 1, Probability of 1: 0.6266968326 */
/* 0 vs 2, Probability of 2: 0.6886189258 */

/* d) */
/* H0: F_1 = F_2 vs H1: F_1 != F_2 */
/* F being CDF, and 1 := (TRT = 1) etc. */

DATA WEEK2_Q2_d; 
	set WEEK2_Q2; 
	if TRT < 2; /* NE is the same as != */
	PROC NPAR1WAY; 
		class TRT;
		var BW;
		/* exact ks/mc; */
	RUN;
	title 'TRT=0 vs TRT=1'; 
RUN;

DATA WEEK2_Q2_d; 
	set WEEK2_Q2; 
	if TRT NE 1; /* NE is the same as != */
	PROC NPAR1WAY; 
		class TRT;
		var BW;
		/* exact ks/mc; */
	RUN;
	title 'TRT=0 vs TRT=2'; 
RUN;

DATA WEEK2_Q2_d; 
	set WEEK2_Q2; 
	if TRT > 0; /* NE is the same as != */
	PROC NPAR1WAY; 
		class TRT;
		var BW;
		/* exact ks/mc; */
	RUN;
	title 'TRT=1 vs TRT=2'; 
RUN;

/* Table: */
/* TRT a vs b | Statistic (S) | P-value */
/* 01 			0.3620			0.1700  */
/* 02 			0.4041			0.0227  */
/* 12 			0.2642			0.6566  */

/* TRT=2 has a different birth weight  */
/* as the p-value < 0.05 */

/* e) */
DATA WEEK2_BWBOXCOX; 
	set WEEK2_Q2;
	BWPLUS2 = (0.5)*(BW**(2) -1);
RUN;

DATA WEEK2_Q2_e;
	set WEEK2_BWBOXCOX;
	if TRT < 2;
	PROC TTEST;
		class TRT;
		var BWPLUS2;
	title "TRT = 0 vs TRT = 1";
RUN;

DATA WEEK2_Q2_e;
	set WEEK2_BWBOXCOX;
	if TRT NE 1;
	PROC TTEST;
		class TRT;
		var BWPLUS2;
	title "TRT = 0 vs TRT = 2";
RUN;

DATA WEEK2_Q2_e;
	set WEEK2_BWBOXCOX;
	if TRT > 0;
	PROC TTEST;
		class TRT;
		var BWPLUS2;
	title "TRT = 1 vs TRT = 2";
RUN;

/* Table: */
/* TRT a vs b | Statistic (S) | P-value */
/* 01 			1.10			0.2771  */
/* 02 			2.15			0.0356  */
/* 12 			0.67			0.5053  */

/* We assume TRT=0 and TRT=2 have a different variance */
/* As the p-value < 0.05 */

/* f) */
/* Perform the tests again but on the whole dataset */

/* Question 2.3 */
/* a) */
%MACRO mann_whitney_u(dataset, class, var);
	ods select none;
	PROC NPAR1WAY data=&dataset wilcoxon;
		var &var;
		class &class;
		ods output WilcoxonScores=OUT_SCORES(
			rename=(SumOfScores=S));
		ods output WilcoxonTest=OUT_TEST
		ods output KruskalWallisTest=OUT_KRUS;
	RUN;
	ods select all;
	
	/* Chi square p-value */
	PROC SQL;
		CREATE TABLE P_KRUS AS
			SELECT Prob FROM OUT_KRUS;
	RUN;
	
	/* Wilcoxon p-value */
	PROC SQL;
		CREATE TABLE P_TABLE AS
			SELECT tProb2 FROM OUT_TEST;
	RUN;
	
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
		merge OUT_N OUT_S P_TABLE P_KRUS;
		P_VALUE = tProb2;
		P_KRUS = Prob;
		U0 = S0 - N0 * (N0+1)/2;
		U1 = S1 - N1 * (N1+1)/2;
		P0 = U0 / (N0*N1);
		P1 = U1 / (N0*N1);
	RUN;
	
	title "Mann Whitney U test";
	PROC PRINT data=RESULT;
		var P_VALUE P_KRUS U0 U1 P0 P1;
		label P_VALUE = "p-value Wilcoxon Test"
		P_KRUS= "p-value Kruskal-Wallis Test"
		U0="statistic (U0)" 
		U1="statistic (U1)" 
		P0="P(class0 > class1)" 
		P1="P(class0 <= class1)";
	RUN;
	title;
	
%MEND;

/* b) */
DATA WEEK2_Q3_b;
	set WEEK2_Q2; 
	if TRT < 2;
RUN;
%mann_whitney_u(WEEK2_Q3_b, TRT, BW);

DATA WEEK2_Q3_b; 
	set WEEK2_Q2; 
	if TRT > 0;
	title 'TRT=0 vs TRT=2'; 
RUN;
%mann_whitney_u(WEEK2_Q3_b, TRT, BW);

DATA WEEK2_Q3_b; 
	set WEEK2_Q2; 
	if TRT NE 1; /* NE is the same as != */
	title 'TRT=1 vs TRT=2'; 
RUN;
%mann_whitney_u(WEEK2_Q3_b, TRT, BW);

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
   input treatment $ low $ high $ count@@; 
   datalines; 
	1 77 23 100 
	2 81 19 100 
 ;
RUN;

PROC FREQ data=WEEK2_Q5;
	tables TREATMENT*LOW /chisq;
	weight HIGH;
RUN;

/* Question 2.6 */

/* a) */






