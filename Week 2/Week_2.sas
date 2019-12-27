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
	PROC PRINT data=OUT_SCORES label noobs;
		var CLASS_ID CLASS; 
		label CLASS_ID="class"
		CLASS="group identifier";
	RUN;

	PROC PRINT data=RESULT label;
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
	input treatment$ high$ count;
	datalines;
	1 1 23
	2 1 19
	1 0 77
	2 0 81
	;
RUN;

PROC FREQ data=WEEK2_Q5;
	tables TREATMENT*HIGH /chisq;
	weight count;
	exact chisq;
RUN;
/* p-value = 0.4822 */

/* Question 2.6 */
DATA WEEK2_Q6;
	set WEEK2;
	keep SEX FIS;
RUN;

/* a) */
PROC FREQ data=WEEK2_Q6;
	tables FIS*SEX /chisq;
	exact chisq;
RUN;
/* p_hat(SEX=0 & FIS=1) = 31.09 */
/* p_hat(SEX=1 & FIS=1) = 36.09 */
/* p_hat(FIS=1) = 33.73 */

PROC TTEST data=WEEK2_Q6;
	class SEX;
	var FIS;
RUN;
/* For the H0: p_1 = p_2 vs H1: p_1 != p_2 */
/* We can assume equal so we take pooled, */
/* test statistic = -0.84 and p-value = 0.4022 */

/* b) */
%MACRO binary_hypothesis(dataset, var, class); 
	PROC MEANS data=&dataset n sum noprint;
		var &var;
		class &class;
		output out=OUT n=N sum=COUNT; 
	RUN;

	DATA OUT0;
		set OUT;
		COUNT0 = COUNT; 
		N0 = N;
		P0 = COUNT0 / N0; 
		where &class = 0; 
		keep COUNT0 N0 P0; 
	RUN;
	
	DATA OUT1;
		set OUT;
		COUNT1 = COUNT; 
		N1 = N;
		P1 = COUNT1 / N1; 
		where &class = 1; 
		keep COUNT1 N1 P1; 
	RUN;
    
    DATA OUT;
		merge OUT0 OUT1;
		P = (COUNT0 + COUNT1) / (N0 + N1);
		STAT = (P0 - P1) / sqrt(P * (1-P) * (1/N0 + 1/N1));
		CHISQ = STAT **2;
		P_VALUE = 2*min(cdf("normal", STAT, 0, 1), 1-cdf("normal", STAT, 0, 1)); 
	RUN;
	
	PROC PRINT data=OUT; 
		var STAT CHISQ P_VALUE; 
	RUN;
%MEND;

/* c) */
/* The gives the same output, */
%binary_hypothesis(WEEK2_Q6, FIS, SEX);

/* d) */
/* With the Chi-Squared test, */
/* We get Test statistic = -0.83774,(T)^2 = 0.70181  */
/* and p-value = 0.40218 */

/* Question 2.7 */
DATA WEEK2_Q7;
	set WEEK2;
	keep FIS TRT;
RUN;

/* a) */
/* Chi-Squared test */
DATA WEEK2_Q7_a;
	SET WEEK2_Q7;
	if TRT < 2;
	PROC FREQ;
		tables FIS*TRT /chisq;
		exact chisq;
RUN;

DATA WEEK2_Q7_a;
	SET WEEK2_Q7;
	if TRT NE 1;
	PROC FREQ;
		tables FIS*TRT /chisq;
		exact chisq;
RUN;

DATA WEEK2_Q7_a;
	SET WEEK2_Q7;
	if TRT > 0;
	PROC FREQ;
		tables FIS*TRT /chisq;
		exact chisq;
RUN;


/* Fisher test */
DATA WEEK2_Q7_a;
	SET WEEK2_Q7;
	if TRT < 2;
	PROC FREQ;
		tables FIS*TRT;
		exact fisher;
RUN;

DATA WEEK2_Q7_a;
	SET WEEK2_Q7;
	if TRT NE 1;
	PROC FREQ;
		tables FIS*TRT;
		exact fisher;
RUN;

DATA WEEK2_Q7_a;
	SET WEEK2_Q7;
	if TRT > 0;
	PROC FREQ;
		tables FIS*TRT;
		exact fisher;
RUN;
/* Groups | Chi Statistic | P-value | Fisher Stat | P-value */
/* 0 1	  | 2.8533		  | 0.0912  | 2.8533	  | 0.1202  */
/* 0 2	  | 4.9708		  | 0.0258  | 4.9708 	  | 0.0347  */
/* 1 2	  | 0.0458		  | 0.8306  | 0.0458 	  | 0.8563  */

/* b) */
PROC FREQ data=WEEK2_Q7;
	tables FIS*TRT /chisq;
	exact chisq;
RUN;
/* Test statistic = 5.7210 and p-value = 0.0572 */
/* We can't reject H0 */

/* Question 2.8 */
DATA WEEK2_Q8;
	input BATCH OUTPUT @@; /* Skip $ for not numeric */
	datalines;
1 102 1 104 1 102 1 97 1 99 1 101 1 103 1 98 1 96 1 97
2 99 2 97 2 99 2 100 2 99 2 96 2 99 2 98 2 97 2 98
	;
RUN;

/* a) */
/* H0: σ_2(BATCH=1) = σ_2(BATCH=2) */

/* F-Test */
/* According to F statistic we assume unequal variance */
PROC TTEST data=WEEK2_Q8; 
	class BATCH;
	var OUTPUT;
RUN;
/* Test statistic = 5.36 and p-value = 0.0199 */

/* Bartlett's test */
PROC GLM data=WEEK2_Q8; 
	class BATCH;
	model OUTPUT = BATCH;
	means BATCH / hovtest=BARTLETT;
RUN;
/* Test statistic = 5.4128 and p-value = 0.0200 */


/* Levene's test */
PROC GLM data=WEEK2_Q8; 
	class BATCH;
	model OUTPUT = BATCH;
	means BATCH / hovtest=LEVENE;
RUN;
/* Test statistic = 10.86 and p-value =  0.0040 */
/* All tests have approximately the same p-values. */
/* And we don't reject the H0 */

/* b) */

/* T-Test */
PROC TTEST data=WEEK2_Q8; 
	class BATCH;
	var OUTPUT;
RUN;
/* Test statistic = 1.73 and p-value = 0.1080 */

/* c) */
PROC NPAR1WAY data=WEEK2_Q8 wilcoxon; 
	class BATCH;
	var OUTPUT;
RUN;
/* Test statistic = 120.5 and p-value = 0.2503 */

/* d) */
%mann_whitney_u(WEEK2_Q8, BATCH, OUTPUT);
/* Test statistic (U_0) = 65.5 and p-value = 0.345 */

/* e) */
PROC NPAR1WAY data=WEEK2_Q8; 
	class BATCH;
	var OUTPUT;
	exact ks/mc;
RUN;
/* Test statistic = 0.5 and p-value = 0.0964 */
