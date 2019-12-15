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
	where ID <= 100 and TRT=0 OR TRT=1;
RUN;

PROC FREQ data=WEEK2_Q2;
	tables TRT;
RUN;

/* a) */
/* Compare treatment (TRT) on birth weight (BW) */
/* H0: m_1 = m_2 vs H1: m_1 != m_2 */
ods output WilcoxonScores=WRS (keep= Class N SumOfScores); 
PROC NPAR1WAY data=WEEK2_Q2 correct=NO;
   class TRT;
   var BW;
   exact wilcoxon/mc;
RUN;

/* b) */


/* c) */


/* d) */

/* e) */

/* f) */


/* Question 2.3 */

ods output WilcoxonScores=WRS (keep= Class N SumOfScores); 
PROC NPAR1WAY data=WEEK2_Q2 correct=NO;
   class FIS;
   var AGEM;
   exact wilcoxon / mc;
RUN;

%macro mann_whitney_u(dataset, class, var); 

ods select none;
PROC NPAR1WAY data=&dataset;
	var &var;
	class &class;
	exact wilcoxon / mc;
	ods output WilcoxonScores=OUT_SCORES; 
	ods output WilcoxonTest=OUT_TEST;
RUN;

PROC IML;
use WRS;
read all var{N SumOfScores}; close WRS;
G={1 , 2}; U=SumOfScores-N#(N+1)/2; P=U/prod(N);
 A=G||N||U||P;
create MWU from A [colname={'Group' 'N' 'U' 'P'}]; append from A;
close MWU;
QUIT;

ods select all;

%mend;

%mann_whitney_u(WEEK2_2, FIS, AGEM);

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
PROC TTEST data=BW_HEAVY;
	class heavy;
	var GA;
RUN;

/* H0: μ (BW>4000) =  μ (BW<4000) */
/* F-Test  */
/* Test statistic: 5.30 with p-value: .0001 */
/* Thus we can not assume equal variance */
/* Unequal, Test statistic: -9.07 with p-value:	.0001 */
/* Therefore we reject H0 */

/* c) */
DATA GA_LATE;
	SET WEEK2_4;
	if GA > 41 then late=1;
	else late=0;
RUN;

/* d) */
PROC TTEST data=GA_LATE;
	class late;
	var BW;
RUN;

/* H0: μ (GA>41) =  μ (GA<41) */
/* F-test */
/* With T statistic: 1.77 and p-value: 0.0407 */
/* So we reject the assumption of equal variance */
/* Under unequal variance assumption  */
/* T statistic: -6,55 with p-value: 0.0001 */
/* We reject the H0 */

/* e) */
/* TODO */

/* Question 2.5 */
/* TODO */

/* Question 2.7 */
/* a) */
DATA WEEK2_7;
	set WEEK2;
	keep FIS TRT;
RUN;




