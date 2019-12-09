LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

DATA WEEK2;
	set SASDATA.IVF;
	where PER=4;
	drop IMP PER AGE;
run;

/* Question 2.1 */
/* Assume mother's age (AGEM) is normally distributed */
/* FIS */

/* a) */
/* null hypothesis: σ_2(FIS=0) = σ_2(FIS=1)*/

/* T-Test and F-Test */
/* https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_ttest_a0000000113.htm */
PROC TTEST data=WEEK2; 
	class FIS;
	var AGEM; 
run;
/* test statistic:  1.03  */
/* p-value: 		0.8490 */

/* Bartlett's test */
PROC GLM data=WEEK2; 
	class FIS;
	model AGEM = FIS;
	means FIS / hovtest=BARTLETT;
run;
/* test statistic: 0.0288 */
/* p-value: 	   0.8653 */

/* Levene's test */
PROC GLM data=WEEK2; 
	class FIS;
	model AGEM = FIS;
	means FIS / hovtest=levene;
run;
/* test statistic: 0.0200  */
/* p-value: 	   0.8813 */


/* All tests have approximately the same values
   as p-values. */

/* b) */

/* null hypothesis: μ_2(F, FIS=0) = μ_2(F, FIS=1)*/
/* T-Test */
/* ods graphics off; */
PROC TTEST data=WEEK2 alpha=0.05 test=diff ci=equal; 
	class FIS;
	var AGEM;
run;
/* Since we didn’t reject the null hypothesis in (a),  */
/* we’ll assume that the variances in both groups are  */
/* equal. */
/* Thus use the Pooled t- test. For this,  */
/* Test Statistic = −0.83, p-value = 0.4065.  */
/* Therefore, we don’t reject the null-hypothesis. */

/* Question 2.2 */
DATA WEEK2_2;
	set WEEK2;
	where ID <= 100;
run;

PROC FREQ DATA=WEEK2_2;
	tables TRT;
run;

/* Compare treatment (TRT) on birth weight (BW) */
/* a) */
/* TODO */
ods output WilcoxonScores=WRS (keep= Class N SumOfScores); 
proc npar1way data=WEEK2_2 correct=NO;
   class TRT;
   var BW;
   exact wilcoxon / mc;
run;

/* b-f */
/* TODO */

/* Question 2.3 */

ods output WilcoxonScores=WRS (keep= Class N SumOfScores); 
proc npar1way data=WEEK2_2 correct=NO;
   class FIS;
   var AGEM;
   exact wilcoxon / mc;
run;

%macro mann_whitney_u(dataset, class, var); 

ods select none;
proc npar1way data=&dataset;
	var &var;
	class &class;
	exact wilcoxon / mc;
	ods output WilcoxonScores=OUT_SCORES; 
	ods output WilcoxonTest=OUT_TEST;
run;

PROC IML;
use WRS;
read all var{N SumOfScores}; close WRS;
G={1 , 2}; U=SumOfScores-N#(N+1)/2; P=U/prod(N);
 A=G||N||U||P;
create MWU from A [colname={'Group' 'N' 'U' 'P'}]; append from A;
close MWU;
quit;

ods select all;

%mend;

%mann_whitney_u(WEEK2_2, FIS, AGEM);

/* Question 2.4 */
DATA WEEK2_4;
	SET WEEK2;
	keep GA BW;
run;

/* a) */
DATA BW_HEAVY;
	set WEEK2_4;
	if BW > 4000 then heavy=1;
	else heavy=0;
run;

/* b) */
PROC TTEST data=BW_HEAVY;
	class heavy;
	var GA;
run;

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
run;

/* d) */
PROC TTEST data=GA_LATE;
	class late;
	var BW;
run;

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


