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
   class TRT;
   var BW;
   exact wilcoxon / mc;
run;

PROC IML;
use WRS;
read all var{N SumOfScores}; 
close WRS;

G={1 , 2}; 
U=SumOfScores-N#(N+1)/2; 
P=U/prod(N);

A=G||N||U||P;
create MWU from A [colname={'Group' 'N' 'U' 'P'}]; 
append from A;
close MWU;
quit;


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

