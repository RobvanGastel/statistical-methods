LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 6.1 */
/* a) */
DATA WEEK6_Q1;
	set SASDATA.IVF;
	where PER=4;
RUN;

PROC FREQ data=WEEK6_Q1;
	tables FIS*TRT /chisq;
	exact chisq;
RUN;
/* Comparing 3 treatments yields */
/* 1 vs 0, 1 vs 2, 0 vs 2 */
/* 𝛼 = 0.05/3 = 0.0166666667 */
/* Now the p-value = 0.0572 - 0.0166 = 0.406 */
/* So we would reject the H0 now */

/* b) */
DATA WEEK6_Q1_b;
	input treatment high count drug;
	datalines;
	1 1 12 0
	0 1 11 0
	1 0 7 0
	0 0 70 0
	1 1 10 1
	0 1 13 1
	1 0 4 1
	0 0 73 1
	;
RUN;

PROC FREQ data=WEEK6_Q1_b;
	where drug = 0;
	tables high*treatment;
	weight count;
	exact mcnem;
RUN;
/* p-value = 0.4807 */
/* The p-value is already insignificant without */
/* adjusting for multiple tests */

PROC FREQ data=WEEK6_Q1_b;
	where drug = 1;
	tables high*treatment;
	weight count;
	exact mcnem;
RUN;
/* p-value = 0.49 */
/* We still reject H0 */

/* c) */
DATA WEEK6_Q1_c;
	input supplement$ iq@@; 
	datalines;
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
	;
RUN;

PROC MIXED data=WEEK6_Q1_c method=TYPE3 cl; 
	class supplement;
	model iq = supplement /solution cl;
	lsmeans supplement/diff=control("suppC") cl adjust=dunnett;
RUN;
/* H0: μA = μC, p-value = 0.0016 */
/* H0: μB = μC, p-value = 0.1295 */

/* Question 6.2 */
DATA WEEK6_Q2;
	set SASDATA.IVF;
RUN;

/* a) */
/* y_ij = mu + a_i + b_j + e_ij */

/* b) */
PROC MIXED data=WEEK6_Q2 method=TYPE3 cl; 
	class TRT ID;
	model IMP = TRT /solution cl;
	random ID;
	lsmeans TRT;
RUN;
/* as the covariance param of σ_G is negative we set it to 0 */

PROC IML;
	/* Estimate COV of random effect */
	S_g = 0; 
	/* Estimate COV of residuals */
	S_res = 56.7474;
	ICC = S_g / (S_g + S_res);
	
	A = ICC;
	create ICC from A [colname={'ICC'}]; 
	append from A;
	close ICC;
RUN;
/* as we can see the ICC = 0 */
/* This tells us that the variance between subjects  */
/* is negligible compared to the variance within subjects. */

/* c) */
/* H0: mu_1 = mu_2 = ... = mu_4 = 0 */
/* H1: mu_1 = mu_2 = .. = mu_4 != 0 */
/* F Value = 5.47, gives p-value = 0.045 */

/* d) */
PROC MIXED data=WEEK6_Q2 method=TYPE3 cl; 
	class TRT ID;
	model IMP = TRT /solution cl;
	random ID;
	lsmeans TRT/diff adjust=tukey cl;
RUN;
/* Comparison | p-value adj */
/* 0 vs 1	  | 0.9999 */
/* 0 vs 2     | 0.0033 */
/* 1 vs 2     | 0.0100 */

/* e) */
/* Scheffe's method is more convservative if this method's CI */
/* Doesn't contain 0 than Tukey's CI surely doesn't. */

/* f) */
/* Scheffe's method has a wider CI so we're not certain if we */
/* should deny the test using Scheffe's method. */

/* g) */
/* We can use the Kruskal-Wallis test */
/* Careful! Apperantly Kruskal-Wallis doesn't apply as it */
/* doesn't take random effects into account. */

/* Question 6.3 */
DATA WEEK6_Q3;
	set SASDATA.RCT;
	where time = 1;
RUN;

/* You can assume that the hemoglobin values (RESP)  */
/* can be modeled with an ANOVA model. We will  */
/* investigate the effect of the medical center CENTER  */
/* in the RCT dataset taking into account the effect of TRT. */

/* a) */
/* Y_ij = mu + a_i + b_j + e_ij */
/* b being for the centers */
/* a being the trt */

/* b) */
PROC MIXED data=WEEK6_Q3 method=TYPE3 cl; 
	class CENTER TRT;
	model RESP = TRT CENTER /solution cl;
	lsmeans TRT/diff=CONTROL adjust=TUKEY cl;
RUN;

/* c) */
/* We observe F value = 6.31 and p-value < 0.0001 */

/* d) */
PROC MIXED data=WEEK6_Q3 method=TYPE3 cl; 
	class CENTER;
	model RESP = CENTER /solution cl;
RUN;

/* e) */
PROC MIXED data=WEEK6_Q3 method=TYPE3 cl; 
	class CENTER TRT;
	model RESP = TRT CENTER /solution cl;
	lsmeans TRT/diff=control adjust=tukey cl;
	lsmeans CENTER/diff=control adjust=tukey cl;
RUN;

/* f) */
PROC MIXED data=WEEK6_Q3 method=TYPE3 cl; 
	class CENTER TRT;
	model RESP = TRT CENTER /solution cl;
	lsmeans TRT/diff=control adjust=dunnet cl;
RUN;

/* g) */
/* Not Kruskal-Wallis again */

/* Question 6.4 */
/* a) A Bonferroni or Sidak experiment 
	wise correction could be applied */
/* b) Least square difference */
/* c) Tukey’s studentized range test 
	> Bonferroni Experiment-wise correction 
	> F-test.
/* d) Experiment wise Bonferroni-Holm correction */
/* e) Dunnett’s many to one adjustment.
/* f) Scheffe's method, is most conservative */




