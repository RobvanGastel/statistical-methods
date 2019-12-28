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

/* a) */
/* Fixed model, */
/* 𝒀_ij = 𝝁 + 𝜶_i + 𝒆_ij */

/* b) */
/* H0:𝛼_1=𝛼_2=⋯=𝛼_m=0 vs H1:∃𝑖,𝛼_𝑖≠0 */
PROC MIXED data=WEEK4_Q1 method=TYPE3 cl; 
	class supplement;
	model iq = supplement /solution cl; 
RUN;
/* For 0 the estimate is -4.1550 and p-value = 0.0747 */
/* Thus, we can't reject H0. */

/* c) */
/* TODO: CI don't match ANOVA */
PROC IML;
	Y_bar = 107.68;
	std = 1.7120;
	alpha = 0.05;
	n = 40;
	qt = quantile("Normal", alpha/2, n);

	LPL = Y_bar - qt * (std/sqrt(n));
	UPL = Y_bar + qt * (std/sqrt(n));
	
	A = LPL||UPL;
	create CI_c from A [colname={'UPL','LPL'}]; 
	append from A;
	close CI_c;
QUIT;

/* d) */
/* TODO */
PROC TTEST data=WEEK4_Q1;
	class supplement;
	var iq; 
RUN;
/* Test statistic = -1.83 and p-value = 0.0747 */
/* So this is the same as the ANOVA */

/* e) */
/* TODO */

/* f) */
/* TODO */

/* g) */
/* TODO */
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

PROC MIXED data=WEEK4_Q1_g method=TYPE3 cl; 
	class supplement;
	model iq = supplement /solution cl; 
RUN;

/* h) */
/* TODO */

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

/* a) */
PROC MIXED data=COAG method=TYPE3 cl;
	class Patient;
	model K = Patient /solution cl;
RUN;





