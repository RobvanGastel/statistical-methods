LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 7.1 */
DATA WEEK7_Q1;
	set SASDATA.RCT;
	where ID < 76;
RUN;

/* a) */
PROC MIXED data=WEEK7_Q1 method=TYPE3 cl; 
	class TRT ID;
	model RESP = TRT /solution cl ddfm=satterthwaite;
	random ID;
RUN;

/* b) */
PROC MIXED data=WEEK7_Q1 method=TYPE3 cl; 
	class TRT;
	model RESP = TRT /solution cl ddfm=satterthwaite;
RUN;
/* TODO: Compare */

/* c) */
PROC MIXED data=WEEK7_Q1 method=TYPE3 cl; 
	class TRT TIME ID;
	model RESP = TRT /solution cl ddfm=satterthwaite;
	random ID TIME;
RUN;

/* d) */
/* TODO */
PROC FREQ data=WEEK7_Q1;
	tables TIME*TRT*RESP /cmh2 scores=rank noprint;
RUN;

/* e) */
/* TODO */

/* Question 7.2 */









