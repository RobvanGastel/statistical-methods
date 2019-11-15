LIBNAME SASDATA "/folders/myfolders/applied-stats/data";

DATA WEEK1;
	set SASDATA.IVF;
	where PER=4;
run;

/* PROC PRINT data=WEEK1 NOOBS; */

/* Question 1.1 */
PROC FREQ data=WEEK1;
	tables TRT;
Run;

/* Question 1.2 */
/* a) */
PROC MEANS DATA=WEEK1 mean var std;
	VAR AGEM;
RUN;

/* confidence intervals */
ods select BasicIntervals;
PROC UNIVARIATE data=WEEK1 cibasic(alpha=0.05);
   var AGEM;
run;
/* TODO: 95% prediction interval for a single new observation of AGEM. */

/* b) How many mothers were 40 years old or older when they became pregnant? */
/* Answer: 3 */
PROC FREQ data=WEEK1;
	tables AGEM;
	where AGEM >= 40;
Run;

/* Is AGEM really normal distributed? */
/* It seems like it is. */
ods select HISTOGRAM, QQPLOT, PPPLOT;
PROC UNIVARIATE data=WEEK1 ;
   VAR AGEM;
   HISTOGRAM AGEM/NORMAL;
   QQPLOT    AGEM/NORMAL;
   PPPLOT 	 AGEM/NORMAL;
RUN;

/* Question 1.3 */
/* a) */
ods select QUANTILES;
PROC UNIVARIATE data=WEEK1 cipctldf;
  var BW;
run;

/* p25 p50 p75 quantiles in MEANS */
PROC MEANS data=WEEK1 qrange;
	var BW;
run;

/* b) */
PROC MEANS data=WEEK1 p25 p75 qrange qmethod=p2;
	var BW;
run;




