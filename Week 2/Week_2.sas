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
/* T-Test and F-Test */
PROC TTEST data=WEEK2; 
	class FIS;
	var AGEM; 
run;

/* Bartlett's test */
PROC GLM data=WEEK2; 
	class FIS;
	model AGEM = FIS;
	means FIS / hovtest=BARTLETT;
run;

/* Levene's test */
PROC GLM data=WEEK2; 
	class FIS;
	model AGEM = FIS;
	means FIS / hovtest=levene;
run;
