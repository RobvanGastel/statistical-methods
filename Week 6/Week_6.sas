LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 6.1 */
DATA WEEK6_Q1;
	set SASDATA.IVF;
	where PER=4;
RUN;

/* a) */