LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

DATA bchc;
	set SASDATA.BCHC;
RUN;

PROC MEANS data=bchc;
RUN;

PROC UNIVARIATE data=bchc;
RUN;

/* Mortality has a positive skew */
PROC UNIVARIATE data=bchc normal;
	histogram mortality/normal;
RUN;

/* Box Cox transform */
/* λ ∈ {−2,−1/2,0,1/2,2} */
DATA bchc_boxcox; 
	set bchc;
	MTMINUS2 = (-1/2)*(mortality**-2 -1); 
	MTINUS1 = (-1)*(mortality**-1 -1); 
	MTMINUS12 = (-2)*(mortality**-(0.5)-1); 
	MT0 = log(mortality);
	MTPLUS12 = (2)*(mortality**(1/2) -1); 
	MTPLUS2 = (0.5)*(mortality**(2) -1);
RUN;

ods select histogram;
PROC UNIVARIATE data=bchc_boxcox; 
   	histogram Mortality /normal;
	histogram MTMINUS2 /normal; 
	histogram MTINUS1 /normal;
   	histogram MTMINUS12 /normal;
   	histogram MT0 /normal;
  	histogram MTPLUS12 /normal;
   	histogram MTPLUS2 /normal;
RUN;
/* Looks normal at MT0 */

