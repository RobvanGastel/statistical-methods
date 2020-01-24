***W7-1***;

*Friedman test;
PROC FREQ DATA=MEANS;
TABLES SUBJECT*ONCOLOGIST*VOLUME / CMH2 SCORES=RANK NOPRINT;
RUN;

*Ex 6.2;
LIBNAME SAS "/folders/myfolders";

DATA IVF_DATASET;
	SET SAS.IVF;
RUN;

PROC MIXED DATA=IVF_DATASET METHOD=TYPE3;
	CLASS ID TRT;
	MODEL IMP = TRT /SOLUTION  OUTP=RES;
	RANDOM ID;
RUN;

PROC MIXED DATA=IVF_DATASET METHOD=TYPE3;
	CLASS ID TRT;
	MODEL IMP = TRT /SOLUTION  OUTP=RES;
	RANDOM ID(TRT);
RUN;

*New nested balanced dataset;

DATA IVF2;
SET IVF_DATASET; 
if cmiss(of _ALL_) then delete;
run;

proc sort data=IVF2;
by ID trt;
run;

proc transpose data=IVF2 out=IVF_wide prefix=per; /* Convert to wide format */
	by ID trt;
	id per;
	var IMP;
run;

proc sort data=IVF_wide;
by TRT;
run;

DATA IVF_wide;
SET IVF_wide; 
if cmiss(of _ALL_) then delete;
run;

data IVFB;
  set IVF_WIDE;
  count + 1;
  by trt;
  if first.trt then count = 1;
run;

proc sort data=IVFB;
by ID TRT COUNT;
run;

proc transpose data=IVFB out=IVF_long; /* Convert to long format */
	by ID TRT COUNT;
run;

data IVF_long;
  set IVF_long;
  per=input(substr(_name_, 4), 5.);
  drop _name_;
run; 

data IVF_balanced;
set IVF_long;
where count<51;
run;

PROC MIXED DATA=IVF_balanced METHOD=TYPE3;
	CLASS ID TRT;
	MODEL IMP = TRT /SOLUTION  OUTP=RES;
	RANDOM ID;
RUN;

*Why it�s necessary to model nesting?;
PROC MIXED DATA=IVF_balanced METHOD=TYPE3;
	CLASS ID TRT;
	MODEL IMP = TRT /SOLUTION  OUTP=RES;
	RANDOM ID;
RUN;

*Interaction model;
PROC SGPLOT DATA=MEANS;
	SERIES X=ONCOLOGIST Y=VOLUME/ GROUP=SUBJECT;
	SCATTER X=ONCOLOGIST Y=VOLUME/ GROUP=SUBJECT;
RUN;

PROC SORT DATA=INPUT;
	BY SUBJECT ONCOLOGIST REPEAT;
RUN;

PROC MIXED DATA=INPUT METHOD=TYPE3 COVTEST CL;
	CLASS SUBJECT ONCOLOGIST;
	MODEL VOLUME=ONCOLOGIST/SOLUTION DDFM=SAT;
	RANDOM SUBJECT SUBJECT*ONCOLOGIST;
RUN;

***W7-2***;
DATA ANALYSIS;
SET SAS.example_invivo_data;
run;

PROC MIXED DATA=ANALYSIS METHOD=TYPE3 CL;
  CLASS DOSES RUNS CAGES;
  MODEL LOGRESP = DOSES/SOLUTION DDFM=SAT CL;
  RANDOM RUNS CAGES(RUNS) DOSES*RUNS DOSES*CAGES(RUNS);
RUN;



