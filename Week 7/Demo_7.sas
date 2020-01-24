
LIBNAME SAS "/folders/myfolders";

DATA IVF_DATASET;
	SET SAS.IVF;
RUN;

*B;
PROC MIXED DATA=IVF_DATASET METHOD=TYPE3;
CLASS ID TRT;
MODEL IMP = TRT / solution;
Random ID;
LSMEANS TRT/DIFF ADJUST=Tukey;
run;

