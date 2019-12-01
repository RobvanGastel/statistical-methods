
*W2-1;
DATA WEEK2;
 SET IVF;
 WHERE TRT=0 AND PER=4;
 IMPSHIFT=IMP-10;
 label IMPSHIFT = "Neurological score - 10";
 IF AGEM>35 then IMPSO=IMP-10;
 else IMPSO=IMP;
 label IMPSO = "Neurological score AGEM>35 - 10";
 IMPU = 2*IMP-79.1363636;
 label IMPU = "Neurological score variance*4";
 RUN;

proc sgplot data=WEEK2;
   histogram IMP / fillattrs=graphdata1 transparency=0.7 binwidth=2;
   histogram IMPSO / fillattrs=graphdata2 transparency=0.5 binwidth=2;
   yaxis grid; xaxis values=(60 to 90 by 1);
run;



data CU;
	input BATCH$ MEASUREMENT@@;
	datalines;
1 102 1 104 1 102 1 97 1 99
1 101 1 103 1 98 1 96 1 97
2 99 2 97 2 99 2 100 2 99
2 96 2 99 2 98 2 97 2 98
	;
run;

proc freq data=cu;
table batch;
run;

proc ttest data=CU;
	class BATCH;
	var MEASUREMENT;
run;

data Fval;
 qU2=quantile('F',0.975,9,9);
 qL2=quantile('F',0.025,9,9);
 qU=quantile('F',0.95,9,9);
 qL=quantile('F',0.05,9,9);
 pU=1-cdf('F',5.360,9,9);
 pL=cdf('F',5.360,9,9);
 test=quantile('CHISQ',0.95,1);
run;

proc glm data=CU;
	class BATCH;
	model MEASUREMENT = BATCH;
	means BATCH / hovtest=BARTLETT;
run;


proc glm data=CU;
	class BATCH;
	model MEASUREMENT = BATCH;
	means BATCH / hovtest=Levene;
run;

proc glm data=CU;
	class BATCH;
	model MEASUREMENT = BATCH;
	means BATCH / hovtest=BF;
run;

ods output WilcoxonScores=WRS (keep= Class N SumOfScores);
proc npar1way data=CU correct=NO;
	class BATCH;
	var MEASUREMENT;
	exact wilcoxon;
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


proc npar1way data=CU correct=NO;
	class BATCH;
	var MEASUREMENT;
	exact wilcoxon /mc;
run;


proc npar1way DATA=CU;
	class BATCH;
	var MEASUREMENT;
	exact KS;
run;

