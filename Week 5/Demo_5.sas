*W4-4;

data COAG ;
input Patient C K@@;
datalines ;
1 120 132 8 145 133 15 117 123
2 114 116 9 120 123 16 125 108
3 129 135 10 129 116 17 136 131
4 128 115 11 126 127 18 151 119
5 155 134 12 136 140 19 130 129
6 105 56 13 135 140 20 136 124 
7 114 114 14 125 114 21 113 112
;
run;

proc sort data=COAG;
by Patient;
run;

proc transpose data=COAG out=W44(RENAME= (COL1=RESP _NAME_=TRT));
by Patient;
run;

*A/D/F;
ods output SolutionR=EBLUP;
proc mixed data=W44 method=type3;
class Patient TRT;
model RESP = /solution outp=PRED;
random Patient /solution;
run;

*B;

*sigma_G^2 = 122.21;
*sigma_E^2 = 123.81;


proc iml;
ICC=122.21/(122.21+123.81);
print(ICC);
run;

*ICC=0.4967482;

*C;
*EBLUP_8=9.6088, so mean patient 8 = 9.6088+124.52;

*D;
*This is the order the data has been collected;
proc sort data=EBLUP;
by Patient;
run;

*E;
%RUNSCUC(data=EBLUP,var=Estimate,alpha=0.05);

*runs=12;
*p-value= 0.8141239;
*Thus randomness (independence) of random effects (alpha_i's) cannot be rejected;

*F;
*This is the order the data has been collected;
proc sort data=PRED;
by Patient TRT;
run;


%RUNSCUC(data=PRED,var=RESID,alpha=0.05);

*runs=21;
*p-value=0.7547058;
*Thus randomness (independence) of residuals (e_ij's) cannot be rejected;
