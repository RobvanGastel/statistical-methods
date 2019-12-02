PROC IML;

/* Confidence interval */
use WEEK1;
read all var{agem}; close WEEK1;
alpha=0.05;
Ybar=mean(agem); 
s=var(agem); 
n=nrow(agem); 
qT=quantile('t',alpha/ 2,n-1); 
UCL=Ybar-qT*sqrt(s/n); 
LCL=Ybar+qT*sqrt(s/n);

A=Ybar||LCL||UCL;
create DATA from A[colname={'mean' 'LCL' 'UCL'}]; append from A;
close DATA;
quit;