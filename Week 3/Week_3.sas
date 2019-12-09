LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";


%Macro SIM_Gum(alpha=, nsim=, seed=); 
proc iml;
call streaminit(&seed);
alpha=&alpha;
do i=1 to &nsim by 1; 
U1=rand('Uniform'); 
U2=rand('Uniform');

start Func(x) global(U1,U2,alpha); 
return(Exp(-((-Log(x))**alpha + (-Log(U1))**alpha)**(1/alpha)) *((-Log(x))**alpha + (-Log(U1))**alpha)**(-1 + 1/alpha)* ((-Log(U1))**(alpha -1))/U1-U2);
finish;

intervals = {0.00001 1};
U2C = froot("Func", intervals);
X=X//U1; Y=Y//U2C; YI=YI//U2; 
end;

Total=X||Y||YI;
create GumC from Total [colname={'X','Y','YI'}]; append from Total;
close GumC;
quit;
%mend SIM_Gum;

%SIM_Gum(nsim=1000, alpha=5, seed=12345);

/* Question 3.1 */
/* a) */
proc corr data=GUMC kendall; 
	var X Y;
run;

proc corr data=GUMC spearman; 
	var X Y;
run;