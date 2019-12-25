LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 3.1 */
%Macro SIM_Gum(alpha=, nsim=, seed=); 
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		do i=1 to &nsim by 1; 
		U1=rand('Uniform'); 
		U2=rand('Uniform');
		
		start Func(x) global(U1, U2, alpha); 
		return(Exp(-((-Log(x))**alpha + (-Log(U1))**alpha)**(1/alpha)) *((-Log(x))**alpha + (-Log(U1))**alpha)**(-1 + 1/alpha)* ((-Log(U1))**(alpha -1))/U1-U2);
		finish;
		
		intervals = {0.00001 1};
		U2C = froot("Func", intervals);
		X=X//U1; 
		Y=Y//U2C; 
		YI=YI//U2; 
		end;
		
		Total=X||Y||YI;
		create GumC from Total [colname={'X','Y','YI'}]; 
		append from Total;
		close GumC;
	QUIT;
%mend SIM_Gum;

%SIM_Gum(nsim=1000, alpha=5, seed=12345);

/* a) */
/* From the plot we can already see there is some */
/* correlation. */
PROC CORR data=GUMC pearson kendall spearman 
		plots=scatter(ellipse=none); 
	var X Y;
RUN;
/* Kendall's Tau = 0.81311 and p-value < 0.001 */
/* Spearman's Rho = 0.95131 and p-value < 0.001 */

/* b) */
/* 𝜏(𝑥𝑖,𝑦𝑖) = 1−1/𝛼, We know 𝛼 = 5 for the simulation */
/* so, 𝜏 = 4/5 */

/* c) */
%MACRO SpearmanRho(rho=);
	PROC IML;
		pi = constant("pi");
		tau=&rho;
		
		start innerGum(y) global(alpha, x); 
		   return(Exp(-((-Log(x))**alpha + (-Log(y))**alpha)**(1/alpha)));
		finish; 
		
		start outerGum(par) global(x,alpha); 
			x=par;
		   yinterval = 0 || 1;
		   /** evaluate inner integral for the parameter value, a=x **/ 
		   call quad(w, "innerGum", yinterval);
		   return (w);
		finish; 
		
		start finalGum(param) global(alpha, tau);
		alpha=param;
		xinterval= {0 1};
		call quad(v, "outerGum", xinterval); /** outer integral **/ 
		return(12*v-(3+tau));
		finish;
		
		intervalsGum = {1 100};        
		SGum = froot("finalGum", intervalsGum);
		print(SGum);
		
		start innerClay(y) global(alpha, x); 
			return((x**(-alpha)+y**(-alpha)-1)**(-1/alpha));
		finish; 
		
		start outerClay(par) global(x, alpha); 
			x=par;
		
		if(alpha>0) then yinterval= 0||1;
		else yinterval= (1-x**(-alpha))**(-1/alpha)||1;
		   /** evaluate inner integral for the parameter value, a=x **/ 
		   call quad(w, "innerClay", yinterval);
		   return (w);
		finish; 
		
		start finalClay(param) global(alpha, tau);
		alpha=param;
		xinterval= {0 1};
		call quad(v, "outerClay", xinterval); /** outer integral **/ 
		return(12*v-(3+tau));
		finish;
		                 
		intervalsClay = {-1 10};        
		SClay = froot("finalClay", intervalsClay);
		print(SClay);
		
		SGau=2*sin(pi*tau/6);
		print(SGau);
		
		SFGM=3*tau;
		print(SFGM);
		
		start innerFrk(y) global(alpha, x); 
		return(-(1/alpha)*Log(1+(Exp(-alpha*x)-1)*(Exp(-alpha*y)-1)/(Exp(-alpha)-1)));
		finish; 
		
		start outerFrk(par) global(x, alpha); 
			x=par;
		   yinterval = 0 || 1;
		   /** evaluate inner integral for the parameter value, a=x **/ 
		   call quad(w, "innerFrk", yinterval);
		   return (w);
		finish; 
		
		start finalFrk(param) global(alpha, tau);
		alpha=param;
		xinterval= {0 1};
		call quad(v, "outerFrk", xinterval); /** outer integral **/ 
		return(12*v-(3+tau));
		finish;
         
		intervalsFrk = {-30 30};        
		SFrk = froot("finalFrk", intervalsFrk);
		print(SFrk);
		
		CPAR=SGum||SClay||SFrk||SGau||SFGM;
		
		create EstSpearman from CPAR [colname={'Gumbel alpha','Clayton alpha','Frank alpha','Gaussian rho','FGM alpha'}]; 
		append from CPAR;       
		close EstSpearman;
	QUIT;
%MEND SpearmanRho;

%MACRO KendallTau(tau=);
	PROC IML;
		pi = constant("pi");
		tau=&tau;
		
		SGum=1/(1-tau);
		print(SGum);
		SClay=2*tau/(1-tau);
		print(SClay);
		SGau=sin(pi*tau/2);
		print(SGau);
		SFGM=9*(tau/2);
		print(SFGM);
		
		start D(y);
		return(y/(Exp(y)-1));
		finish;
		
		*IF alpha>0 / tau>0;
		start FC(x) global(tau);
		dinterval=0||x;
		call quad(w, "D", dinterval);
		return(1-(4/x)*(1-(1/x)*w)-tau);
		finish;
		
		intervals = {0.00001 20};        
		SFrk = froot("FC", intervals);
		print(SFrk);
		
		
		CPAR=SGum||SClay||SFrk||SGau||SFGM;
		
		create EstKendall from CPAR [colname={'Gumbel alpha','Clayton alpha','Frank alpha','Gaussian rho','FGM alpha'}]; 
		append from CPAR;       
		close EstKendall;
	QUIT;
%MEND KendallTau;

%SpearmanRho(rho=0.95131);
/* Estimates an alpha of 5.411538 */

%KendallTau(tau=0.81311);
/* Estimates an alpha of 5.350741 */

/* d) */
%SpearmanRho(rho=0.95131);
/* Estimates an alpha of 18.473158 */

%KendallTau(tau=0.81311);
/* Estimates an alpha of 19.607393 */
		
/* e) */ 
%MACRO SIM_Frk(alpha=, nsim=, seed=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		do i=1 to &nsim by 1;
		U1=rand('Uniform'); 
		U2=rand('Uniform');
		
		start Func(x) global(U1,U2,alpha);
		return((Exp(alpha)*(-1 + Exp(alpha*x)))/(-Exp(alpha) + Exp(alpha*(1+x)) - Exp(alpha*(U1+x)) + Exp(alpha*(1 + U1)))-U2);
		finish;
		
		intervals = {0.00001 1};        
		U2C = froot("Func", intervals);
		
		X=X//U1;
		Y=Y//U2C;
		end;
		
		Total=X||Y;
		
		create FrkC from Total [colname={'X','Y'}]; 
		append from Total;       
		close FrkC;
	QUIT;
%MEND SIM_Frk;

/* Plot the results */
%SIM_Gum(nsim=1000, alpha=5, seed=12345);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbell's Copula, alpha = 5";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=1000, alpha=5.38, seed=12345);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbell's Copula, alpha = 5.38";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Frk(nsim=1000, alpha=19, seed=6789);
PROC SGPLOT data=FrkC aspect=1;
	title "Frank's Copula, alpha = 19";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

/* Question 3.2 */
DATA WEEK3_Q2;
	set SASDATA.IVF;
	IMP = IMP + (ranuni(1)-0.5); 
RUN;

/* a) */
/* Creates the Wide dataset */
PROC TRANSPOSE out=WEEK3_Q2_W(drop = _NAME_ _LABEL_) 
		data=WEEK3_Q2 prefix=IMP;
	by ID;
	id PER;
	var IMP; 
RUN;

/* Removes the missing values */
DATA WEEK3_Q2_W;
	set WEEK3_Q2_W;
	if cmiss(of _all_) then delete; 
RUN;

/* Creates uniform marginals, percentiles */
PROC RANK data=WEEK3_Q2_W out=WEEK3_Q2_R;
      var IMP4 IMP18;
      ranks rank_IMP4 rank_IMP18;
RUN;

PROC MEANS data=WEEK3_Q2_R N; 
	var rank_IMP4 rank_IMP18; 
RUN;

/* Marginals, use these to look at correlation */
DATA WEEK3_Q2_M;
	set WEEK3_Q2_R; 
	U_IMP4=rank_IMP4/237; 
	U_IMP18=rank_IMP18/237;
RUN;

/* Plot of the marginals and wide datasets */
ods graphics on;
proc corr data=WEEK3_Q2_W plots=scatter(ellipse=none); 
	var IMP4 IMP18;
run;

proc corr data=WEEK3_Q2_M plots=scatter(ellipse=none); 
	var U_IMP4 U_IMP18;
run;
ods graphics off;

/* b) */
PROC CORR data=WEEK3_Q2_W kendall spearman; 
	var IMP4 IMP18;
RUN;
/* Kendall's Tau = 0.03411 and p-value = 0.4352 */
/* Spearman's Rho = 0.05106 and p-value = 0.4349 */
/* We can't reject H0: p = 0, as p-value > 0.05 */

/* c) */
/* Estimate alphas parameters for the copulas */
%SpearmanRho(rho=0.05106);
%KendallTau(tau=0.03411);

/* As the estimates of don't differ too much we take */
/* the spearman's Rho as input */

/* Simulations of Copula's */
/* Addition with marginal dataset and uniform variable */
/* To be able to make a more accurate scatter plot */

/* Simulation of Gumbel's copula */
%MACRO SIM_Gum(alpha=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		use &dataset;
		read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
		U1=U[i];
		U2=rand('Uniform');
		
		start Func(x) global(U1,U2,alpha);
		return(Exp(-((-Log(x))**alpha + (-Log(U1))**alpha)**(1/alpha))*((-Log(x))**alpha + (-Log(U1))**alpha)**(-1 + 1/alpha)*((-Log(U1))**(alpha-1))/U1-U2);
		finish;
		
		intervals = {0.00001 1};        
		U2C = froot("Func", intervals);
		
		X=X//U1;
		Y=Y//U2C;
		YI=YI//U2;
		end;
		
		Total=X||Y||YI;
		
		create GumC from Total [colname={'X','Y','YI'}]; 
		append from Total;       
		close GumC;
	QUIT;
%MEND SIM_Gum;

/* Simulation of the Gaussian copula */
%MACRO SIM_GC(rho=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		rho=&rho;
		
		use &dataset;
		read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
		U1=U[i];
		U2=rand('Uniform');
		
		start Func(x) global(U1,U2,rho);
		return(CDF('Normal',quantile('NORMAL', x),rho*quantile('NORMAL',U1),(1-rho**2))-U2);
		finish;
		
		intervals = {0.00001 0.99999};        
		U2C = froot("Func", intervals);
		
		X=X//U1;
		Y=Y//U2C;
		YI=YI//U2;
		end;
		
		Total=X||Y||YI;
		
		create GC from Total [colname={'X','Y','YI'}]; 
		append from Total;       
		close GC;
	QUIT;
%MEND SIM_GC;

/* Simulation of Clayton's copula */
%MACRO SIM_Clay(alpha=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		use &dataset;
		read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
		U1=U[i];
		U2=rand('Uniform');
		
		start Func(x) global(U1,U2,alpha);
		return(U1**(-1 -alpha)*(x**(-alpha) + U1**(-alpha)-1)**(-1 - 1/alpha)-U2);
		finish;
		
		intervals = {0.001 1};        
		U2C = froot("Func", intervals);
		
		X=X//U1;
		Y=Y//U2C;
		end;
		
		Total=X||Y;
		
		create CC from Total [colname={'X','Y'}]; 
		append from Total;       
		close CC;
	QUIT;
%MEND SIM_Clay;

/* Simulation of Frank's copula */
%MACRO SIM_Frk(alpha=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		use &dataset;
		read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
		U1=U[i];
		U2=rand('Uniform');
		
		start Func(x) global(U1,U2,alpha);
		return((Exp(alpha)*(-1 + Exp(alpha*x)))/(-Exp(alpha) + Exp(alpha*(1+x)) - Exp(alpha*(U1+x)) + Exp(alpha*(1 + U1)))-U2);
		finish;
		
		intervals = {0.00001 1};        
		U2C = froot("Func", intervals);
		
		X=X//U1;
		Y=Y//U2C;
		end;
		
		Total=X||Y;
		
		create FrkC from Total [colname={'X','Y'}]; 
		append from Total;       
		close FrkC;
	QUIT;
%MEND SIM_Frk;

/* Simulation of FGM's copula */
%MACRO SIM_FGM(alpha=, nsim=, seed=, dataset=, uvar=);
	PROC IML;
		call streaminit(&seed);
		alpha=&alpha;
		
		use &dataset;
		read all var{&uvar};
		close &dataset;
		U=&uvar;
		
		do i=1 to &nsim by 1;
		U1=U[i];
		U2=rand('Uniform');
		
		start Func(x) global(U1,U2,alpha);
		return(x*(1 + alpha*(1 - x)*(1 - U1)) - alpha*(1 - x)*x*U1-U2);
		finish;
		
		intervals = {0.00001 1};        
		U2C = froot("Func", intervals);
		
		X=X//U1;
		Y=Y//U2C;
		YI=YI//U2;
		end;
		
		Total=X||Y||YI;
		
		create FGMC from Total [colname={'X','Y','YI'}]; 
		append from Total;       
		close FGMC;
	QUIT;
%MEND SIM_FGM;

%SIM_GC(nsim=236, rho=0.0534635, seed=6789, dataset=WEEK3_Q2_M, uvar=U_IMP4);
PROC SGPLOT data=GC aspect=1;
	title "Gaussian Copula";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=236, alpha=1.0353378, seed=6789, dataset=WEEK3_Q2_M, uvar=U_IMP4);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Clay(nsim=236, alpha=0.0705201, seed=6789, dataset=WEEK3_Q2_M, uvar=U_IMP4);
PROC SGPLOT data=CC aspect=1;
	title "Clayton's Copula";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Frk(nsim=236, alpha=0.3067442, seed=6789, dataset=WEEK3_Q2_M, uvar=U_IMP4);
PROC SGPLOT data=FrkC aspect=1;
	title "Frank's Copula";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_FGM(nsim=236, alpha=0.15318, seed=6789, dataset=WEEK3_Q2_M, uvar=U_IMP4);
PROC SGPLOT data=FGMC aspect=1;
	title "FGM's Copula";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

/* Original scatter plot */
ods graphics on;
proc corr data=WEEK3_Q2_M plots=scatter(ellipse=none); 
	title "Original margin's";
	var U_IMP4 U_IMP18;
run;
ods graphics off;

/* Question 3.3 */
DATA WEEK3_Q3;
	do i = 1 to 1000;
	   X = rand('Uniform');
	   Y = rand('Uniform');
	   output;
	end;
RUN;

PROC PRINT data=WEEK3_Q3;  
	var X Y;  
RUN;

/* a) */
/* if X and Y are independent, */
/* P(X =< 0.7, Y =< 0.7) = F_1(0.7)*F_2(0.7) */
/* 1.2 * 1.2 = 1.44 */

PROC CORR data=WEEK3_Q3_a kendall spearman; 
	var X Y;
RUN;

/* b) */
PROC SGPLOT data=WEEK3_Q3_a aspect=1;
	title "Gumbel's Copula, real data";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=1000, alpha=1, seed=6789);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula, alpha=1";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=1000, alpha=2, seed=6789);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula, alpha=2";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=1000, alpha=5, seed=6789);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula, alpha=5";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=1000, alpha=10, seed=6789);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula, alpha=10";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;
/* Most similair to the Gumbell copula with alpha = 1 */

/* c) */
%SIM_Clay(nsim=1000, alpha=2, seed=6789);
PROC SGPLOT data=CC aspect=1;
	title "Clayton's Copula, alpha=2";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Clay(nsim=1000, alpha=5, seed=6789);
PROC SGPLOT data=CC aspect=1;
	title "Clayton's Copula, alpha=5";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Clay(nsim=1000, alpha=10, seed=6789);
PROC SGPLOT data=CC aspect=1;
	title "Clayton's Copula, alpha=10";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

/* d) */

/* Question 3.4 */
DATA WEEK3_Q4;
	set SASDATA.RCT;
RUN;

/* a) */















