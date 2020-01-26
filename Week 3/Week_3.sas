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
PROC CORR data=WEEK3_Q2_W plots=scatter(ellipse=none); 
	var IMP4 IMP18;
RUN;

PROC CORR data=WEEK3_Q2_M plots=scatter(ellipse=none); 
	var U_IMP4 U_IMP18;
RUN;
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

/* a) */
/* if X and Y are independent, */
/* P(X =< 0.7, Y =< 0.7) = F_1(0.7)*F_2(0.7) */
/* 0.7 * 0.7 = 0.49 */

PROC CORR data=WEEK3_Q3 kendall spearman 
	plots=scatter(ellipse=none); 
	var X Y;
RUN;

/* Creates uniform marginals, percentiles */
PROC RANK data=WEEK3_Q3 out=WEEK3_Q3_R;
      var X Y;
      ranks X Y;
RUN;

PROC MEANS data=WEEK3_Q3_R N; 
	var X Y; 
RUN;

/* Marginals, use these to look at correlation */
DATA WEEK3_Q3_M;
	set WEEK3_Q3_R; 
	U_X=X/1000; 
	U_Y=Y/1000;
RUN;

/* b) */
PROC SGPLOT data=WEEK3_Q3_M aspect=1;
	title "real data";
	scatter x=U_X y=U_Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=1000, alpha=1, seed=6789, dataset=WEEK3_Q3_M, uvar=U_X);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula, alpha=1";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=1000, alpha=2, seed=6789, dataset=WEEK3_Q3_M, uvar=U_X);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula, alpha=2";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=1000, alpha=5, seed=6789, dataset=WEEK3_Q3_M, uvar=U_X);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula, alpha=5";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Gum(nsim=1000, alpha=10, seed=6789, dataset=WEEK3_Q3_M, uvar=U_X);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula, alpha=10";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;
/* Most similair to the Gumbel's copula with alpha = 1 */

/* c) */
PROC SGPLOT data=WEEK3_Q3_M aspect=1;
	title "real data";
	scatter x=U_X y=U_Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Clay(nsim=1000, alpha=2, seed=6789, dataset=WEEK3_Q3_M, uvar=U_X);
PROC SGPLOT data=CC aspect=1;
	title "Clayton's Copula, alpha=2";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Clay(nsim=1000, alpha=5, seed=6789, dataset=WEEK3_Q3_M, uvar=U_X);
PROC SGPLOT data=CC aspect=1;
	title "Clayton's Copula, alpha=5";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

%SIM_Clay(nsim=1000, alpha=10, seed=6789, dataset=WEEK3_Q3_M, uvar=U_X);
PROC SGPLOT data=CC aspect=1;
	title "Clayton's Copula, alpha=10";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

/* d) - e) */
/* TODO */

/* Question 3.4 */
DATA WEEK3_Q4;
	set SASDATA.RCT;
RUN;

/* a) */
PROC TRANSPOSE out=WEEK3_Q4_W(drop = _NAME_ _LABEL_) 
		data=WEEK3_Q4 prefix=RESP;
	by ID;
	id TIME;
	var RESP; 
RUN;

/* If we remove missing values the results don't match */
/* the answers. */
/* DATA WEEK3_Q4_W; */
/* 	set WEEK3_Q4_W; */
/* 	if cmiss(of _all_) then delete;  */
/* RUN; */

/* b) */
PROC CORR data=WEEK3_Q4_W plots=scatter(ellipse=none) 
		spearman pearson kendall;
	var RESP1 RESP2;
RUN;
/* Pearson rho = 0.47800 and p-value < 0.001 */
/* Kendall tau = 0.35941 and p-value < 0.001 */
/* Spearman rho = 0.50477 and p-value < 0.001 */
/* Thus, we reject the H0 */

/* Creating margins to see the correlation plot */
PROC RANK data=WEEK3_Q4_W out=WEEK3_Q4_R;
	var RESP1 RESP2;
    ranks rank_RESP1 rank_RESP2;
RUN;

PROC MEANS data=WEEK3_Q4_R N; 
	var rank_RESP1 rank_RESP2;
RUN;

DATA WEEK3_Q4_M;
	set WEEK3_Q4_R; 
	U_RESP1=rank_RESP1/716; /*694; After cleaning */
	U_RESP2=rank_RESP2/716; /*694; */
RUN;

PROC CORR data=WEEK3_Q4_M plots=scatter(ellipse=none);
	var U_RESP1 U_RESP2;
RUN;

/* c) */
/* These are already given. */

/* d) */
/* Bias adjustment = no is important */
PROC CORR data=WEEK3_Q4_W pearson fisher(biasadj=no);
	var RESP1 RESP2;
RUN;
/* The CI is given by (0.419252, 0.532766) */

/* e) */
%SpearmanRho(rho=0.50477);
%KendallTau(tau=0.35941);
/* Both alpha = 1.51431 and alpha = 1.617345 are outside */
/* of the supported range for FGM's Copula */

/* f) */
/* The alpha = 1.1221218 for Clayton's Copula using tau */

/* g) */
/* The alpha = 1.0935697 for Clayton's Copula using rho */
/* This is within reasonable range of f. */

/* h) */
/* Taking the averages for Clayton's alpha = 1.10784575 */
/* and frank's alpha = 3.5580449 */
%SIM_Clay(nsim=700, alpha=1.10784575, seed=6789, dataset=WEEK3_Q4_W, uvar=RESP1);
%SIM_Frk(nsim=700, alpha=3.5580449, seed=6789, dataset=WEEK3_Q4_W, uvar=RESP1);
/* The simulation fails after 324 data points */

PROC SGPLOT data=WEEK3_Q4_M aspect=1;
	title "Original data";
	scatter x=U_RESP1 y=U_RESP2 / markerattrs=(color='blue' size=6);
RUN;

PROC SGPLOT data=FrkC aspect=1;
	title "Frank's Copula";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;

PROC SGPLOT data=CC aspect=1;
	title "Clayton's Copula";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;
/* We can see it resembles Frank's Copula better, as Clayton's Copula */
/* has a heavier one sided correlation tail. */

/* Add noise to the discrete variables, to have a better overview */
DATA WEEK3_Q4_M_h;
	SET WEEK3_Q4_M;
	U_RESP1=U_RESP1 + 0.1*(ranuni(1)-0.5); 
	U_RESP2=U_RESP2 + 0.1*(ranuni(1)-0.5); 
RUN;

PROC SGPLOT data=WEEK3_Q4_M_h aspect=1;
	title "Original data";
	scatter x=U_RESP1 y=U_RESP2 / markerattrs=(color='blue' size=6);
RUN;

/* Question 3.6 */
DATA WEEK3_Q6;
	set SASDATA.RCT;
	WHERE center = 1;
RUN;

/* a) */
PROC TRANSPOSE out=WEEK3_Q6_W(drop = _NAME_ _LABEL_) 
		data=WEEK3_Q6 prefix=RESP;
	by ID;
	id TIME;
	var RESP; 
RUN;

/* Remove observations with missing values */
/* DATA WEEK3_Q6_W; */
/* 	set WEEK3_Q6_W; */
/* 	if cmiss(of _all_) then delete;  */
/* RUN; */

DATA WEEK3_Q6_W;
	set WEEK3_Q6_W;
	drop RESP1 RESP2 RESP3 RESP4;
	Z_DIFF = RESP6 - RESP5;
	Z_LDIFF = log(RESP6) - log(RESP5);
	Z_RATIO = RESP6/RESP5;
RUN;

ods select histogram;
PROC UNIVARIATE data=WEEK3_Q6_W;
	histogram Z_DIFF/normal;
    histogram Z_LDIFF/normal;
    histogram Z_RATIO/normal;
RUN;

/* b) */
/* Testing for H0: μ = 0 vs μ != 0 */
ods select TestsForLocation;
PROC UNIVARIATE data=WEEK3_Q6_W normal;
	var Z_DIFF;
RUN;

ods select TestsForLocation;
PROC UNIVARIATE data=WEEK3_Q6_W normal;
	var Z_LDIFF;
RUN;

/* Testing for H0: μ = 1 vs μ != 1 */
ods select TestsForLocation;
PROC UNIVARIATE data=WEEK3_Q6_W normal MU0=1;
	var Z_RATIO;
RUN;
/* So we can't reject H0 for all 3 */

/* c) */
/* For the paired t-test, we assume normality.  */
/* Normality does not seem apparent for the ratio. */

/* For the sign test, no assumptions are made. */

/* For the Wilcoxon signed rank test, we need symmetric  */
/* distributions. The his- tograms suggest slightly  */
/* left-skewed distributions, but nothing extreme. */

/* d) */
/* Based on power, */
/* t-test, Wilcoxon-signed rank test, Sign test. */

/* e) */
/* The t-test is not directly appropriate in this case.  */
/* However, if the sample size is large one could rely on  */
/* the CLT.  */

/* The assumption of symmetric differences is  */
/* violated so the Wilcoxon-signed rank test cannot be used  */
/* to test equality in median of both groups.  */

/* The sign-test does not rely on any distributional  */
/* assumptions, thus can be used here to test for equality  */
/* in medians of the groups. */

/* Question 3.7 */
/* a) */
DATA WEEK3_Q7;
	set SASDATA.IVF;
	where PER = 4 or PER = 18;
RUN;

PROC TRANSPOSE out=WEEK3_Q7_wide(drop = _NAME_ _LABEL_) 
		data=WEEK3_Q7 prefix=IMP;
	by ID;
	id PER;
	var IMP; 
RUN;

DATA WEEK3_Q7_wide;
	set WEEK3_Q7_wide;
	DIFF = IMP4 - IMP18;
	RATIO = IMP4/IMP18;
	RATIO_0 = RATIO-1;
	LDIFF = log(RATIO);
RUN;

PROC UNIVARIATE data=WEEK3_Q7_wide normaltest;
	var DIFF LDIFF RATIO;
	histogram /normal;
    ods select histogram;
RUN;
/* LDIFF and RATIO seem to approximate a normal  */
/* distribution. DIFF does not approximate a normal */
/* distribution and is not symmetric. */

/* b) */
PROC UNIVARIATE data=WEEK3_Q7_wide mu0=0 0 1;
	var DIFF LDIFF RATIO;
RUN;
/* c) */
/* For DIFF only the Sign-test is reliable which  */
/* rejects that the medians are equal */
/* For LDIFF and RATIO all 3 tests are appropriate */
/* and indicate that the medians and means are not equal */

/* Question 3.8 */
DATA COAG; 
	input Patient C K@@;
	datalines;
	1 120 132 8 145 133 15 117 123
	2 114 116 9 120 123 16 125 108
	3 129 135 10 129 116 17 136 131 
	4 128 115 11 126 127 18 151 119
	5 155 134 12 136 140 19 130 129 
	6 105 56 13 135 140 20 136 124
	7 114 114 14 125 114 21 113 112
	;
RUN;

/* a) */
PROC CORR data=COAG plots=scatter(ellipse=none) 
		spearman pearson kendall;
	var C K;
RUN;
/* Under H0 Rho = 0, */
/* Pearson's rho = 0.59084 and p-value = 0.0048 */
/* Spearman's rho = 0.66134 and p-value = 0.0011 */
/* Kendall's tau = 0.51346 and p-value = 0.0014 */
/* So we reject the H0 */

/* b) */
/* Creates uniform marginals, percentiles */
PROC RANK data=COAG out=COAG_R;
      var C K;
      ranks rank_C rank_K;
RUN;

PROC MEANS data=COAG_R N; 
	var rank_C rank_K; 
RUN;

/* Marginals, use these to look at correlation */
DATA COAG_M;
	set COAG_R; 
	U_C=rank_C/21; 
	U_K=rank_K/21;
RUN;

PROC SGPLOT data=COAG_M aspect=1;
	title "Actual data";
	scatter x=U_C y=U_K / markerattrs=(COLOR='blue' size=6);
RUN;

%SIM_Gum(nsim=21, alpha=1.0353378, seed=6789, dataset=COAG_M, uvar=U_C);
PROC SGPLOT data=GumC aspect=1;
	title "Gumbel's Copula";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;
%SIM_FGM(nsim=21, alpha=0.15318, seed=6789, dataset=COAG_M, uvar=U_C);
PROC SGPLOT data=FGMC aspect=1;
	title "FGM's Copula";
	scatter x=X y=Y / markerattrs=(color='blue' size=6);
RUN;
/* Looks more alike to the FGM copula */

/* c) */
DATA COAG;
	set COAG;
	Z_DIFF = C - K;
	Z_LDIFF = log(C) - log(K);
	Z_RATIO = C/K;
	Z0_RATIO = C/K - 1;
RUN;

ods select histogram;
PROC UNIVARIATE data=COAG normal;
	var Z_DIFF Z_LDIFF Z_RATIO;
	histogram Z_DIFF/normal;
    histogram Z_LDIFF/normal;
    histogram Z_RATIO/normal;
RUN;
/* Z_DIFF looks most like a normal distribution */
/* But the evidence is not that convincing. */

/* Testing for H0: μ = 0 vs μ != 0 */
ods select TestsForLocation;
PROC UNIVARIATE data=COAG normal;
	var Z_DIFF;
RUN;

ods select TestsForLocation;
PROC UNIVARIATE data=COAG normal;
	var Z_LDIFF;
RUN;

/* Testing for H0: μ = 1 vs μ != 1 */
ods select TestsForLocation;
PROC UNIVARIATE data=COAG normal MU0=1;
	var Z_RATIO;
RUN;

/* d) */
DATA COAG;
	set COAG;
	TT_C = (C<120);
	TT_K = (K<120);
RUN;

/* a Pearson correlation coefficient estimated for two */
/* binary variables will return the phi coefficient. */
PROC CORR data=COAG pearson;
	var TT_C TT_K;
RUN;
/* We get Test statistic = 0.41957 and p-value = 0.0583 */
/* We can't reject H0 */

/* e) */
PROC FREQ data=COAG;
      tables TT_C*TT_K; 
      exact mcnem;
RUN;

/* Question 3.9 */
DATA WEEK3_Q9;
	set SASDATA.IVF;
RUN;

/* Creates the wide dataset */
PROC TRANSPOSE out=WEEK3_Q9_W(drop = _NAME_ _LABEL_) 
		data=WEEK3_Q9 prefix=IMP;
	by ID;
	id PER;
	var IMP; 
RUN;

DATA WEEK3_Q9_W;
	set WEEK3_Q9_W;
/* 	if IMP10<85 or IMP18<85 then delete; */
	if cmiss(of IMP10 IMP18) then delete;
RUN;

/* a) */
DATA WEEK3_Q9_a;
	set WEEK3_Q9_W;
	IMP_DIFF =  IMP18 - IMP10;
/* 	IMP_RATIO = IMP18 / IMP10; */
/* 	IMP_0RATIO = IMP_RATIO - 1; */
/* 	IMP_LDIFF = log(IMP_RATIO); */
/* 	IMP_SIGN = (IMP18 >= IMP10); */
RUN;

/* Sign Test */
/* 2 ⋅ min(𝑃(𝑆 ≤ 𝑠),𝑃(𝑆 ≥ 𝑠)), for UPL and LPL */
PROC IML;
	use WEEK3_Q9_a;
	read all var{IMP_DIFF};
	CLOSE WEEK3_Q9_a;
	
	plus = sum((IMP_DIFF>0));
	min = sum((IMP_DIFF<0));
	zero = sum((IMP_DIFF=0));
	total = plus+min;
	print(total||min||plus||zero);
	
	LPL = cdf("Binom", plus, 0.5, total);
	UPL = 1 - cdf("Binom", plus-1, 0.5, total);

	print(LPL||UPL);
	A=LPL||UPL;
	
	create SignT from A [colname={'LPL' 'UPL'}]; 
	append from A;       
	close SignT;
QUIT;

PROC FREQ data=WEEK3_Q9_a;
	tables SIGN/binomial(P=0.5);
	exact binomial;
RUN;

/* b) */
PROC FREQ data=WEEK3_Q9_W;
      tables NP10*NP18; 
      exact mcnem;
RUN;

/* c) */
DATA WEEK3_Q9_W; 
	set WEEK3_Q9_W;
	if cmiss (of IMP10 IMP18 ) then delete;
RUN;

PROC FREQ data=WEEK3_Q9_W;
      tables NP10*NP18; 
      exact mcnem;
RUN;

/* d) */
/* - */

/* Question 3.11 */
/* Drug A = 0, Drug B = 1 */
DATA WEEK3_Q11;
	input treatment high count drug;
	datalines;
	1 1 12 0
	0 1 11 0
	1 0 7 0
	0 0 70 0
	1 1 10 1
	0 1 13 1
	1 0 4 1
	0 0 73 1
	;
RUN;

/* a) */
/* 2 sided, 1 sided use tables A*B/AGREE; */
PROC FREQ data=WEEK3_Q11;
	where drug = 0;
	tables high*treatment;
	weight count;
	exact mcnem;
RUN;
/* Test statistic = 0.8889 and p-value = 0.4807 */

/* b) */
/* - */

/* c) */
PROC FREQ data=WEEK3_Q11;
	where drug = 0;
	tables high*treatment /chisq;
	weight count;
	exact chisq;
RUN;
/* The performed Chi-square test would test  */
/* if the probability to have a high blood  */
/* pressure after the treatment would be  */
/* independent of the blood pressure before  */
/* the treatment. */
