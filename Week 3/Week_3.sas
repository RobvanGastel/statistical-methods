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
PROC CORR data=GUMC kendall; 
	var X Y;
RUN;
/* Kendall's Tau = 0.81311 */

PROC CORR data=GUMC spearman; 
	var X Y;
RUN;
/* Spearman's Rho = 0.95131 */

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
PROC TRANSPOSE out=WEEK3_Q2_W(drop = _NAME_ _LABEL_) 
		data=WEEK3_Q2 prefix=IMP;
	by ID;
	id PER;
	var IMP; 
RUN;

DATA WEEK3_Q2_W;
	set WEEK3_Q2_W;
	if cmiss(of _all_) then delete; 
RUN;

/* Create marginals */
PROC RANK data=WEEK3_Q2_W out=ranked_a;
      var IMP4 IMP10;
      ranks rank_IMP4 rank_IMP10;
RUN;

PROC MEANS data=ranked_a N; 
	var rank_IMP4 rank_IMP10; 
RUN;

DATA marginals;
	set ranked_a; 
	U_IMP4=rank_IMP4/236; 
	U_IMP10=rank_IMP10/236;
RUN;

/* b) */
ods graphics on;
PROC CORR data=WEEK3_Q2_W plots=scatter(ellipse=none); 
	var IMP10 IMP18;
RUN;
ods graphics off;

PROC CORR data=marginals kendall; 
	var IMP10 IMP18;
RUN;
/* Kendall's Tau = 0.09282 */

PROC CORR data=marginals spearman; 
	var IMP10 IMP18;
RUN;
/* Spearman's Rho = 0.13739 */

/* c) */




