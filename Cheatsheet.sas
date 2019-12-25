/* Week 1: Intervals */

/* Week 2: Independent Samples */

/* Mann-Whitney U Test */
%MACRO mann_whitney_u(dataset, class, var);
	ods select none;
	PROC NPAR1WAY data=&dataset wilcoxon;
		var &var;
		class &class;
		ods output WilcoxonScores=OUT_SCORES(
			rename=(SumOfScores=S));
		ods output WilcoxonTest=OUT_TEST
		ods output KruskalWallisTest=OUT_KRUS;
	RUN;
	ods select all;
	
	/* Chi square p-value */
	PROC SQL;
		CREATE TABLE P_KRUS AS
			SELECT Prob FROM OUT_KRUS;
	RUN;
	
	/* Wilcoxon p-value */
	PROC SQL;
		CREATE TABLE P_TABLE AS
			SELECT tProb2 FROM OUT_TEST;
	RUN;
	
	DATA OUT_SCORES;
		set OUT_SCORES; 
		CLASS_ID = _N_ - 1;
	RUN;
	
	PROC TRANSPOSE data=OUT_SCORES
			out=OUT_N(drop=_NAME_) prefix=N;
		id CLASS_ID;
		var N;
	RUN;
	
	PROC TRANSPOSE data=OUT_SCORES
			out=OUT_S(drop=_NAME_ _LABEL_) PREFIX=S;
		id CLASS_ID;
		var S;
	RUN;
	
	DATA RESULT;
		merge OUT_N OUT_S P_TABLE P_KRUS;
		P_VALUE = tProb2;
		P_KRUS = Prob;
		U0 = S0 - N0 * (N0+1)/2;
		U1 = S1 - N1 * (N1+1)/2;
		P0 = U0 / (N0*N1);
		P1 = U1 / (N0*N1);
	RUN;
	
	title "Mann Whitney U test";
	PROC PRINT data=OUT_SCORES label noobs;
		var CLASS_ID CLASS; 
		label CLASS_ID="class"
		CLASS="group identifier";
	RUN;

	PROC PRINT data=RESULT label;
		var P_VALUE P_KRUS U0 U1 P0 P1;
		label P_VALUE = "p-value Wilcoxon Test"
		P_KRUS= "p-value Kruskal-Wallis Test"
		U0="statistic (U0)" 
		U1="statistic (U1)" 
		P0="P(class0 > class1)" 
		P1="P(class0 <= class1)";
	RUN;
	title;
%MEND;

%MACRO binary_hypothesis(dataset, var, class); 
	PROC MEANS data=&dataset n sum noprint;
		var &var;
		class &class;
		output out=OUT n=N sum=COUNT; 
	RUN;

	DATA OUT0;
		set OUT;
		COUNT0 = COUNT; 
		N0 = N;
		P0 = COUNT0 / N0; 
		where &class = 0; 
		keep COUNT0 N0 P0; 
	RUN;
	
	DATA OUT1;
		set OUT;
		COUNT1 = COUNT; 
		N1 = N;
		P1 = COUNT1 / N1; 
		where &class = 1; 
		keep COUNT1 N1 P1; 
	RUN;
    
    DATA OUT;
		merge OUT0 OUT1;
		P = (COUNT0 + COUNT1) / (N0 + N1);
		STAT = (P0 - P1) / sqrt(P * (1-P) * (1/N0 + 1/N1));
		CHISQ = STAT **2;
		P_VALUE = 2*min(cdf("normal", STAT, 0, 1), 1-cdf("normal", STAT, 0, 1)); 
	RUN;
	
	PROC PRINT data=OUT; 
		var STAT CHISQ P_VALUE; 
	RUN;
%MEND;

/* Week 3: Dependence */

/* Copula simulation macro's */
/* Addition with marginal dataset and uniform variable */
/* To be able to make a more accurate scatter plot */

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

/* Estimate parameters copula's */
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
