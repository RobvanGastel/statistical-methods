/* Week 1: Intervals */
/* DESC: Used to approximate the power of test with ns runs and n samples */
%MACRO samples(dataset=, ns=, n=);
	PROC SURVEYSELECT data=&dataset noprint method=urs n=&n out=FINAL;
	RUN;
	
	DATA FINAL;
		set FINAL; 
		sampleno = 1;
	RUN;
	
	%do sn = 2 %to &ns;
		PROC SURVEYSELECT data=&dataset noprint
		method=urs n=&n out=SAMPLEI;
	RUN;
	
	DATA SAMPLEI;
		set SAMPLEI;
		sampleno = &sn;
	RUN;
	
	DATA FINAL;
		set Final SAMPLEI;
	RUN;
	%end;
	
	PROC DATASETS library=work noprint;
		delete SAMPLEI;
	RUN;
%MEND;

/* Week 2: Independent Samples */
/* DESC: Two sample binary comparison using Chi-Square test */
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

/* Mann-Whitney U Test */
/* DESC: Test extension on WRS */
%MACRO Mann_Whitney_U(dataset, class, var);
	ods select none;
	PROC NPAR1WAY data=&dataset wilcoxon correct=NO;
		var &var;
		class &class;
		ods output WilcoxonScores=OUT_SCORES(
			rename=(SumOfScores=S));
		ods output WilcoxonTest=OUT_TEST;
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
			SELECT Prob2 FROM OUT_TEST;
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
		P_VALUE = Prob2;
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

/* Week 4: ANOVA and randomness */

/* Important! Before a Runs test sort the input */
/* Runs test conditional and unconditional */
%MACRO Runs_test(data=, var=, alpha=);
	PROC IML;
		use &data;
		read all var {&var};
		close &data;
		
		X=&var;
		n=nROW(X);
		MED=median(X);
		
		XC=X;
		do i=1 to n by 1;
			if (XC[i] >= MED) then XC[i]=1;
			else XC[i]=0;
		end;
		
		n1C=sum(XC);
		n2C=n-n1C;
		
		RC=1;
		do i=2 to n by 1;
			if(XC[i] ^= XC[i-1]) then RC=RC+1;
		end;
		
		MUC=1+(2*n1C*n2C)/(n1C+n2C);
		VARC=2*n1C*n2C*(2*n1C*n2C-n1C-n2C)/((n1C+n2C-1)*(n1C+n2C)**2);
		
		SC=(RC-MUC)/sqrt(VARC);
		TC=quantile('NORMAL',&alpha/2);
		TCU=quantile('NORMAL',1-&alpha/2);
		PC=(1-cdf('NORMAL',abs(SC)))*2;
		
		XUC=repeat(0,n-1,1);
		TIES=0;
		do i=1 to (n-1) by 1;
			if (X[i+1] > X[i]) then XUC[i]=1;
			if (X[i+1] < X[i]) then XUC[i]=0;
			if (X[i+1] = X[i]) then XUC[i]=XUC[i-1];
			if (X[i+1] = X[i]) then TIES=TIES+1;
		end;
		
		RUC=1;
		do i=2 to (n-1) by 1;
			if(XUC[i] ^= XUC[i-1]) then RUC=RUC+1;
		end;
		
		MUUC=(2*(n-TIES)-1)/3;
		VARUC=(16*(n-TIES)-29)/90;
		
		SUC=(RUC-MUUC)/sqrt(VARUC);
		TUC=quantile('NORMAL',&alpha/2);
		TUCU=quantile('NORMAL',1-&alpha/2);
		PUC=(1-cdf('NORMAL',abs(SUC)))*2;
		
		A = RC||MUC||sqrt(VARC)||PC||SC||TC||TCU||n;
		create CondRuns from A [colname={'runs','mean runs','std runs','p-value','Normalized statistic','Critical values L','Critical values U','n'}];
		append from A; 
		
		B = RUC||MUUC||sqrt(VARUC)||PUC||SUC||TUC||TUCU||n-TIES;
		create UncondRuns from B [colname={'runs','mean runs','std runs','p-value','Normalized statistic','Critical values L','Critical values U','n'}]; 
		append from B;
		
		C = ties;
		create Colties from C [colname={'Ties'}];
		append from C;
		close UncondRuns CondRuns Colties;
				
		PROC PRINT data=CondRuns noobs;
			title "Median based (conditional) runs test";
		RUN;
		
		PROC PRINT data=UncondRuns noobs;
			title "(unconditional) runst test for serial randomness";
		RUN;
		
		PROC PRINT data=Colties noobs;
			title "Number of ties for both";
		RUN;
	QUIT;
%MEND;

/* Week 5: Normality and Outliers */

/* DESC: Skewness and Kurtosis tests for normality */
%MACRO Skewness_Kurtosis_test(skewness=, kurtosis=, n=);
	DATA approx;
		N=&n;
		G1=&skewness;
		G2=&kurtosis;
		b1=(N-2)*G1/(sqrt(N*(N-1)));
		b2=G2*((N-2)*(N-3))/((N+1)*(N-1))+3*(N-1)/(N+1);

		Cn=(3*(N**2+27*N-70)*(N+1)*(N+3))/((N-2)*(N+5)*(N+7)*(N+9));
		Wn2=-1+SQRT(2*(Cn-1));
		
		Alphan=SQRT(2/(Wn2-1));
		Dn=1/sqrt(log(sqrt(Wn2)));
		Bn=sqrt((N+1)*(N+3)/(6*(N-2)))*b1;
		Ts=Dn*log(Bn/Alphan+sqrt(1+(Bn/Alphan)**2));
		
		Mun=3*(N-1)/(N+1);
		Sigman=sqrt((24*N*(N-2)*(N-3))/((N+3)*(N+5)*(N+1)**2));
		Gamma1n=((6*(N**2-5*N+2))/((N+7)*(N+9)))*sqrt(6*(N+3)*(N+5)/(N*(N-2)*(N-3)));
		An=6+(8/(Gamma1n))*(2/Gamma1n+sqrt(1+4/(Gamma1n**2)));
		Un=(b2-Mun)/Sigman;
		Tk=sqrt(9*An/2)*((9*An-2)/(9*An)-((1-2/An)/(1+Un*sqrt(2/(An-4))))**(1/3));
		K2=Tk**2+Ts**2;
		
		Ps=2*min(cdf('Normal',Ts,0,1),1-cdf('Normal',Ts,0,1));
		Pk=2*min(cdf('Normal',Tk,0,1),1-cdf('Normal',Tk,0,1));
		PK2=1-cdf('chisq',K2,2);
	RUN;
	
	PROC PRINT data=approx noobs;
	RUN;
%MEND;

/* Grubbs test for outliers */
/* For n =< 100 look up the critical value in the table */
/* After sorting if the first value > approx or absolute */
/* value, we reject H0 and determine there is a outlier. */
%MACRO Grubbs_test(dataset, var, id, alpha=0.05);
    ods select none;
	PROC MEANS data=&dataset mean var n;
		var &var;
	    output out=out mean=mean var=var n=n; 
	RUN;
	
	DATA outliers;
		set &dataset(keep=&id &var);
		if _n_=1 then set out;
		/* statistic */
		u = abs((&var - mean) / sqrt(var));
		/* critical value */
		t = quantile("t", &alpha / (2*n), n-2);
		c = (n-1) * sqrt(t**2 / (n * (t**2 + n - 2))); 
		/* check if this is an outlier */
		if(u > c) then outlier = "yes"; 
		else outlier = "no"; /* p-value */
		u_inv = u*sqrt((n-2)*n) / sqrt(1 - (u**2-(n-2))*n); 
		p_value = min(2*n*(1-cdf("t", u_inv, n-2)), 1); 
		keep &id p_value u c outlier;
	RUN;
	
	PROC SORT data=outliers; 
		by descending u;
	RUN;
	ods select all;
	              
	title1 "Grubbs test";
	title2 "Variable: &var";
	PROC PRINT data=outliers(obs=1) label noobs;
		var &id p_value u c outlier; 
		label &id="id" p_value="p-value"
		u="statistic (|u|)" c="critical value" 
		outlier="is outlier?";
    RUN;
    title;
%MEND;

/* Tukey's method */
%MACRO Tukey_method(dataset, var, id);
	ods select none;
	PROC MEANS data=&dataset median p25 p75;
		var &var;
		output out=quartiles p25=p25 p75=p75;
	RUN;
	
	DATA tukey;
		set &dataset(keep=&id &var);
		if _n_=1 then set quartiles;
		iqr = p75-p25;
		lower = p25 - 1.5*IQR;
		upper = p75 + 1.5*IQR;
		if &var >= lower and &var <= upper then delete;
	RUN;
	ods select all;
	
	title1 "Outliers found using tukey’s method."; 
	title2 "Variable: &var";
	PROC PRINT data=tukey noobs label;
		var &id &var lower upper iqr; 
		label lower="Lower limit" upper="Upper limit"
		iqr="Interquartile range (IQR)";
	RUN; 
	title;
%MEND;

/* Doornbos test */
%MACRO Doornbos_test(dataset, var, id, alpha=0.05);
	ods select none;
	PROC MEANS data=&dataset mean var n;
		var &var;
		output out=out mean=mean var=var n=n; 
	RUN;
	
	DATA outliers;
		set &dataset(keep=&id &var);
		if _n_=1 then set out;
		/* leave-one-out variance */ 
		var_loo = ((n - 1) / (n - 2)) * var - (n / ((n - 1)*(n - 2)))*(&var - mean)**2;
		/* statistics */
		w = abs((&var - mean) / sqrt(var_loo * (n - 1) / n)); 
		/* critical value */
		c = quantile("t", 1 - &alpha / (2*n), n-2);
		/* check if this is an outlier */
		if(w > c) then outlier = "yes"; 
		else outlier = "no"; /* p-value */
		p_value = min(2*n*(1-cdf("t", w, n-2)), 1);
		keep &id p_value w c outlier; 
	RUN;
	
	PROC SORT data=outliers; 
		by descending w;
	RUN;
	ods select all;
	
	title1 "Doornbos test";
	title2 "Variable: &var";
	PROC PRINT data=outliers(obs=1) label noobs;
		var &id p_value w c outlier;
		label &id="id" p_value="p-value" 
		w="statistic (|W|)" c="critical value" 
		outlier="is outlier?";
	RUN;
	title; 
%MEND;

/* Hampel's Rule */
%MACRO Hampel(dataset, var, id);
	ods select none;
	PROC MEANS data=&dataset median;
		var &var;
     	output out=median median=median; 
    RUN;
    
    DATA hampel;
    	set &dataset(keep=&id &var); 
    	if _n_=1 then set median; 
    	abs_dev = abs(&var - median);
    RUN;
    
    PROC MEANS data=hampel median;
		var abs_dev;
		output out=abs_dev_median median=abs_dev_median;
	RUN;
	
	DATA hampel;
		set hampel;
		if _n_=1 then set abs_dev_median;
		abs_norm_val = abs_dev / abs_dev_median;
		if abs_norm_val <= 3.5 then delete;
	RUN;
	
	PROC SORT data=hampel;
		by descending abs_norm_val;
	RUN;
	ods select all;
    
    title1 "Outliers found using hampel’s rule.";
    title2 "Variable: &var";
	PROC PRINT data=hampel noobs label;
		var &id &var abs_norm_val;
      	label abs_norm_val="absolute normalized value (z_k)";
	RUN;
	title; 
%MEND;
      
      
      
      
      
      
      
      
      
      
      