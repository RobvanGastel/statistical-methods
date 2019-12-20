LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 3.1 */
%Macro SIM_Gum(alpha=, nsim=, seed=); 
	PROC IML;
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
	QUIT;
%mend SIM_Gum;

%SIM_Gum(nsim=1000, alpha=5, seed=12345);

/* a) */
PROC CORR data=GUMC kendall; 
	var X Y;
RUN;

PROC CORR data=GUMC spearman; 
	var X Y;
RUN;

/* b) */


/* c) */


/* Question 3.2 */
DATA WEEK3_Q2;
	set SASDATA.IVF;
	IMP = IMP + (ranuni(1)-0.5); 
RUN;

/* a) */
PROC TRANSPOSE out=WIDE_IVF(drop = _NAME_ _LABEL_) 
		data=WEEK3_Q2 prefix=IMP;
	by ID;
	id PER;
	var IMP; 
RUN;

/* b) */
ods graphics on;
PROC CORR data=WIDE_IVF plots=scatter(ellipse=none); 
	var IMP10 IMP18;
RUN;
ods graphics off;

PROC CORR data=WIDE_IVF kendall; 
	var IMP10 IMP18;
RUN;

PROC CORR data=WIDE_IVF spearman; 
	var IMP10 IMP18;
RUN;

/* c) */

/* Question 3.3 */

/* a) */


