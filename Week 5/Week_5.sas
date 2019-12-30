LIBNAME SASDATA "/folders/myfolders/statistical-methods/data";

/* Question 5.1 */
/* you assumed that the age of the mother (AGEM) */
/* is normally distributed. In this question you */
/* will verify these assumptions. */

DATA WEEK5_Q1;
	set SASDATA.IVF;
RUN;

/* a) */
/* TODO: Doesn't generate same results as answer */

/* Shapiro-Wilk */
/* Test statistic = 0.992907 and p-value = 0.0011 */
/* Conclusion: */
/* Is sensitive to non symmetric departures from normality. */
/* This distribution seems relatively symmetri. */

/* b) */
/* Kolmogorov-Smirnov */
/* Test statistic = 0.040049 and p-value < 0.0100 */
/* Conclusion: */
/* As it compares the vertical distance (supremum) between */
/* the distributions we can see that there is a relatively */
/* large distance between F_n and F. */

/* c) */
/* Cramer-von Mises */
/* Test statistic = 0.178199 and p-value = 0.0099 */
/* Conclusion: */
/* Puts more weight on the tails as the tails deviate less */
/* it gets a higher score then the other tests as the tails */
/* Seem relatively normal. */

/* d) */
/* Anderson-Darling */
/* Test statistic = 1.271699 and p-value < 0.0050 */
/* Conclusion: */
/* As it looks at the integral of the distance between F_n */
/* and F. This has the smallest p-value of all the tests. */

/* e) */
/* TODO */

/* f) and g) */
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

PROC UNIVARIATE data=WEEK5_Q1;
	var AGEM;
RUN;
/* Skewness = -0.1773114, Kurtosis = -0.1983322 */
/* and n = 759 */

%Skewness_Kurtosis_test(skewness=-0.1773114, kurtosis=-0.198322, n=759);




