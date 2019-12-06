%Macro SIM_PROP(nsim=, seed1=, alpha=);

data sample;
call streaminit(&seed1);
do nA= 10,25,50;
	do pA= 0.1, 0.3, 0.5;
		do pB= 0.1, 0.3, 0.5;
			do n=1 to &nsim by 1;
				do g=1 to 2 by 1;
					if(g=1) then 
						do i = 1 to nA by 1;
						X = rand("Binomial",pA,1);
	 					output;
						end;
					else	do i = 1 to nA by 1;
						X = rand("Binomial",pB,1);
	 					output;
						end;
				end;
			end;
		end;
	end;
end;
run; 

proc sort;
by nA pA nB pB;
run;

ods select none;
ods output ChiSq=powerC;
ods output FishersExact=powerF;
proc freq data=sample;
table g*X / chisq;
*exact chisq;
by nA pA pB n;
RUN;


data powerC;
set powerC;
where Statistic='Chi-Square';
rejectC = (Prob<&alpha);
keep nA pA pB n Prob rejectC;
*drop Table Statistic df;
run;

data powerF;
set powerF;
where Name1='XP2_FISH';
rejectF=(cValue1<&alpha);
keep nA pA pB n cValue1 rejectF;
*drop Table Name1 Label1 nValue1;
run;

data Final;
merge powerC powerF;
by nA pA pB n;
run;

ods select none;
proc means data=Final mean;
var rejectC rejectF;
by nA pA pB;
output out=power(drop= _TYPE_ _FREQ_);
run;

data power;
set power;
where _STAT_='MEAN';
drop _STAT_;
run;


%mend SIM_PROP;

%SIM_PROP(nsim=10000, seed1=1234,alpha=0.05);


