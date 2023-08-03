function [outgI,outgE,numCases] = TableSearch(gE,gI,T,t,A,a)

%find gE/gI values that minimize the norm ||(T,A)-(t-a)||_2 

%compute the norm to find the minimum distance with the pair (T,A)
NormAT=sqrt((T-t).^2+(A-a).^2);

%find the minimum
minNorm=min(min(NormAT));

%find the indexs of the minimum
[rowT,columnT]=find(NormAT==minNorm);

%count repetitions
numCases=length(rowT);

%find the values of gE/gI correspondig to the indexs
outgI=gI(rowT);
outgE=gE(columnT);


%if more than one pair is found, we select the firt one
if numCases>1
    outgI=outgI(1);
    outgE=outgE(1);
end
