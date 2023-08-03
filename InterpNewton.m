function [outgI,outgE] = InterpNewton(gE,gI,T,t,A,a,gE0,gI0,nmax,tol)

%-------------------INTERPOLATE FUNCTION AND DIFFERENTIAL------------------
%initial conditions
gEi=gE0;
gIi=gI0;
 
%solve sistem:
% f=t <=> f-t=0
% g=a <=> g-a=0
T=T-t;
A=A-a;

%calculate the discretization of the gradient
dE=abs(gE(2)-gE(1));
dI=abs(gI(2)-gI(1));

%calculate partial derivatives of T and A with respect to gE and gI
[Te,Ti]=gradient(T,dE,dI);
[Ae,Ai]=gradient(A,dE,dI);

%create a mesh for interpolation
[gIgrid,gEgrid]=ndgrid(gI,gE);

%interpolated functions of T and A: 
% T->f1 // A->f2
f1=griddedInterpolant(gIgrid,gEgrid,T,'spline');
f2=griddedInterpolant(gIgrid,gEgrid,A,'spline');

%interpolated functions of partial derivatives
% Te->d11 // Ti->d12 // Ae->d21 // Ai->d22 
d11=griddedInterpolant(gIgrid,gEgrid,Ti,'spline');
d12=griddedInterpolant(gIgrid,gEgrid,Te,'spline');
d21=griddedInterpolant(gIgrid,gEgrid,Ai,'spline');
d22=griddedInterpolant(gIgrid,gEgrid,Ae,'spline');

%---------------------APPLY CLASSIC NEWTON'S METHOD------------------------


%evaluate function(F) and differential(D) at x0 
F=[f1(gI0,gE0); f2(gI0,gE0)];
D=[d11(gI0,gE0) d12(gI0,gE0); d21(gI0,gE0) d22(gI0,gE0)];

%first step: solve the system 
dk=linsolve(D,-F);

%calculate next x_k
gI1=gI0+dk(1);
gE1=gE0+dk(2);

%update number of iterations
n=2;

%while non-stop conditions are satisfied, iterate to find a solution
while n<nmax && (abs(gI1-gI0)>=tol || abs(gE1-gE0)>=tol)
    %update x_k
    gI0=gI1;
    gE0=gE1;
    
    %evaluate at next x_k
    F=[f1(gI0,gE0); f2(gI0,gE0)];
    D=[d11(gI0,gE0) d12(gI0,gE0); d21(gI0,gE0) d22(gI0,gE0)];
    
    %if the system does not have a solution, finish
    if det(D)==0
        outgI=gIi;
        outgE=gEi;
        break;
    else
        %solve the system
        dk=linsolve(D,-F);
        
        %calculate next x_k
        gI1=gI0+dk(1);
        gE1=gE0+dk(2);
        
        %next iteration
        n=n+1;
    end
end

%if the system does not converge, we preserve initial gE/gI values
if n>=nmax || gI1<=0.23 || gI1>=0.44 || gE1<=0 || gE1>=0.21
    outgI=gIi;
    outgE=gEi;
else
    outgI=gI1;
    outgE=gE1;
end

