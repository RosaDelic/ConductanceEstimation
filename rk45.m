function [wi, ti] = rk45 (RHS, t0, x0, tf, N ,gEv, gIv, param)
%Solve Runge-Kutta45
neqn = length(x0);
ti(1) = t0;
wi(1:neqn, 1) = x0';
i = 2;

if(N <= 0)
    disp( 'N must be positive and different to 0' );
    return;
else
    h = (tf - t0)/N;
end
while(t0 < tf)
    %set the actual gE and gI values
    gE=gEv(i);
    gI=gIv(i);
    
    g=[gE,gI];
    k1 = h * feval(RHS, t0, x0, g, param);
    k2 = h * feval(RHS, t0 + h/2, x0 + k1/2, g, param);
    k3 = h * feval(RHS, t0 + h/2, x0 + k2/2, g, param);
    k4 = h * feval(RHS, t0 + h, x0 + k3, g, param);
    
    fv=k1(1)/h;
%   fprintf(fdesti,'%d \t %d \n' ,t0,fv);
    
    x0 = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;
    t0 = t0 + h;
    ti(i) = t0;
    wi(1:neqn, i) = x0';
	i = i + 1;
    
end
    
    