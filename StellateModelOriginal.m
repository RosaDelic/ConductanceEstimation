function dx=StellateModelOriginal(t,x,g,param)

%Neural parameters
Iapp=param(1); 
gL=param(2); %0.1; 
vL=param(3); %-65;
C=param(4);  %1;
vE=param(5); %0;
vI=param(6); %-80;
gK=11.0; gNa=52.0; gp=0.5; gh=1.5; 
vK=-90.0; vNa=55.0; vh=-20;
gE=g(1); gI=g(2);

%Initial conditions
v=x(1); m=x(2); h=x(3); n=x(4); p=x(5); rf=x(6); rs=x(7);

%Currents:

%leak current
%il=-gL*(v-vL);

%sodium current
%ina=-gNa*m^3*h*(v-vNa);
am=-0.1*(v+23)/(exp(-0.1*(v+23))-1);
bm=4*exp(-(v+48)/18);
minf=am/(am+bm);
taum=1/(am+bm);

ah=0.07*exp(-(v+37)/20);
bh=1/(1+exp(-0.1*(v+7)));
hinf=ah/(ah+bh);
tauh=1/(ah+bh);

%delayed rectifier -- Potassium current
%ik=-gK*(n^2*n^2)*(v-vK);       
an=-0.01*(v+27)/(exp(-0.1*(v+27))-1);
bn=0.125*exp(-(v+37)/80);
ninf=an/(an+bn);
taun=1/(an+bn);

%persistent sodium current (I_Nap)
%iNap=-gp*p*(v-vNa);      
ap=1/(0.15*(1+exp(-(v+38)/6.5)));
bp=exp(-(v+38)/6.5)/(0.15*(1+exp(-(v+38)/6.5)));
pinf=ap/(ap+bp);
taup=1/(ap+bp);

%hyperpolarization-activated current (fast component)
%ih=-gh*(0.65*rf+0.35*rs)*(v-vh);
rfinf= 1/(1+exp((v+79.2)/9.78));
taurf=0.51/(exp((v-1.7)/10) + exp(-(v+340)/52))+1;
rsinf= 1/(1+exp((v+2.83)/15.9))^58;
taurs=5.6/(exp((v-1.7)/14) + exp(-(v+260)/43))+1;

%Gating variables
gionic=gNa*m^3*h+gK*n^2*n^2+gp*p+gh*(0.65*rf+0.35*rs);
gsyn=gL+gE+gI;
gT=gionic+gsyn;
Iion=gNa*m^3*h*vNa+gK*n^2*n^2*vK+gp*p*vNa+gh*(0.65*rf+0.35*rs)*vh;
Isyn=gL*vL+vE*gE+vI*gI;
iT=Iion+Isyn+Iapp;

dx(1)=1/C*(-gT*v+iT);        %deriv(Voltage)
dx(2)=(minf-m)/taum;         %deriv(m)
dx(3)=(hinf-h)/tauh;         %deriv(h)
dx(4)=(ninf-n)/taun;         %deriv(n)
dx(5)=(pinf-p)/taup;         %deriv(p)
dx(6)=(rfinf-rf)/taurf;      %deriv(rf)
dx(7)=(rsinf-rs)/taurs;      %deriv(rs)

