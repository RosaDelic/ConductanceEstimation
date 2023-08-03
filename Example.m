%initialize parameters
clear all
clc

format long;
lleva=20;
lleva2=60;

TOL = 1e-2;
t0 = 0;
tf = 750;
N=(tf-t0)/TOL;

C=1; Iapp=8;
gL=0.5; vL=-65;
vE=0; vI=-80;
param=[Iapp,gL,vL,C,vE,vI];
x0=[-65.0, 0.1, 0.9, 0.1, 2.0, 0.4, 0.6];

%length of conductance vectors
long=N+2;

%to find the maximum and minimum
llevar=25000;

%load gE,gI,A,T
load('Results.mat')
%% ---------------------------CHOOSE AN OPTION-----------------------------
%variable to indicate the end of the simulation
final=0;

while ~final
    clc
    disp('Indicate which example would you like to illustrate:')
    disp('Type 0 for constant conductances example');
    disp('Type 1 for simple frequency conductances example');
    disp('Type 2 for double frequency conductances example');
    disp('Type 3 for in silico conductances example');
    option=input('Case: ');





    %% --------DEFINE THE gE/gI CORRESPONDING TO THE CHOSEN OPTION --------
    if option==0
        %CONSTANT CASE

        %Introduce a constant gE
        errorgE=1;
        while errorgE
            disp('Choose constant gE value in interval (0.01,0.20)')
            cntgE=input('gE: ');
            if cntgE<=0.01 || cntgE>=0.20
                disp('You have introduced a wrong gE value');
            else
                errorgE=0;
            end
        end
        
        %Introduce a constant gI
        errorgI=1;
        while errorgI
            disp('Choose constant gI value in interval (0.24,0.43)')
            cntgI=input('gI: ');
            if cntgI<=0.24 || cntgI>=0.43
                disp('You have introduced a wrong gI value');
            else
                errorgI=0;
            end
        end

        %define constant conductances (gE and gI) vectors
        vectorgE=cntgE*ones(long,1);
        vectorgE=transpose(vectorgE);
        vectorgI=cntgI*ones(long,1);
        vectorgI=transpose(vectorgI);
        
        %solve the system and estimate the conductances
        [x, t]=rk45('StellateModelOriginal', t0, x0, tf, N , vectorgE, vectorgI, param);
        [estimatedgE, estimatedgI] = EstimationProcedure(A,T,gE,gI,x(1,:), t);

    elseif option==1 || option==2
        %SINUSOIDAL CASE

        %inicialize gE,gI,tE vectors
        vectorgE=ones(long,1);
        vectorgE=transpose(vectorgE);
        vectorgI=ones(long,1);
        vectorgI=transpose(vectorgI);
        tE=0:TOL:(tf+TOL);

        %choose period and option
        T0=300; %period(ms)---->slow frequency
        T1=75;  %period(ms)---->fast frequency
        %A=0.095 maximum allowable amplitude
    
        %introduce the phase shift between the excitatory and inhibitory conductances 
        desfase=input('Introduce the phase shift between the excitatory conductance and the inhibitory one: ');
        
        %define the conductance (gE and gI) vectors
        if option==1
            %SIMPLE FREQUENCY
            vectorgE=0.055*cos(2*pi*tE/T0)+0.105;
            vectorgI=0.075*cos(2*pi*(tE/T0)+desfase)+0.335;
        else
            %DOUBLE FREQUENCY
            vectorgE=0.060*cos(2*pi*tE/T0)+0.035*cos(2*pi*tE/T1)+0.105;
            vectorgI=0.060*cos(2*pi*(tE/T0)+desfase)+0.035*cos(2*pi*(tE/T1)+desfase)+0.335;
        end

        %solve the system and estimate the conductances
        [x, t]=rk45('StellateModelOriginal', t0, x0, tf, N , vectorgE, vectorgI, param);
        [estimatedgE, estimatedgI] = EstimationProcedure(A,T,gE,gI,x(1,:), t);


    elseif option==3
        %IN SILICO CONDUCTANCES CASE

        %read the conductances (gE and gI) from file

        %open input file
        finput = fopen('sd93-01.dat', 'r');
        %dimensions of input table
        column=3;
        row=10000;
        %read matrix of values
        G = fscanf(finput, '%g %g', [column row]);
        %close the file
        fclose(finput);
        At=transpose(G);
        %rename the columns of matrix A
        vectorgE=(At(:,1)+At(:,3))*0.001;
        vectorgI=(At(:,2))*0.001;

        %define the conductance (gE and gI) vectors
        tE=t0:1:tf+TOL;
        vectorgE=vectorgE(1:length(tE));
        vectorgI=vectorgI(1:length(tE));
        vectorgE=spline(tE,vectorgE);
        vectorgI=spline(tE,vectorgI);
        tE=t0:TOL:tf+TOL;
        vectorgE= transpose(ppval(vectorgE,tE));
        vectorgI= transpose(ppval(vectorgI,tE));

        %solve the system and estimate the conductances
        [x, t]=rk45('StellateModelOriginal', t0, x0, tf, N , vectorgE, vectorgI, param);
        [estimatedgE, estimatedgI] = EstimationProcedure(A,T,gE,gI,x(1,:), t);

    else
        disp('You have introduced a wrong number');
    end

    %plot actual conductances
    figure(1)
    hold on;
    set(gca,'FontSize',24);
    plot(t(floor(lleva/TOL):end-floor(lleva2/TOL)),vectorgE(floor(lleva/TOL):end-floor(lleva2/TOL)),'--k','DisplayName','g_{E,Actual}','LineWidth',2);
    plot(t(floor(lleva/TOL):end-floor(lleva2/TOL)),vectorgI(floor(lleva/TOL):end-floor(lleva2/TOL)),'--b','DisplayName','g_{I,Actual}','LineWidth',2);
    ylabel('Conductances (µS/cm^2)','FontSize',24.4);
    xlabel('t(ms)');
    lgd=legend();
    set(lgd,'Orientation','horizontal','FontSize',10,'Location','northoutside');
    hold off;
    
    %allow the user to choose if he/she wants to continue with another example
    error=1;
    while error
        disp('Do you want to illustrate another example?');
        disp('YES-->Type 0');
        disp('NO-->Type 1');
        final=input('Choose one option: ');
        if final==0 || final==1
            error=0;
        else
            disp('You have introduced a wrong option')
        end
    end
    if final==0
        close all
    end
end
close all