function [gE,gI,A,T] = TablesCreator(model,d)
    %--------------------------INITIALIZE PARAMETERS---------------------------
    format long;
    
    TOL = 1e-2;
    t0 = 0;
    tf = 750;
    N=(tf-t0)/TOL;
    
    C=1; Iapp=8;
    gL=0.5; vL=-65;
    vE=0; vI=-80;
    param=[Iapp,gL,vL,C,vE,vI];
    x0=[-65.0, 0.1, 0.9, 0.1, 2.0, 0.4, 0.6];
    %--------------------------------------------------------------------------
    %Length of conductance vectors gE, gI
    long=N+2; 
    
    %Initial and final conductance values (gE,gI)
    xini=0.01;
    xfin=0.20;
    yini=0.24;
    yfin=0.43;
    
    %Total number of conductance values that we take
    n=round((xfin-xini)/d);
    
    %Useful variable to treat the conductances
    llevar=25000;
    
    %Create meshgrid for the conductances vectors
    gE=linspace(xini,xfin,n); 
    gI=linspace(yini,yfin,n); 
    
    
    %for each pair of conductances
    for j=1:length(gE)
      for i=1:length(gI)
        
        %Define constant gE, gI vectors
        gEe=gE(j)*ones(long,1);
        gEe=transpose(gEe);
        gIi=gI(i)*ones(long,1);
        gIi=transpose(gIi);
        
        %Execute Runge-kutta 45 to solve the neuronal model
        [x, t]=rk45(model, t0, x0, tf, N , gEe, gIi, param);
        
        %Compute maximum and minimum of the solution v(t)
        maxim=max(x(1,end-llevar:end));
        minim=min(x(1,end-llevar:end));
        
        %Compute the amplitud
        dif=abs(maxim-minim);
           
        %Check if the pair of conductances (gE,gI) produces spikes
        if dif<=80
          %if it does not produce spikes, assign zero value to the amplitude and period
          A(i,j)=0;
          T(i,j)=0;
        else
          %if it produces spikes, assign the corresponding value to the amplitude
          A(i,j)=dif;
          
          %Calculate peaks and indexs of v(t)
          [pks,locs]=findpeaks(x(1,:));
          
          %Eliminate peaks that are not high enough
          k=find(pks<20);
          pks(k)=[];
          locs(k)=[];
          
          %Calculate the period
          periode=t(locs(end))-t(locs(end-1));
          
          T(i,j)=periode;
        end
      end
    end

    save('Results.mat','gE','gI','A','T');
end