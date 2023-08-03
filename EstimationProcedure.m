function [EstimatedgE, EstimatedgI] = EstimationProcedure(A,T,gE,gI,V, t)
%calculate peaks and indexs of the voltage V
[pks,locs]=findpeaks(V);

%eliminate non-spikes
k=find(pks<20);
pks(k)=[];
locs(k)=[];
if length(pks)<2
    error('The membrane potential does not present a spiking regime')
end

%search amplitudes and periods
Tinterspike=zeros(1,length(locs)-1);
Aspike=zeros(1,length(locs));

for i=1:(length(locs)-1)
    vmin=min(V(locs(i):locs(i+1)));
    Aspike(i)=abs(pks(i)-vmin);
    Tinterspike(i)=abs(t(locs(i+1))-t(locs(i)));
end

Aspike(length(locs))=abs(pks(end)-min(V(locs(end):end)));
Aspike=Aspike(2:end);

%initialize vectors for TableSearch
gEtable=zeros(1,length(Aspike));
gItable=zeros(1,length(Aspike));
gEnewton=zeros(1,length(Aspike));
gInewton=zeros(1,length(Aspike));
numCases=zeros(1,length(Aspike));

%estimate discretized conductances 
for i=1:length(Aspike)
    [gItable(i),gEtable(i),numCases(i)] = TableSearch(gE,gI,T,Tinterspike(i),A,Aspike(i));
    [gInewton(i),gEnewton(i)] = InterpNewton(gE,gI,T,Tinterspike(i),A,Aspike(i),gEtable(i),gItable(i),2000,0.000001);
end

%define the time of the estimated conductances
tstr=zeros(1,length(locs));
for i=1:length(locs)
    tstr(i)=t(locs(i));
end
tstr=tstr(2:end);


%interpolate estimated conductances
EstimatedgE=spline(tstr,gEnewton,t);
EstimatedgI=spline(tstr,gInewton,t);

%plot estimated conductances
lleva=20;
lleva2=60;
TOL = 1e-2;

figure(1);
hold on; 
set(gca,'FontSize',24);
plot(t(floor(lleva/TOL):end-floor(lleva2/TOL)),EstimatedgE(floor(lleva/TOL):end-floor(lleva2/TOL)),'-k','DisplayName','g_{E,Estimated}','LineWidth',2);
plot(t(floor(lleva/TOL):end-floor(lleva2/TOL)),EstimatedgI(floor(lleva/TOL):end-floor(lleva2/TOL)),'-b','DisplayName','g_{I,Estimated}','LineWidth',2);
ylabel('Conductances (µS/cm^2)','FontSize',24.4);
xlabel('t(ms)');
hold off;
lgd=legend();
set(lgd,'Orientation','horizontal','FontSize',10,'Location','northoutside');
end