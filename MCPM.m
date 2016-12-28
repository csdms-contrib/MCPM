clear;close all;clc;

%basic parameters
rho=1030;%water density
L=400;
dx=0.2;dx2=dx^2;%cell discretization [m]
N=75/dx;%transect length
x=[0:N-1]*dx;

%flow
kr=0.2;kro=10^-8; %kr is the bed elevation variability [m]. kro is for numerical stability
Cb=0.005; %bed drag coefficient
ws=0.5/1000;%in the paper I use this 0.5/1000;% settlign velocity [m/s]
me=0.0001; %sediment erodability [kg/m2/s]
taucr=0.1; % critical shear stress [Pa]
K=3.65/365/24; %m2/hr creeping 
Kveg=0.5/365/24;
slr=2.6/1000/365/24; %relative sea level rise [m/hr]
Co=20/1000;% reference seidment concentration in the flood water [g/l] = [kg/m3]

rmud=2500; %kg/m3 %for the mud (dry bulk density with 40% mud 60 water. already very consolidated)
Wmud=0.4;%rmud=1000; %kg/m3 %for the mud (dry bulk density with 40% mud 60 water. already good consolidated)
rorg=1200;%density of organic matter * 1 - water content
Worg=0.1;

%tide
amp=2.7/2;%mean tidal range of Plum Island south. 
T=12.5;%tidal period, hours

%time paramters
mor=1;
dt=5/60;%hours
tint=floor(24*365/dt);% [hours]. Dislay interval, does not affect simulation
numyr=10;tmax=numyr*24*365;
time=[dt:dt:tmax];
y=amp*sin(2*pi*time/T);%water level. Can be any function, can be asymmetric
dydti=[y(3:end)-y(1:end-2)]/(2*dt)/3600;dydti=[dydti(1) dydti dydti(end)];%m/s %dydti=2*pi/T*amp*cos(2*pi*t/T)/3600;

%marsh parameters
BMax=2.5; %kg/m2 % (to g/cm3) 1.5 from morris ocenaoraphy
nuGp=0.0138;%1/day%rate at which organic matter is stored. to convert to AMC, total mass of organic accumalted
chiref=0.158;%refractory fraction
Dmin=0;Dmax=amp+0.237*(amp*2)-0.092;
AA=0.25*(Dmin+Dmax)*(Dmax-3*Dmin);

%specify which boundary type.
%case '1'%. Ci=0 conserve mass
%case '2';% Ci=C for ebb and flood. no variation along axis channel
%case '3';% Ci=Cstm (during flood). Ci=C (during ebb)
boundary='3';

%initial topography and mass
h=(amp*0.9)*ones(N,1);h(1:1/dx)=-amp;
M=zeros(N,1);%suspended sediment mass

%initialization
A1=zeros(N);A2=eye(N);Aeye=find(A2>0);A3=zeros(N);
Ctide=NaN*time;Tau=NaN*time;Q=NaN*time;qx1=NaN*time;Ux1=NaN*time;
ER=h*0;DEP=h*0;CR=h*0;Fcr=h*0;
hh=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 

for t=1:length(time)
d=y(t)-h;

eta=0.5*(1+erf(2*d/kr));%fraction of wetted area
d=eta.*d+kr/4/sqrt(pi)*exp(-4*(d/kr).^2);d(d<kro)=kro;

%SALT MARSH GROWTH
dm=amp-h;%marsh depth relative to MHW
Bpeak=BMax*(Dmax-dm).*(dm-Dmin)/AA;Bpeak(Bpeak<=0)=0;
org=(Bpeak*365/2*nuGp)*chiref/(rorg*Worg)/(365*24); %in m/hour

%vegetation characteristics as function of biomass
hs=0.0609*(1000*Bpeak).^0.1876;
ns=250*(1000*Bpeak).^0.3032;
ds=0.0006*(1000*Bpeak).^0.3;
as=0.25*(1000*Bpeak).^0.5;
Cdo=1;Cv=1/2*as.*d.*Cdo;%vegetation drag coefficent

%FLOW
[Ux,qx,Uy,qy,Uref,tau]=flow(N,Cb,Cv,L,dx,d,eta,dydti(t),y,Bpeak,rho,kr);
if y(t)<=(h(end)-kr);tau=tau*0;end

%vegeation induced sedimentation
neta=0.224*(Uref.*ds/10^-6).^0.718.*(100*10^-6./ds).^2.08;
wsB=Uref.*neta.*ds.*ns.*min(hs,d);wsB(Bpeak<=0)=0;
WS=ws+wsB;
WS=WS*2;%double the settlign velocity (same thing as doubling the sediment concentration near the bed, used to parametrize SSC profile)

%EROSION
er=3600*me*max(0,tau-taucr)/taucr /(Wmud*rmud);
M=M+er*(Wmud*rmud)*dt;

%ADVECTION_DIFFUSION_SEDIMENTATION

%Diffusion
Diff=0.13*d.*sqrt(taucr/rho).*eta;
dsR=([Diff(2:end); Diff(end)]+Diff)/2.*([d(2:end); d(end)]+d)/2.*(d>0).*([d(2:end); d(end)]>0);
dsL=([Diff(1); Diff(1:end-1)]+Diff)/2.*([d(1); d(1:end-1)]+d)/2.*(d>0).*([d(1); d(1:end-1)]>0);

A3(Aeye(2:end-1))=3600*dt/dx2*((dsR(2:end-1)+dsL(2:end-1))./d(2:end-1));
A3(1,1)=3600*dt/dx2*(dsR(1)./d(1));
A3(N,N)=3600*dt/dx2*(dsL(N)./d(N)); 
A3(Aeye(1:end-1)+N)=-3600*dt/dx2*(dsR(1:end-1)./d(2:end));
A3(Aeye(2:end)-N)=-3600*dt/dx2*(dsL(2:end)./d(1:end-1));

%Advection
if dydti(t)>0 %(Uy>0) (flood, to the right) 
A1(Aeye(end))=0;A1(Aeye(1:end-1)+N)=0;
A1(Aeye(1:end-1))=Uy(1:end-1)/dx*dt*3600;
A1(Aeye(2:end)-N)=-Uy(1:end-1)/dx*dt*3600;
else %Uy(1)<0; (ebb, to the left)
A1(Aeye(1))=0;A1(Aeye(2:end)-N)=0;
A1(Aeye(2:end))=-Uy(2:end)/dx*dt*3600;
A1(Aeye(1:end-1)+N)=Uy(2:end)/dx*dt*3600;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch boundary
  case '1'%. Ci=0 conserve mass
    A2(Aeye)=(1 +dt*3600./d.*WS.*eta);

  case '2';% Ci=C for ebb and flood. no variation along axis channel.
    A2(Aeye)=(1 +dt*3600./d*ws.*eta -qx/L*dt*3600./d);

  case '3';% Ci=Cstm (during flood). Ci=C (during ebb)
    if dydti(t)>0 %(flood, to the right) 
     M=M+Co*qx/L*dt*3600;
     A2(Aeye)=(1 +dt*3600./d.*WS.*eta);% );
    else % (ebb, to the left) 
     A2(Aeye)=(1 +dt*3600./d.*WS.*eta -qx/L*dt*3600./d);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%update amount of suspended sediment
if d(1)>kr
MM=(A1+A2+A3)\M;MM(isnan(MM))=0;
dep=3600*WS.*MM./d/(Wmud*rmud).*eta;
else;MM=M;MM(isnan(MM))=0;dep=0*h;end

%redistributed sediment of bed that got lost
a=find(d<=kro);dep(a)=dep(a)+MM(a)/(Wmud*rmud)/dt;MM(a)=0;M=MM;

C=M./d;C(y(t)<=h(1))=NaN;

%CREEP
FCR=K*ones(N,1);FCR(Bpeak>0)=Kveg;
FR=(h-[h(2:end); h(end)]).*FCR;FL=([h(1); h(1:end-1)]-h).*[FCR(1); FCR(1:end-1)];
F=-(FR+FL)/2/dx;
cr=(FL-FR)/dx^2;

%UPDATE BED
h=h+(org+dep-er+cr-slr)*dt*mor;


ER=ER+er;DEP=DEP+dep;CR=CR+cr;Fcr=Fcr+F;%This is just for graphic purpose
if mod(t,tint)==0;
ER=ER/tint;DEP=DEP/tint;CR=CR/tint;Fcr=Fcr/tint;%This is just for graphic purpose
org(org==0)=NaN;%This is just for graphic purpose
xe=x(end)/3;%This is just for graphic purpose

subplot(4,1,1)
plot(x,h,'k',x,ones(N,1)*(amp-(Dmax-Dmin)/2),'--r',x,ones(N,1)*(amp-Dmax),'--r',x,ones(N,1)*(amp),'--r',x,ones(N,1)*(-amp),'--r')
ylabel('d');axis([x(1) xe -4.7 amp+0.1])
title(mor*time(t)/24/365)

subplot(4,1,2)
s=1000*365*24;plot(x,ER*s,'r',x,DEP*s,'m',x,-CR*s,'c',x,org*s,'g',x,ones(N,1)*slr*s,'b');
xlim([x(1) xe])

subplot(4,1,3)
s=1000*365*24;plot(x,ER*s,'r',x,DEP*s,'m',x,-CR*s,'c',x,org*s,'g',x,ones(N,1)*slr*s,'b');
xlim([x(1) xe]);ylim([-10 20])

subplot(4,1,4)
plot(x,Fcr*365*24,'-')
xlim([0 xe])

pause(0.001)
end
end