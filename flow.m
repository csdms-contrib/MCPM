function [Ux,qx,Uy,qy,Uref,tau]=flow(N,Cb,Cv,L,dx,d,eta,dht,y,Bpeak,rho,kr);
g=9.8;
if abs(dht)<10^-10;Ux=0*d;Uy=0*d;qx=0*d;qy=0*d;Uref=0*d;tau=0*d;else
    
qt=dht*L*eta;
Q=sum(qt);

Cd=Cb+Cv;
Ux=sqrt(d)./sqrt(Cd).*eta;%this is an advection velocity, it already includes eta
Ux=Ux*Q/sum(Ux.*d);

tau=rho*Cb.*Ux.^2;

qx=Ux.*d;% discharge per unit of width
qy=cumsum(-qt+qx)/L*dx;
Uy=qy./d;

Uref=sqrt(Ux+Uy)/2;Uref(Uref<=0)=0;
end


