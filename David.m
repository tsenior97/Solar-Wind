syms u(r)

y=1.221; %polytropic index
ro=696340; %solar radius

ue=400; %Speed at Earth 
pe=11.7*10^(-24); %density at Earth
BrE=5*y; %Magnetic field at Earth
TE=2*10^(5); %Temprature at Earth
PE=4.7*10^(18)*11.7*10^(-24); %Pressure near Earth

ra=24.3*ro; 
ua=3.32*10^(7);

rc=3*ro;
uc=1.82*10^(7);

rf=24.6*ro;
uf=3.33*10^(7);

F=9.02*10^(14)*1.05*10^(11);
Om=3*10^(-6); 
GM=1.33*10^(26); %gravitaional force times mass of the Sun

m=1.6735575*10^(-27); %hydrogen mass
k=1.38064852*10^(-23); %boltzmann constat

PA1=(2*k*TE/m)*(BrE^2/4*pi*pe*ue^2)^(y-1); %Pressure at critical point over density at critical point
MA=(u*r^2)/(ua*ra^2);

ode1=diff(u,r)==u/r*(((2*y*PA1/MA^(y-1))-GM/r)*(MA-1)^(3)+Om^2*r^2*(u/ua-1)*((MA+1)*(u/ua)-3*MA+1))*((u^2-(y*PA1/MA^(y-1)))*(MA-1)^3-Om^2*r^2*MA*(ra^2/r^2-1)^(2))^(-1);
cond1=u(ra)==ua;
cond2=u(rc)==uc;
cond3=u(rf)==uf;
conds=[cond1 cond2 cond3];
opts = odeset('RelTol',1e-3);
usol(r)=dsolve(ode1,opts);



