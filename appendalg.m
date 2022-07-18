%function[p,q,s,t] = transalgor()
%function[p,q,s,t] = trans1( p0 , q0, s0 ,Ca,Cb,Tend)

t1=tic;
Db=18; %diffusion coefficient in base metal (18*10^-6 m^2 s^-1)

Da=500; %*10^(-4); %diffusion coefficient in the interlayer (500*10^-6 m^2 s^-1)

%for all I at j=0 p0=0.0 this is the initial condition within the domain of the base metal

%for all I at j=0 q0=19.0 this is the initial condition within the domain of the interlayer

%for all i=n-1   Ca=0.166 this is the boundary condition at the side of the base metal

%for all i=0  Cb=10.23  this is the boundary condition at the side of the interlayer

p01=19.0;%0.0;
q01=0.0;%19.0; 

Cb=0.166;%10.23;
Ca=10.23;%0.166;


l2=25; %0.000025*0.5; %this is the thickness of the interlayer

l1=6025; %0.006025; %this is the length of the base metal
%Tend=100000000000;
Tend=0.5;
dt=0.0001;

t=0:dt:Tend;
t(1)=0;

nt=length(t);

nu=10;
%du=l1/(2*(nu-1));
%du=1/((nu-1));
l=(0.5*(l1));%+l2));
%p0=zeros(1,nu);

p=zeros(1,nu);
%%p(2:nt,nu)=Ca;

p(1,:)=p01;

piplusend=p01;
pjsig(1,1:nu)=p01;

nv=10;
%dv=1/((nv-1));
q=zeros(1,nv);
%%q(2:nt,1)=Cb;

q(1,:)=q01;

qiplusbeg=q01;

qjsig(1,1:nv)=q01;

u=[0 0.3010 0.4771 0.6021 0.6990 0.7782 0.8451 0.9031 0.9542 1.0000];
v=[0 0.0458 0.097 0.155 0.2219 0.3011 0.398 0.5229 0.699 1];

%{
%end
u=zeros(1,nu); v=zeros(1,nv);
%estimating the length at each coordinate
  %check this later
for i=1:nu
u(i)=(i-1)*du;
end
for i=1:nv
v(i)=(i-1)*dv;
end

%}


s=zeros(1,nt);
%sjsig0=zeros(1,nt);


s(1)=l2*0.5;
sjsig=s(1);
%t(1)=0;
s1=zeros(1,nt);
j=1;
rgnu=1;
while 1
    if t(1)<0.00000001
    t(j+1)=(j)*dt;
    elseif t(1)>0.0001
        t(j+1)=t(j)+dt;
    end
    
    
while 1

    s(j+1)=((((Db*dt)/(l-sjsig))*((qjsig(2)-Cb)/(v(2)))-((Da*dt)/(sjsig))*((Ca-pjsig(nu-1))/(1-u(nu-1))))/...
    ((0.5*(1+u(nu-1))*piplusend)+(0.5*(1-u(nu-1))*Ca)-((1-0.5*v(2))*qiplusbeg)-((0.5*v(2))*Cb)))+s(j);


sjsig=s(j+1);

%sjsig0(j+1)=s(j+1);%usefulness still not visible

    nu1=nu-1; nv1=nv-1; %number of unknown is less than total segment because of the known boundary values
a=zeros(nu1,nu1);
b=zeros(nv1,nv1);
L11=zeros(1,nu1);M11=zeros(1,nu1);R11=zeros(1,nu1);%,W11,Ohm11
L12=zeros(1,nv1);M12=zeros(1,nv1);R12=zeros(1,nv1);%,W11,Ohm11

if s(j+1)>s(j)

%case a: s(j+1)>s(j)
for i=2:nu1
L11(i)=((s(j+1)*(u(i+1)-u(i))+((dt*Da)/(s(j+1)*(u(i+1)-u(i))))+((dt*Da)/(s(j+1)*(u(i)-u(i-1))))+((s(j+1)-s(j))*u(i))))/...
    ((dt*Da)/(s(j+1)*(u(i)-u(i-1))));
M11(i)=(((dt*Da)/(s(j+1)*(u(i+1)-u(i))))+(s(j+1)-s(j))*u(i))/((dt*Da)/(s(j+1)*(u(i)-u(i-1))));
R11(i)=(s(j)*(u(i+1)-u(i)))/((dt*Da)/(s(j+1)*(u(i)-u(i-1))));
W11=(((dt*Da)/(s(j+1)*(u(2))))+(s(j+1)-s(j))*u(2))/((dt*Da)/(s(j+1)*(u(2)))+(s(j+1)*(u(2))));
Ohm11=((s(j)*u(2))/((dt*Da)/(s(j+1)*(u(2)))+(s(j+1)*(u(2)))));
end

for i=2:nv1
L12(i)=((l-s(j+1))*(v(i+1)-v(i))+((dt*Db)/((l-s(j+1))*(v(i+1)-v(i))))+((dt*Db)/((l-s(j+1))*(v(i)-v(i-1))))+...
    (s(j+1)-s(j))*(1-v(i)))/((dt*Db)/((l-s(j+1))*(v(i)-v(i-1))));
M12(i)=(((dt*Db)/((l-s(j+1))*(v(i+1)-v(i))))+(s(j+1)-s(j))*(1-v(i)))/((dt*Db)/((l-s(j+1))*(v(i)-v(i-1))));
R12(i)=((l-s(j))*(v(i+1)-v(i))/((dt*Db)/((l-s(j+1))*(v(i)-v(i-1)))));

W12=((l-s(j+1))*(-v(nv))+((dt*Db)/((l-s(j+1))*(v(nv)-v(nv-1))))+(s(j+1)-s(j))*(1-v(nv)))/...
    ((dt*Db)/((l-s(j+1))*(v(nv)-v(nv-1))));

Ohm12=(l-s(j+1))*(-v(nv))/((dt*Db)/((l-s(j+1))*(v(nv)-v(nv-1))));
end

a(1,1:2)=[1 -W11];
a(length(a),(length(a)-1):length(a))=[1 -L11(length(a))];
 
for k=1:nu1-2
    a(k+1,k:k+2)=[1 -L11(k+1) M11(k+1)];
end

b(1,1:2)=[-L12(2) M12(2)];
b(length(b),(length(b)-1):length(b))=[1 -W12];
 
for k=1:nv1-2
    b(k+1,k:k+2)=[1 -L12(k+2) M12(k+2)];
end
c =zeros(nu1,1);
d=zeros(nv1,1);
for k=2:nu1-1
c(k)=-R11(k)*p(1,k);
end
for k=2:nv1-1
    d(k)=-R12(k+1)*q(1,k+1);
end

c(1)=Ohm11*p(1,1);  c(nu1)=-R11(nu1)*p(1,nu1)-M11(nu1)*Ca;
d(1)=-R12(2)*q(1,2)-Cb;  d(nv1)=-Ohm12*q(1,nv1+1);

% solve for concentrations
%%p(j+1,1:nu1)=a\c;
p(1,1:nu1)=a\c;
q(1,2:nv1+1)=b\d;



elseif s(j+1)<s(j)
%case b: s(j+1)<s(j)
for i=2:nu1
L11(i)=(s(j+1)*(u(i)-u(i-1))+((dt*Da)/(s(j+1)*(u(i+1)-u(i))))+((dt*Da)/(s(j+1)*(u(i)-u(i-1))))-(s(j+1)-s(j))*u(i))/...
    ((s(j+1)-s(j))*u(i-1)-(dt*Da)/(s(j+1)*(u(i)-u(i-1))));
M11(i)=((dt*Da)/(s(j+1)*(u(i+1)-u(i))))/((s(j+1)-s(j))*u(i-1)-((dt*Da)/(s(j+1)*(u(i)-u(i-1)))));

R11(i)=((s(j)*(u(i)-u(i-1))))/((s(j+1)-s(j))*u(i-1)-((dt*Da)/(s(j+1)*(u(i)-u(i-1)))));

W11=((dt*Da)/(s(j+1)*(u(2))))/(((dt*Da)/(s(j+1)*(u(2)))+(s(j+1)*(u(1)))-(s(j+1)-s(j))*u(1)));

Ohm11=(s(j)*u(1))/(((dt*Da)/(s(j+1)*(u(2)))+(s(j+1)*(u(1)))-(s(j+1)-s(j))*u(1)));
end

%for case a: u(0.5)=u(1) for case b: u(0.5)=u(0)
for i=2:nv1
L12(i)=((l-s(j+1))*(v(i)-v(i-1))+((dt*Db)/((l-s(j+1))*(v(i+1)-v(i))))+((dt*Db)/((l-s(j+1))*(v(i)-v(i-1))))-...
    (s(j+1)-s(j))*(1-v(i-1)))/((s(j+1)-s(j))*(1-v(i-1))-((dt*Db)/((l-s(j+1))*(v(i+1)-v(i)))));

M12(i)=((dt*Db)/((l-s(j+1))*(v(i+1)-v(i))))/((s(j+1)-s(j))*(1-v(i-1))-((dt*Db)/((l-s(j+1))*(v(i+1)-v(i)))));

R12(i)=((l-s(j))*(v(i)-v(i-1))/((s(j+1)-s(j))*(1-v(i-1))-((dt*Db)/((l-s(j+1))*(v(i+1)-v(i))))));

W12=(((l-s(j+1))*(-v(nv-1))+((dt*Db)/((l-s(j+1))*(v(nv)-v(nv-1)))))/...
    ((s(j+1)-s(j))*(1-v(nv-1))-((dt*Db)/((l-s(j+1))*(v(nv)-v(nv-1))))));

Ohm12=((l-s(j))*(-v(nv-1))/((s(j+1)-s(j))*(1-v(nv-1))-((dt*Db)/((l-s(j+1))*(v(nv)-v(nv-1))))));
end


a(1,1:2)=[1 -W11];
a(length(a),(length(a)-1):length(a))=[1 L11(length(a))];
for k=1:nu1-2
    a(k+1,k:k+2)=[1 L11(k+1) -M11(k+1)];
end
b(1,1:2)=[L12(2) -M12(2)];
b(length(b),(length(b)-1):length(b))=[1 W12];
for k=1:nv1-2
    b(k+1,k:k+2)=[1 L12(k+2) -M12(k+2)];
end
c = zeros (nu1,1);
d=zeros(nv1,1);


for k=2:nu1-1
c(k)=R11(k)*p(1,k);
end
for k=2:nv1-1
    d(k)=R12(k+1)*q(1,k+1);
end

c(1)=Ohm11*p(1,1);  c(nu1)=R11(nu1)*p(1,nu1)+M11(nu1)*Ca;
d(1)=R12(2)*q(1,2)-Cb;  d(nv1)=Ohm12*q(1,nv1+1);


% solve for concentrations
p(1,1:nu1)=a\c;
q(1,2:nv1+1)=b\d;

    


end
s1(j+1)=((((Db*dt)/(l-sjsig))*((qjsig(2)-Cb)/(v(2)))-((Da*dt)/(sjsig))*((Ca-pjsig(nu-1))/(1-u(nu-1))))/...
    ((0.5*(1+u(nu-1))*piplusend)+(0.5*(1-u(nu-1))*Ca)-((1-0.5*v(2))*qiplusbeg)-((0.5*v(2))*Cb)))+s(j);


if abs((s(j+1)-s1(j+1)))<0.001; %0.0000000001;     
    break;
end

clear s1 %
clear L11 M11 R11 %,W11,Ohm11
clear L12 M12 R12 %,W11,Ohm11
clear a b c d


end

%if t(j)>1; %s(j)<=0.01;
 %   break;
%end

if s(j+1)>s(j)

p(1,nu)=Ca;
q(1,1)=Cb;
piplusend=p(1,nu);
pjsig=p(1,1:nu);
qjsig=q(1,1:nv);
qiplusbeg=q(1,2);
    
    
elseif s(j+1)<s(j)
    
 p(1,nu)=Ca;
 q(1,1)=Cb;
piplusend=p(1,nu-1);
pjsig=p(1,1:nu);
qjsig=q(1,1:nv);
qiplusbeg=q(1,1);
    
    
end

if s(j+1)<=0.1; %abs((s(j+1)-s1(j+1)))<0.01; %0.0000000001;     
   % xlswrite('C:\Users\Muideen AJIBOYE\Documents\storedtest1.xlsx',[(t(1:j+1))' (s(1:j+1))'],1,xlscol((3*(rgnu-j)/(nt-1))+1))
 semilogx(t,s,'LineWidth',2)    
       hold on
 
   break;
elseif t(j+1)>=10 %00000;
    %xlswrite('C:\Users\Muideen AJIBOYE\Documents\storedtest1.xlsx',[(t(1:j+1))' (s(1:j+1))'],1,xlscol((3*(rgnu-j)/(nt-1))+1))
       semilogx(t,s,'LineWidth',2)    
       
       hold on
    break;
elseif j==nt-1;
    %xlswrite('C:\Users\Muideen AJIBOYE\Documents\storedtest1.xlsx',[t(1:nt-1)' (s(1:nt-1))'],1,xlscol((3*rgnu/(nt-1))-2))
    if t(j+1)<0.51;
    semilogx(t(1000:end),s(1000:end),'LineWidth',2)    
       hold on
    s(1)=s(j+1);
    t(1)=t(j+1);
    j=0;
    else
 semilogx(t,s,'LineWidth',2)    
       hold on
    s(1)=s(j+1);
    t(1)=t(j+1);
    j=0;
    end
    
end

j=j+1;
rgnu=1+rgnu;


end




toc(t1)

%end
grid
xlabel('t(secs)')
ylabel('s')
title('s(t)')
