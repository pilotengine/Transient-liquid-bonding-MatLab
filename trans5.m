t1=tic;


%function[p,q,s,t] = transabs(~,~,~,~,~,~)

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
Tend=0.0001;
dt=0.0000001;

t=0:dt:Tend;
t(1)=0;
nt=length(t);

nu=500; %number of coordinate point (length of base metal 6025/elemental thickness 25)+1
%nu=51; % look as been work
%nu=76;
du=1/((nu-1));
l=(l1/2);%+l2;
p0=zeros(1,nu);
%pjsig0=zeros(1,nu);
%piplus=zeros(1,nu);
p=zeros(nt,nu);
p(2:nt,nu)=Ca;

p(1,:)=p01;

%for i=1:nu
p0(1,1:nu)=p01;
piplusend=p01;
pjsig(1,1:nu)=p0;
%pjsig=pjsig0;
%piplus=pjsig0;

%end
nv=500; % look as been work
dv=1/((nv-1));
q0=zeros(1,nv);
qjsig0=zeros(1,nv);
%qiplus=zeros(1,nv);
q=zeros(nt,nv);
q(2:nt,1)=Cb;

q(1,:)=q01;

%for i=1:nv
q0(1,1:nv)=q01;
qiplusbeg=q01;

qjsig0(1,1:nv)=q0;
qjsig=qjsig0;
u=zeros(1,nu);
v=zeros(1,nv);

for i=1:nu
u(i)=(i-1)*du;
end
for i=1:nv
v(i)=(i-1)*dv;
end

  
%time taken t=Tend


 %  moving boundary position initiation

s=zeros(1,nt);
s(1)=25*0.5; %0.000025;



sjsig=s(1);
s1=zeros(1,nt);

j=1;
while 1
  
t(j+1)=(j)*dt;
while 1
    
%while abs(s(j+1)-s(j))>0.00000001
%if j==1
s(j+1)=((((Db*dt)/(l-sjsig))*((qjsig(2)-Cb)/(v(2)))-((Da*dt)/(sjsig))*((Ca-pjsig(nu-1))/(1-u(nu-1))))/...
    ((0.5*(1+u(nu-1))*piplusend)+(0.5*(1-u(nu-1))*Ca)-((1-0.5*v(2))*qiplusbeg)-((0.5*v(2))*Cb)))+s(j);

sjsig=s(j+1);



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
c(k)=-R11(k)*p(j,k);
end
for k=2:nv1-1
    d(k)=-R12(k+1)*q(j,k+1);
end

c(1)=Ohm11*p(j,1);  c(nu1)=-R11(nu1)*p(j,nu1)-M11(nu1)*Ca;
d(1)=-R12(2)*q(j,2)-Cb;  d(nv1)=-Ohm12*q(j,nv1+1);

% solve for concentrations
p(j+1,1:nu1)=a\c;
q(j+1,2:nv1+1)=b\d;
%p(j+1,1:nu1)=pinv(a)*c;
%q(j+1,2:nv1+1)=pinv(b)*d;
%p(j+1,1:nu1)=bicgstab(a,c);
%q(j+1,2:nv1+1)=bicgstab(b,d);


piplusend=p(j+1,nu);
pjsig=p(j+1,1:nu);
qjsig=q(j+1,1:nv);
qiplusbeg=q(j+1,2);

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
c(k)=R11(k)*p(j,k);
end
for k=2:nv1-1
    d(k)=R12(k+1)*q(j,k+1);
end

c(1)=Ohm11*p(j,1);  c(nu1)=R11(nu1)*p(j,nu1)+M11(nu1)*Ca;
d(1)=R12(2)*q(j,2)-Cb;  d(nv1)=Ohm12*q(j,nv1+1);


% solve for concentrations
p(j+1,1:nu1)=a\c;
q(j+1,2:nv1+1)=b\d;
%p(j+1,1:nu1)=pinv(a)*c;
%q(j+1,2:nv1+1)=pinv(b)*d;

%p(j+1,1:nu1)=bicgstab(a,c);
%q(j+1,2:nv1+1)=bicgstab(b,d);

piplusend=p(j+1,nu-1);
pjsig=p(j+1,1:nu);
qjsig=q(j+1,1:nv);
qiplusbeg=q(j+1,1);


end


s1(j+1)=((((Db*dt)/(l-sjsig))*((qjsig(2)-Cb)/(v(2)))-((Da*dt)/(sjsig))*((Ca-pjsig(nu-1))/(1-u(nu-1))))/...
    ((0.5*(1+u(nu-1))*piplusend)+(0.5*(1-u(nu-1))*Ca)-((1-0.5*v(2))*qiplusbeg)-((0.5*v(2))*Cb)))+s(j);





%{if j>=2;
 %   break;
%elseif abs((s(j+1)-s1(j+1)))<0.0000000001;     
 %   break;
%end %}

if abs((s(j+1)-s1(j+1)))<0.000001; %0.0000000001;     
    break;
end





end

if t(j)>0.1;
    break;
end
  j=j+1;

end


%xlswrite('C:\Users\Muideen AJIBOYE\Desktop\restoredefault\Manitoba work\report\malikexcel.xlsx',t',1,'A4')
%xlswrite('C:\Users\Muideen AJIBOYE\Desktop\restoredefault\Manitoba work\report\malikexcel.xlsx',s',1,'B4')
figure
subplot(2,1,2)

plot(t,s)
grid
xlabel('t')
ylabel('s')
title('s(t)')

 
%axis



%end

 toc(t1)
