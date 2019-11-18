clc;
clear;
m=2105;
m_fuel=88;
a=0.53;
b=-1.43;
d=1;
J_22=1883;
F0=2450;

arefa=m*b-m*a;
beta=m*a;
gama=((1/m)-(d/b)*((1/m)+(1/m_fuel)));
omega=(b-(J_22/b)*((1/m)+(1/m_fuel))-a);

num_1=((a*F0)/(beta*omega-a*arefa));
num_2=-((F0*omega)/(omega*beta-arefa*arefa));
num_3=(a-beta*gama)/(beta*omega-a*arefa);
num_4=(gama*arefa-omega)/(omega*beta-arefa*a);

A=[0 1 0 0;0 0 num_1 0;0 0 0 1;0 0 num_2 0];
B=[0;num_3;0;num_4];
C=[1 0 0 0];
D=[0];

[Numcl,Dencl] = ss2tf(A,B,C,D);
G0 = tf(Numcl,Dencl);
pole(G0);
zero(G0);
target_poles=[-0.5;-1.3;-11;-12];
K=place(A,B,target_poles)
%eigvalue=eig(A-(B*K))
Kapp=[-410000 -1188400 -14100 212000];
nbar=-1/(C*inv(A-(B*Kapp))*B);
Gp = ss(A-(B*Kapp),nbar*B,C,D);
step(Gp);
stepinfo(Gp)% overshoot=0; settletime=9
L=place(A',C',[-20 -23 -24 -25]');
Lapp=[100;3200;646900;3676900];
% eigvaluec=eig(A-(Lapp*C));
Acl = [A-(B*Kapp) B*Kapp; zeros(size(A)) A-(Lapp*C)];
Bcl = nbar*[B;zeros(length(A),1)];
Ccl = [C zeros(1,length(A))];
Dcl = D;
[Numcl1,Dencl1] = ss2tf(Acl,Bcl,Ccl,Dcl);

Gcl = tf(Numcl1,Dencl1);
step(Gcl)
stepinfo(Gcl)