clc;
clear;
m=2105;
m_fuel=88;
a=0.53;
b=-1.43;
d=1;
J_22=1883;
F0=2450;

alpha=m*b-m*a;
beta=m*a;
gama=((1/m)-(d/b)*((1/m)+(1/m_fuel)));
omega=(b-(J_22/b)*((1/m)+(1/m_fuel))-a);

%phi2=gama*f*alpha-omega*f-F0*phi*omega/(omega*beta-alpha*a);
A=[0 1;-F0/(omega*beta-alpha*a) 0];
B=[0;gama*alpha-omega/(omega*beta-alpha*a)];
C=[1 0];
D=[0];
Gp=ss(A,B,C,D);
step(Gp);
[Numcl,Dencl] = ss2tf(A,B,C,D);
G0 = tf(Numcl,Dencl);
pole(G0)
zero(G0)
target_poles=[-0.5;-11];
K=place(A,B,target_poles);
%eigvalue=eig(A-(B*K))
Kapp=[-0.1484 -0.3184];
nbar=-1/(C*inv(A-(B*Kapp))*B);
Gp = ss(A-(B*Kapp),nbar*B,C,D);
step(Gp);
stepinfo(Gp);% overshoot=0; settletime=7.9
% 
% y=[1 0]*x;