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
B=[0;0;0;0];
C=[0 0 1 0];
D=[0];

% Gp=ss(A,B,C,D);
% step(Gp);
[Np,Dp]= ss2tf(A,B,C,D);
G0 = tf(Np,Dp);
pole(G0)
% T = [0:0.001:1];
% U = zeros(1,length(T));
% x0 = [0;0;0;0];
% lsim(Gp,U,T,x0)
target_poles = [-11;-2;-0.5;-12];
K = place(A,B,target_poles);
% eig(A-(B*K));
Nbar = -1/(C*inv(A-(B*K))*B)
Gcl = ss(A-(B*K),Nbar*B,C,D);
step(Gcl)
stepinfo(Gcl)