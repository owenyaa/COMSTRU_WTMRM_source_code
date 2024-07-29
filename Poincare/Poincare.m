clear;clc;close all;
global A 
global k13 k14 k15 k16 k23 k24 k25 k26
global c11 c12 c13 c14 c15 c16 c21 c22 c23 c24 c25 c26

A1=0.2;A2=0.02;
b=40*10^(-9);h=20*10^(-9);E=76*10^9;E_s=1.22;L=30*h;
tau_s=0.570761;
H_s=2*tau_s*b;E_I=E*b*h^3/12+E_s*h^3/6+E_s*b*h^2/2;
beta=H_s*L^2/(E_I);mu=2*E_s*(h/b+3)/(E*h);

k11=12.3596-0.8588*beta+12.3596*mu;k12=0.0068+11.7444*beta+0.0068*mu;k13=-16.3892-2.1571*beta+40.4251*mu;
k14=1528.0120+41.2266*beta-306.8542*mu;k15=-12783.2313-147.2032*beta+2483.9228*mu;k16=8597.5497+192.7130*beta-2176.2470*mu;

k21=0.0002-1.8739*beta+0.0002*mu;k22=485.4811+13.2938*beta+485.4811*mu;k23=-57.8499-3.4640*beta-102.2985*mu;
k24=-872.6730-32.0639*beta+2484.1643*mu;k25=-15292.1786-44.4417*beta-6530.1215*mu;k26=-56843.5188+53.2269*beta+13415.5523*mu;

c11=4.5968;c12=-3.5964;c13=-7.1928;c14=12.2359;c15=25.1741;c16=-22.1929;
c21=-3.5963;c22=25.1733;c23=12.2356;c24=-44.3844;c25=22.1922;c26=144.7195;


% M=[1,0;0,ep];
% C=[ep*lamda,-ep*lamda;-ep*lamda,ep*lamda];K=[1,0;0,0];
A=[zeros(2),eye(2);-[k11,k12;k21,k22],zeros(2)];
% odex=[x(1,1);x(2,1);dx(1,1);dx(2,1)];
odex=[A1;A2;0;0];
dt=0.0002;
Tdata=0:dt:1000;
% Tdata=0:0.004:50;
options=odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,num]=ode45('odehomo',Tdata,odex,options);


plot3(num(1:end,1),num(1:end,3),num(1:end,2)) %画出三维图
%poincare截面图，即用Z=Z0这个截面，去截上面的三维图，截面留下的痕迹点就是poincare截面图
num1=num(:,1);
num2=num(:,3);
num3=num(:,2);
% 设定x2的截面，理论上可以任意选取
section_num=0.01;
section=num3-section_num;
j=1;
for i=1:length(section)-1
    %前后两个点乘积小于0，表示此时穿过截面
    if section(i,1)*section(i+1,1)<0
        %看前后哪个点离截面更近，保存该点的位置
        %真正的poincare图，应该是用两点的插值，此处用的其实是近似点,要求时间步长很小才会很准
        if abs(section(i,1))<abs(section(i+1,1))
            jj(j)=i;
            j=j+1;
        else
            jj(j)=i+1;
            j=j+1;
        end
    end
end

for i=1:length(jj)
    x1(i)=num1(jj(i));
    x2(i)=num2(jj(i));
end
figure;
plot(x1(1:3:end),x2(1:3:end),'k.')
%画出poincare截面图
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



