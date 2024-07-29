%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2016-11-03 17:10
% Last Revised : GUANG_LIU ,2016-11-03
% Remark : ������Ľṹ�����ܺͱ�����һ��˵��%

clear;clc;close all;
%% 
global A 
global k13 k14 k15 k16 k23 k24 k25 k26
global c11 c12 c13 c14 c15 c16 c21 c22 c23 c24 c25 c26

A1=0.2;A2=0.02;
b=40*10^(-9);h=20*10^(-9);E=76*10^9;E_s=1.22;L=50*h;
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
dt=0.001;
% Tdata=0:dt:30;
Tdata=0:0.004:50;
options=odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,num]=ode45('odehomo',Tdata,odex,options);


figure;
subplot(2,1,1);
plot(Tdata,num(:,1),'r-','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
subplot(2,1,2);
plot(Tdata,num(:,2),'b-','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



figure;
plot(Tdata,num(:,1),'k-');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
plot(Tdata,num(:,2),'r-');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% 
figure;
plot(num(1:end,1),num(1:end,3),'k-');
hold on
plot(num(1:end,2),num(1:end,4),'r-');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


% figure;
% plot(tt(1:50:end),num(1:50:end,1),'k.');
% hold on
% plot(tt(1:50:end),num(1:50:end,2),'r.');
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
% plot(num(70000:end,1),num(70000:end,3),'k-');
% hold on
% plot(num(70000:end,2),num(70000:end,4),'r-');
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


% residual=-M\(C*num(:,3:4)'+K*num(:,1:2)'+[-ep*f1*cos(w*tt)+ep*k_n*(num(:,1)'-num(:,2)').^3;ep*k_n*(num(:,2)'-num(:,1)').^3]);
% figure;
% plot(tt,residual(1,:)-ddx(1,:),'k-');
% hold on
% plot(tt,residual(2,:)-ddx(2,:),'r-');
% 
% figure;
% plot(tt,ddx(1,:),'r-');
% hold on
% plot(tt,ddx(2,:),'k-');

% dt=tt(2)-tt(1);
% data=num(1:end,1);
% N=length(data);
% N_fft=2^16;
% Y=fft(data,N_fft);
% Pyy=Y.*conj(Y)/N_fft;
% f=1/dt*(0:N_fft/2)/N_fft;
% figure��
% plot(f,Pyy(1:(N_fft/2+1))/10000,'r-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% save('F:\matlab.mat','dt','tt','num','-v7.3');


% % fs=1/dt;%����Ƶ��
% % % ����Ƶ����ʱ����֮��Ĺ�ϵ�� fs=1/dt
% % % ��������������ǣ�����Ƶ��Ҫ�����ź�Ƶ�ʵ������� 
% % N=2^23;  %��������2^17
% % % N�������㣬����FFT֮�󣬾Ϳ��Եõ�N�����FFT�����Ϊ�˷������FFT���㣬ͨ��Nȡ2�������η���
% % % Ҫ��ȷ��xHz������Ҫ��������Ϊ1/x����źţ�����FFT��
% % % Ҫ���Ƶ�ʷֱ��ʣ�����Ҫ���Ӳ�������
% % n=0:N-1;
% % t=n/fs;  % dt=1/fs ��ʾʱ����   fs=1/dt
% % y=fft(num(1:end,1),N);  % ����fft�任
% % % �������Ƶ��ΪFs���ź�Ƶ��F����������ΪN����ôFFT֮��������һ��ΪN��ĸ�����
% % % ÿһ����Ͷ�Ӧ��һ��Ƶ�ʵ㡣������ģֵ�����Ǹ�Ƶ��ֵ�µķ������ԡ�
% % % y % ���y����fft֮��Ľ����
% % m=abs(y(1:N/2))*2/N; % ���źŵ���ʵ��ֵ
% % f=2*pi*n*fs/N;  % ע�� ��2*pi �Ͳ���2*pi������%m=log10(m);
% % figure;
% % % plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
% % plot(f(1:N/2),log(m(1:N/2)),'r-','LineWidth',1);
% % xlim([0,15]);
% % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% % h1=legend('$$Iteration steps$$');
% % set(h1,'Interpreter','latex','FontSize',15);
% % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% % 
% % y=fft(num(1:end,2),N);  % ����fft�任
% % % �������Ƶ��ΪFs���ź�Ƶ��F����������ΪN����ôFFT֮��������һ��ΪN��ĸ�����
% % % ÿһ����Ͷ�Ӧ��һ��Ƶ�ʵ㡣������ģֵ�����Ǹ�Ƶ��ֵ�µķ������ԡ�
% % % y % ���y����fft֮��Ľ����
% % m=abs(y(1:N/2))*2/N; % ���źŵ���ʵ��ֵ
% % f=2*pi*n*fs/N;  % ע�� ��2*pi �Ͳ���2*pi������%m=log10(m);
% % hold  on;
% % % plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
% % plot(f(1:N/2),m(1:N/2),'k-','LineWidth',1);
% % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% % h1=legend('$$Iteration steps$$');
% % set(h1,'Interpreter','latex','FontSize',15);
% % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% hold on;
% plot(tt(1:30:end),num(1:30:end,1),'k.','MarkerSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
% plot(num(1:30:end,1),num(1:30:end,2),'k.','MarkerSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


% aa=ddx+e*(1-x.^2).*dx+w0_2*x+a*x.^3-f1*cos(w1*Tdata)-f2*cos(w2*Tdata);
% figure;
% plot(num(1:end,1),num(1:end,2),'k-','LineWidth',1.5);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
