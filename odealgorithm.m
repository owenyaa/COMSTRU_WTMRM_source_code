%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2016-11-03 17:10
% Last Revised : GUANG_LIU ,2016-11-03
% Remark : 本程序的结构、功能和变量做一下说明%

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
% figure。
% plot(f,Pyy(1:(N_fft/2+1))/10000,'r-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% save('F:\matlab.mat','dt','tt','num','-v7.3');


% % fs=1/dt;%采样频率
% % % 采样频率与时间间隔之间的关系： fs=1/dt
% % % 采样定理告诉我们，采样频率要大于信号频率的两倍。 
% % N=2^23;  %采样点数2^17
% % % N个采样点，经过FFT之后，就可以得到N个点的FFT结果。为了方便进行FFT运算，通常N取2的整数次方。
% % % 要精确到xHz，则需要采样长度为1/x秒的信号，并做FFT。
% % % 要提高频率分辨率，就需要增加采样点数
% % n=0:N-1;
% % t=n/fs;  % dt=1/fs 表示时间间隔   fs=1/dt
% % y=fft(num(1:end,1),N);  % 进行fft变换
% % % 假设采样频率为Fs，信号频率F，采样点数为N。那么FFT之后结果就是一个为N点的复数。
% % % 每一个点就对应着一个频率点。这个点的模值，就是该频率值下的幅度特性。
% % % y % 输出y看看fft之后的结果。
% % m=abs(y(1:N/2))*2/N; % 求信号的真实幅值
% % f=2*pi*n*fs/N;  % 注意 乘2*pi 和不乘2*pi的区别%m=log10(m);
% % figure;
% % % plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
% % plot(f(1:N/2),log(m(1:N/2)),'r-','LineWidth',1);
% % xlim([0,15]);
% % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% % h1=legend('$$Iteration steps$$');
% % set(h1,'Interpreter','latex','FontSize',15);
% % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% % 
% % y=fft(num(1:end,2),N);  % 进行fft变换
% % % 假设采样频率为Fs，信号频率F，采样点数为N。那么FFT之后结果就是一个为N点的复数。
% % % 每一个点就对应着一个频率点。这个点的模值，就是该频率值下的幅度特性。
% % % y % 输出y看看fft之后的结果。
% % m=abs(y(1:N/2))*2/N; % 求信号的真实幅值
% % f=2*pi*n*fs/N;  % 注意 乘2*pi 和不乘2*pi的区别%m=log10(m);
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
