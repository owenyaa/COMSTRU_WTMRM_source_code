clear;clc;close all;
global Tdata parameter_a N_dof N_harm N_w0 index_global N_min alpha
%% different initial values
% ��ͬϵͳϵ��Ҫ���ĵĲ���
N_dof=2;N_harm=72;N_w0=2;%��Ƶ����
N_min=35; alpha=0.02;
% w1=0.85137;
% w=1;
% % w0=zeros(1,2*N_dof);
% w0(1,1:N_w0)=[w1,w];
Tdata=0:0.05:400;
%% �����Ƶ�����ϵ��
index_global=[1 0;0 1];
for N=3:2:11
    temp=[N 0];
    for n=N-1:-1:1
        m=N-n;
        temp=[temp;n m;n -m];
    end
    temp=[temp;0 N];
    index_global=[index_global;temp];
end
index_global=[index_global,index_global];

% load 'w_1_N_11_Nmax_72_5_10_4.mat';
load 'w_1_N_11_Nmax_72_5_10_4.mat';
% load '2_10_7_128.mat';
% Tdata=0:0.1:40000;
% residual=cal_residual(parameter_a);
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% ��X,Y�������ɶ�
fundamental_w=parameter_a(1,1:N_w0);
vector_w=index_global(:,1:N_w0)*fundamental_w';
Harm_parameter_a=parameter_a(2:end,:);
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
        dx(j,:)=dx(j,:)-vector_w(i)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i)*Tdata)+vector_w(i)*Harm_parameter_a(i,2*j)*cos(vector_w(i)*Tdata);
        ddx(j,:)=ddx(j,:)-(vector_w(i))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)-(vector_w(i))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
    end
end
residua=cal_residual(parameter_a);

figure;
plot(Tdata,residua(:,1),'r-','LineWidth',1);
hold on;
plot(Tdata,residua(:,2),'k-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



figure;
plot(Tdata,x(1,:),'k.','LineWidth',1);
hold on;
plot(Tdata,x(2,:),'b.','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

figure;
plot(x(1,:),dx(1,:),'k-','LineWidth',1);
hold on;
plot(x(2,:),dx(2,:),'b-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
