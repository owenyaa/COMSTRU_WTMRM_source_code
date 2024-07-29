clear;clc;close all;
global Tdata parameter_a N_dof N_harm N_w0 index_global R_dof vector_w index
%% different initial values
% 不同系统系需要更改的参数
N_dof=2;N_w0=2;R_dof=6;%基频个数
%N_min=30; alpha=0.02;

w1=3.510628883212891;wd=1.498777127226156;
% w1=3.51067997274465;wd=1.49483282433107;
% w0=zeros(1,2*N_dof);
w0(1,1:N_w0)=[w1,wd];
Tdata=0:0.004:30;
index=25;
%% 计算基频的组合系数
%% 计算基频的组合系index_global=[1,-1,1,2;3,-1,1,2...]
for i=1:2:index
    index_global((i+1)/2,1)=i;
    index_global((i+1)/2,2:4)=[-1,1,2];
end
for i=1:2:index_global(end,1)
    temp_vector_w((i+1)/2,1)=index_global((i+1)/2,1)*w0(1,1);
    temp_vector_w((i+1)/2,2:4)=index_global((i+1)/2,1)*w0(1,1)+w0(1,2)*index_global((i+1)/2,2:4);
end
size_temp_vector_w=size(temp_vector_w);
N_harm=size_temp_vector_w(1,1)*size_temp_vector_w(1,2);

%% 第一行存储频率，后面存储谐波系数,每个自由度两列
% parameter_a=[w0,  0,    0,   0;
%             C_11,S_11,C_21,S_21;
%             C_12,S_12,C_22,S_22;
%             C_13,S_13,C_23,S_23;...];
parameter_a=zeros(N_harm+1,2*N_dof);
load 'NS_ini_parameter_a.mat';
% parameter_a(1,1)=3.5;
residua=cal_residual(parameter_a);

figure;
plot(Tdata,residua(:,1),'r-','LineWidth',1);
hold on;
plot(Tdata,residua(:,2),'k-','LineWidth',1);
hold on;
plot(Tdata,residua(:,3),'k-','LineWidth',1);
hold on;
plot(Tdata,residua(:,4),'k-','LineWidth',1);
hold on;
plot(Tdata,residua(:,5),'k-','LineWidth',1);
hold on;
plot(Tdata,residua(:,6),'k-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


% Tdata=0:0.01:20;
% 有X,Y两个自由度
%% 计算两个基频率的组合，注意，频率组合要去掉负数频率
fundamental_w=parameter_a(1,1:N_w0);

for i=1:2:index_global(end,1)
    temp_vector_w((i+1)/2,1)=index_global((i+1)/2,1)*fundamental_w(1,1);
    temp_vector_w((i+1)/2,2:4)=index_global((i+1)/2,1)*fundamental_w(1,1)+fundamental_w(1,2)*index_global((i+1)/2,2:4);
end
size_temp_vector_w=size(temp_vector_w);
vector_w=[];
for i=1:size_temp_vector_w(1,1)
    vector_w=[vector_w;temp_vector_w(i,:)'];
end

Harm_parameter_a=parameter_a(2:end,:);
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有三个自由度
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
        dx(j,:)=dx(j,:)-vector_w(i)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i)*Tdata)+vector_w(i)*Harm_parameter_a(i,2*j)*cos(vector_w(i)*Tdata);
        ddx(j,:)=ddx(j,:)-(vector_w(i))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)-(vector_w(i))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
    end
end
figure;
plot(Tdata,x(1,:),'r-','LineWidth',1);
hold on;
plot(Tdata,x(2,:),'k-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
plot(x(1,:),dx(1,:),'r-','LineWidth',1);
hold on;
plot(x(2,:),dx(2,:),'k-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);




load matlab.mat;
figure;
plot(Tdata,x(1,:)'-num(:,1),'r-','LineWidth',1);
hold on;
plot(Tdata,x(2,:)'-num(:,2),'k-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);











