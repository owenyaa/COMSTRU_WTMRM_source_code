clc;clear;close all;
global tf h N_dof N_harm N_w0 Tdata index_global index %w
%% parameters
N_dof=2;%N_harm=8;
N_w0=2;%基频个数
% N_min=30; alpha=0.02;
N_dof1=6;




w1=3.510628883212891;wd=1.498777127226156;
% w0=zeros(1,2*N_dof);
w0(1,1:N_w0)=[w1,wd];
% 其他参数
tf=2*pi/min(w0(1,1:2));h=2*pi/(min(w0(1,1:2))*1000);Tdata=(0:h:tf);
%% 计算基频的组合系数
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
% parameter_a(2,:)=[0.4,0.1,0.1,0.1,0.2,0.1,];% 一阶谐波初值

parameter_a=zeros(N_harm+1,2*N_dof);
load 'Copy_of_NS_ini_parameter_a.mat';
% % parameter_a(2:end,:)=-sqrt(parameter_a(2:end,:).^2/2);
parameter_a(1,:)=[w1,wd,0,0];
parameter_a(2:end,:)=parameter_a1(2:end,:);

%% 计算残差
residual=cal_residual(parameter_a);
%% 绘制残差曲线
% figure;
% plot(Tdata,residual(:,1),'r-','LineWidth',1);
% hold on;
% plot(Tdata,residual(:,2),'k-','LineWidth',1);
% hold on;
% plot(Tdata,residual(:,3),'b-','LineWidth',1);
% h1=legend('$$h$$','$$\alpha$$','$$\beta$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

%% 此处为灵敏度验证程序，本质为用差分代替灵敏度，可以验证灵敏度是否算错
% 计算正问题之后，将某个参数减去一个小量，用新参数再算一次
% 两次所得结果做差再除以小量即为灵敏度，验证差分的灵敏度和直接计算的灵敏度曲线是否重合
% 但是要注意，差分结果肯定是对的(即下文差量结果)，可能出错的是原程序的内容，x_cal中内容
% 另外关于对应，差分结果(parameter_a1)中的位移(x1(1,:))，速度，加速度对应到原程序(parameter_a)中的灵敏度内容x_cal(2,:)
%
% for i=1:2*N_harm*N_dof
for i=1:4
    ddt=0.000001;
    %初始化参数
    w0=zeros(1,2*N_dof);
    w0(1,1:N_w0)=[w1,wd];
    parameter_a(1,:)=w0;
    parameter_a(2:end,:)=parameter_a1(2:end,:);
    
    % 参数系数矩阵，1代表求此处参数的差分
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    temp_real_w0=zeros(1,2*N_dof);
    sensitivity_parameter_a=[temp_real_w0;sensitivity_parameter_a];
    parameter_a=parameter_a+ddt*sensitivity_parameter_a;% 增量加上去了
    
    residual1=cal_residual(parameter_a);
    aaaa=(residual1-residual)/ddt;
    figure;
    plot(Tdata,residual(:,N_dof1*(i+N_w0)+1),'r-')
    hold on
    plot(Tdata,residual(:,N_dof1*(i+N_w0)+2),'r-')
    hold on
    plot(Tdata,residual(:,N_dof1*(i+N_w0)+3),'r-')
    hold on
    plot(Tdata,residual(:,N_dof1*(i+N_w0)+4),'r-')
    hold on
    plot(Tdata,residual(:,N_dof1*(i+N_w0)+5),'r-')
    hold on
    plot(Tdata,residual(:,N_dof1*(i+N_w0)+6),'r-')
    hold on
    
    plot(Tdata,aaaa(:,1),'k-')
    hold on
    plot(Tdata,aaaa(:,2),'k-')
    hold on
    plot(Tdata,aaaa(:,3),'k-')
    hold on
    plot(Tdata,aaaa(:,4),'k-')
    hold on
    plot(Tdata,aaaa(:,5),'k-')
    hold on
    plot(Tdata,aaaa(:,6),'k-')
end


% for i=1:2
%     % i=8;
%     ddt=0.0000000001;
%     w0=zeros(1,2*N_dof);
%     w0(1,1:N_w0)=[w1,wd];
%     w0(1,i)=w0(1,i)+ddt;
%     parameter_a(1,:)=w0;
%     parameter_a(2:end,:)=parameter_a1(2:end,:);
% 
%     residual1=cal_residual(parameter_a);
%     aaaa=(residual1-residual)/ddt;
%     figure;
%     plot(Tdata,residual(:,i*N_dof1+1),'r-')
%     hold on
%     plot(Tdata,residual(:,i*N_dof1+2),'r-')
%     hold on
%     plot(Tdata,residual(:,i*N_dof1+3),'r-')
%     hold on
%     plot(Tdata,residual(:,i*N_dof1+4),'r-')
%     hold on
%     plot(Tdata,residual(:,i*N_dof1+5),'r-')
%     hold on
%     plot(Tdata,residual(:,i*N_dof1+6),'r-')
%     hold on
% 
%     plot(Tdata,aaaa(:,1),'k-')
%     hold on
%     plot(Tdata,aaaa(:,2),'k-')
%     hold on
%     plot(Tdata,aaaa(:,3),'k-')
%     hold on
%     plot(Tdata,aaaa(:,4),'k-')
%     hold on
%     plot(Tdata,aaaa(:,5),'k-')
%     hold on
%     plot(Tdata,aaaa(:,6),'k-')
% 
% end

















