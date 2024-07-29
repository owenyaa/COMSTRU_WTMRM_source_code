clear;clc;close all;
global Tdata parameter_a N_dof N_harm N_w0 index_global R_dof vector_w index
%% different initial values
% ��ͬϵͳϵ��Ҫ���ĵĲ���
N_dof=2;N_w0=2;R_dof=6;%��Ƶ����
%N_min=30; alpha=0.02;

w1=3.510628883212891;wd=1.498777127226156;
% w1=3.51067997274465;wd=1.49483282433107;
% w0=zeros(1,2*N_dof);
w0(1,1:N_w0)=[w1,wd];
Tdata=0:0.004:30;
index=25;
%% �����Ƶ�����ϵ��
%% �����Ƶ�����ϵindex_global=[1,-1,1,2;3,-1,1,2...]
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

%% ��һ�д洢Ƶ�ʣ�����洢г��ϵ��,ÿ�����ɶ�����
% parameter_a=[w0,  0,    0,   0;
%             C_11,S_11,C_21,S_21;
%             C_12,S_12,C_22,S_22;
%             C_13,S_13,C_23,S_23;...];
parameter_a=zeros(N_harm+1,2*N_dof);
load 'NS_ini_parameter_a_index25.mat';%NS �����ĳ�ֵ
% load 'NS_ini_parameter_a1.mat';%֮ǰָ��Ƶ����ϵĵ�����ֵ


tic;
% % parameter_a(2:end,:)=-sqrt(parameter_a(2:end,:).^2/2);
% parameter_a(1,:)=[w1,wd,0,0];
% parameter_a(2:end,:)=parameter_a1(1:N_harm,:);
% parameter_a(2,1)=0.5;
%%
ini_parameter_a=parameter_a;
iteration=length(Tdata);
% to indicate where a is in the parametric space
parameter_a_judge=@(parameter_a)(abs(parameter_a(2,:))<1);
gammaT=1.414;rhob=0.5; % parameter for trust-region algorithm
parameter_a_record=parameter_a(1:N_harm+1,1:2*N_dof); % To record the values of parameters during iteration, and for each iteration, a is recorded in a single row of a_record
TR_record=[];  % recording the parameters during trust region
%% �����������ٶȺͼ��ٶ���Ӧʶ�����������ɶ������ȼ���ͬ���߼��ٶ�1%���Ǽ��ٶ�2%��5%
%% Response sensitivity iteration
Nmax=1000;   % maximum number for response sensitivity iteration
Ntr=20;      % maximum number for trust region iteration
%% response sensitivity Solution by ode45
for iii=1:Nmax
    %% ��������w0,�˴�Ϊһ�������ھ���ѡȡ1K����
    % compute response and response sensitivity for each incremental
    Etol=1e-10;  % Relative error tolerance for convergence of the RS algorithm
    %%
    residual_iden=cal_residual2(parameter_a);
    %% SSSΪλ����Ӧ�����Ⱦ��󣬵�һ�к͵ڶ���Ϊ�в���Ӧ��Ƶ�������ȴӵ����п�ʼ
    % ����������Ⱦ����а���S_11����������Ҫȥ��
    SSS=reshape(residual_iden(:,R_dof+1:2*R_dof),R_dof*length(Tdata),1);%w01
    for i=1:2*N_harm*N_dof-1+N_w0
        SSS=[SSS,reshape(residual_iden(:,R_dof*(i+1)+1:R_dof*(i+2)),R_dof*length(Tdata),1)];
    end
    
    dR=-reshape(residual_iden(:,1:R_dof),R_dof*length(Tdata),1);
    [U,s,V]=csvd(SSS);
    lambda_inverse=l_curve(U,s,dR);
    atemp=parameter_a(1:N_harm+1,1:2*N_dof);
    % trust-region algorithm
    for trust=1:Ntr
        %% �����da����Ϊw0,C_11,S_11(ȱ,���������),C_12,S_12,���ղ����������ʽ����da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        temp_real_w0=zeros(1,2*N_dof);% ��������ĵ�һ��
        temp_real_w0(1,1:N_w0)=real_da(1:N_w0,1)';
        temp_real_da=real_da(N_w0+1:end,1);
        da=reshape(temp_real_da,2,N_dof*N_harm);da=da';
        sensitivity_parameter_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
        end
        sensitivity_parameter_da=[temp_real_w0;sensitivity_parameter_da];
        
        %%
        if ~parameter_a_judge(atemp+sensitivity_parameter_da)      % if updated a is not out of the parametric space, then, lambda should be increased until updated a is in ...
            lambda_inverse=lambda_inverse*gammaT; %  update of lambda
            continue;
        end
        %% �����da����Ϊw0,C_11,S_11(ȱ,���������),C_12,S_12,���ղ����������ʽ����da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        temp_real_w0=zeros(1,2*N_dof);% ��������ĵ�һ��
        temp_real_w0(1,1:N_w0)=real_da(1:N_w0,1)';
        temp_real_da=real_da(N_w0+1:end,1);
        da=reshape(temp_real_da,2,N_dof*N_harm);da=da';
        sensitivity_parameter_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            sensitivity_parameter_da=[sensitivity_parameter_da,da(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
        end
        sensitivity_parameter_da=[temp_real_w0;sensitivity_parameter_da];
        
        parameter_a=atemp+sensitivity_parameter_da;

        %         N_real_harm=ceil(N_min+(N_harm-N_min)*exp(-alpha*(iii-1)));
        %         Harm_parameter_a=parameter_a(2:N_harm+1,:);
        %         %% ����ÿһ��Ƶ�ʵ����
        %         for j=1:N_dof
        %             for i=1:N_harm
        %                 vector_amplitude(i,j)=sqrt(Harm_parameter_a(i,2*j-1)^2+Harm_parameter_a(i,2*j)^2);
        %             end
        %         end
        %         % ��������н������У�����¼����֮ǰ��λ��
        %         % descend����ascend����positionԭ��Ϊλ�ã�
        %         [descend_vector_amplitude,position]=sort(vector_amplitude,1,'descend');
        %         temp_zeros=zeros(size(Harm_parameter_a));
        %         for j=1:N_dof
        %             for i=1:N_real_harm
        %                 temp_zeros(position(i,j),2*j-1)=1;
        %                 temp_zeros(position(i,j),2*j)=1;
        %             end
        %         end
        %         %         temp_zeros(4,1:2)=zeros(1,2);temp_zeros(5,1:2)=zeros(1,2);temp_zeros(9,1:2)=zeros(1,2);
        %         %         temp_zeros(11,1:2)=zeros(1,2);temp_zeros(13:32,1:2)=0;
        %         parameter_a(2:N_harm+1,:)=temp_zeros.*parameter_a(2:N_harm+1,:);
        
        %         sensitivity_parameter_da=sensitivity_parameter_da(2:end,:);
        %         sensitivity_parameter_da=temp_zeros.*sensitivity_parameter_da;
        
        ini_da=real_da;

        %% ���µ�parameter_a=atemp+da������Ӧ
        residual_da=cal_residual2(parameter_a);
        dRtemp=-reshape(residual_da(:,1:R_dof),R_dof*length(Tdata),1);
        LdR=SSS*ini_da-dR;
        rhos=(dR'*dR-dRtemp'*dRtemp)/(dR'*dR-LdR'*LdR);  % agreement indicator
        if rhos>=rhob
            break;
        end
        lambda_inverse=lambda_inverse*gammaT;
    end
    tolt=norm(da)/norm(parameter_a);
    parameter_a_record=[parameter_a_record,parameter_a];
    TR_record=[TR_record;lambda_inverse];
    parameter_a
    if tolt<=Etol
        break;
    end
    every_a(iii).parameter_a=parameter_a;
    iii
end
toc;

% parameter_a(1,1:N_w0)=[w1,wd];
residua=cal_residual(parameter_a);

% figure;
% plot(Tdata,residua(:,1),'r-','LineWidth',1);
% hold on;
% plot(Tdata,residua(:,2),'k-','LineWidth',1);
% h1=legend('$$x_1$$','$$x_2$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
subplot(2,1,1);
plot(Tdata,residua(:,1),'r-','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
subplot(2,1,2);
plot(Tdata,residua(:,2),'b-','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


Tdata=0:0.004:30;
% ��X,Y�������ɶ�
%% ����������Ƶ�ʵ���ϣ�ע�⣬Ƶ�����Ҫȥ������Ƶ��
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
% ���������ɶ�
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
        dx(j,:)=dx(j,:)-vector_w(i)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i)*Tdata)+vector_w(i)*Harm_parameter_a(i,2*j)*cos(vector_w(i)*Tdata);
        ddx(j,:)=ddx(j,:)-(vector_w(i))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)-(vector_w(i))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
    end
end

load matlab.mat;
figure;
plot(Tdata,x(1,:),'r-','LineWidth',1);
hold on;
plot(Tdata,x(2,:),'k-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
plot(Tdata,num(:,1),'k-','LineWidth',1);
hold on;
plot(Tdata,num(:,2),'r-','LineWidth',1);
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
hold on;
plot(num(:,1),num(:,3),'k-','LineWidth',1);
hold on;
plot(num(:,2),num(:,4),'r-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
ini_parameter_a

% figure;
% plot(Tdata,x(1,:)'-num(:,1),'r-','LineWidth',1);
% hold on;
% plot(Tdata,x(2,:)'-num(:,2),'k-','LineWidth',1);
% h1=legend('$$x_1$$','$$x_2$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
subplot(2,1,1);
plot(Tdata,x(1,:)'-num(:,1),'r-','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
subplot(2,1,2);
plot(Tdata,x(2,:)'-num(:,2),'b-','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



figure;
subplot(2,1,1);
plot(Tdata,x(1,:),'r-','LineWidth',1);
hold on;
plot(Tdata,num(:,1),'k--','LineWidth',1);
h1=legend('$$WTMRM$$','$$RK$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
subplot(2,1,2);
plot(Tdata,x(2,:),'b-','LineWidth',1);
hold on;
plot(Tdata(1:5:end),num(1:5:end,2),'k.','LineWidth',1);
h1=legend('$$WTMRM$$','$$f$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);









