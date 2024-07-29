clear;clc;close all;
load 'all_results.mat';%NS 迭代的初值
% for j=1:N_dof
%     for i=1:N_harm   % i=1,3,5
%         A_x(j,i)=sqrt(Harm_parameter_a(i,2*j-1)^2+Harm_parameter_a(i,2*j)^2);
%     end
% end
% figure;
% plot(vector_w,A_x(1,:),'k.','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
% plot(vector_w,A_x(2,:),'k.','LineWidth',1);
% % h1=legend('$$x_1$$','$$x_2$$');
% % set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% 
% figure;
% subplot(2,1,1);
% plot(vector_w,A_x(1,:),'k.','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% subplot(2,1,2);
% plot(vector_w,A_x(2,:),'k.','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

load matlab.mat;


figure;
dot1=60;dot2=16;
subplot(2,1,1);
plot(Tdata,x(1,:),'r-','LineWidth',1.5);
hold on;
plot(Tdata(1:dot1:end),num(1:dot1:2000001,1),'k.','MarkerSize',8);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
subplot(2,1,2);
plot(Tdata,x(2,:),'b-','LineWidth',1.5);
hold on;
plot(Tdata(1:dot2:end),num(1:dot2:2000001,2),'k.','MarkerSize',8);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);





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
h1=legend('$$N=5N=15N=25$$','$$x_2$$');
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

figure;
plot(Tdata,x(1,:)'-num(:,1),'r-','LineWidth',1);
hold on;
plot(Tdata,x(2,:)'-num(:,2),'k-','LineWidth',1);
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

