function residual=cal_residual(parameter_a)
global index_global N_dof Tdata N_harm N_w0 index %w
amplification=1;
Harm_parameter_a=parameter_a(2:end,:);
%% 自由度数目，谐波数目，参数个数2*N_harm*N_dof
%系统参数
N_dof1=6;
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
M=[1,0;0,1];
K=[k11,k12;k21,k22];

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
%此时，vector_w的频率和h_mat完全对应


%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有三个自由度
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
        dx(j,:)=dx(j,:)-vector_w(i)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i)*Tdata)+vector_w(i)*Harm_parameter_a(i,2*j)*cos(vector_w(i)*Tdata);
        ddx(j,:)=ddx(j,:)-(vector_w(i))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)-(vector_w(i))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
    end
end
non1=k13*x(1,:).^3+k14*x(1,:).^2.*x(2,:)+k15*x(1,:).*x(2,:).^2+k16*x(2,:).^3+...
    c11*x(1,:).*dx(1,:).^2+c12*x(2,:).*dx(1,:).^2+c13*x(1,:).*dx(1,:).*dx(2,:)+c14*x(2,:).*dx(1,:).*dx(2,:)+c15*x(1,:).*dx(2,:).^2+c16*x(2,:).*dx(2,:).^2;
non2=k23*x(1,:).^3+k24*x(1,:).^2.*x(2,:)+k25*x(1,:).*x(2,:).^2+k26*x(2,:).^3+...
    c21*x(1,:).*dx(1,:).^2+c22*x(2,:).*dx(1,:).^2+c23*x(1,:).*dx(1,:).*dx(2,:)+c24*x(2,:).*dx(1,:).*dx(2,:)+c25*x(1,:).*dx(2,:).^2+c26*x(2,:).*dx(2,:).^2;
residual(1:N_dof,:)=M*ddx+K*x+[non1;non2];

residual(N_dof+1,:)=amplification*(sum(Harm_parameter_a(:,1))-A1);
residual(N_dof+2,:)=amplification*(sum(Harm_parameter_a(:,3))-A2);
residual(N_dof+3,:)=amplification*sum(vector_w(:,1).*Harm_parameter_a(:,2));
residual(N_dof+4,:)=amplification*sum(vector_w(:,1).*Harm_parameter_a(:,4));
%% 计算频率的灵敏度  只需要计算基频的就可以
%生成数组
for i=1:2:index
    temp_index_global_w1((i+1)/2,1:4)=i;
end
size_temp_w1=size(temp_index_global_w1);
index_global_w1=[];
for i=1:size_temp_w1(1,1)
    index_global_w1=[index_global_w1;temp_index_global_w1(i,:)'];
end
%此处vector_w1为一个列向量
%% 计算频率的灵敏度  只需要计算基频w1的就可以
x_w1=zeros(N_dof,length(Tdata));dx_w1=zeros(N_dof,length(Tdata));ddx_w1=zeros(N_dof,length(Tdata));
for ij=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x_w1(ij,:)=x_w1(ij,:)-index_global_w1(i)*Harm_parameter_a(i,2*ij-1)*Tdata.*sin(vector_w(i)*Tdata)+index_global_w1(i)*Harm_parameter_a(i,2*ij)*Tdata.*cos(vector_w(i)*Tdata);
        dx_w1(ij,:)=dx_w1(ij,:)-(index_global_w1(i)*Harm_parameter_a(i,2*ij-1)*sin(vector_w(i)*Tdata)+vector_w(i)*index_global_w1(i)*Tdata*Harm_parameter_a(i,2*ij-1).*cos(vector_w(i)*Tdata))+...
            (index_global_w1(i)*Harm_parameter_a(i,2*ij)*cos(vector_w(i)*Tdata)-vector_w(i)*index_global_w1(i)*Tdata*Harm_parameter_a(i,2*ij).*sin(vector_w(i)*Tdata));
        ddx_w1(ij,:)=ddx_w1(ij,:)-(2*vector_w(i)*index_global_w1(i)*Harm_parameter_a(i,2*ij-1)*cos(vector_w(i)*Tdata)-vector_w(i)^2*index_global_w1(i)*Tdata*Harm_parameter_a(i,2*ij-1).*sin(vector_w(i)*Tdata))-...
            (2*vector_w(i)*index_global_w1(i)*Harm_parameter_a(i,2*ij)*sin(vector_w(i)*Tdata)+vector_w(i)^2*index_global_w1(i)*Tdata*Harm_parameter_a(i,2*ij).*cos(vector_w(i)*Tdata));
    end
end
non1_w1=3*k13*x(1,:).^2.*x_w1(1,:)+k14*(x(1,:).^2.*x_w1(2,:)+2*x(1,:).*x(2,:).*x_w1(1,:))+...
    k15*(x_w1(1,:).*x(2,:).^2+2*x(1,:).*x(2,:).*x_w1(2,:))+3*k16*x(2,:).^2.*x_w1(2,:)+...
    c11*(x(1,:)*2.*dx(1,:).*dx_w1(1,:)+x_w1(1,:).*dx(1,:).^2)+c12*(x(2,:)*2.*dx(1,:).*dx_w1(1,:)+x_w1(2,:).*dx(1,:).^2)+c13*(x_w1(1,:).*dx(1,:).*dx(2,:)+x(1,:).*dx(1,:).*dx_w1(2,:)+x(1,:).*dx(2,:).*dx_w1(1,:))...
    +c14*(x_w1(2,:).*dx(1,:).*dx(2,:)+x(2,:).*dx(1,:).*dx_w1(2,:)+x(2,:).*dx(2,:).*dx_w1(1,:))+c15*(x(1,:)*2.*dx(2,:).*dx_w1(2,:)+x_w1(1,:).*dx(2,:).^2)+c16*(x(2,:)*2.*dx(2,:).*dx_w1(2,:)+x_w1(2,:).*dx(2,:).^2);
non2_w1=3*k23*x(1,:).^2.*x_w1(1,:)+k24*(x(1,:).^2.*x_w1(2,:)+2*x(1,:).*x(2,:).*x_w1(1,:))+...
    k25*(x_w1(1,:).*x(2,:).^2+2*x(1,:).*x(2,:).*x_w1(2,:))+3*k26*x(2,:).^2.*x_w1(2,:)+...
    c21*(x(1,:)*2.*dx(1,:).*dx_w1(1,:)+x_w1(1,:).*dx(1,:).^2)+c22*(x(2,:)*2.*dx(1,:).*dx_w1(1,:)+x_w1(2,:).*dx(1,:).^2)+c23*(x_w1(1,:).*dx(1,:).*dx(2,:)+x(1,:).*dx(1,:).*dx_w1(2,:)+x(1,:).*dx(2,:).*dx_w1(1,:))...
    +c24*(x_w1(2,:).*dx(1,:).*dx(2,:)+x(2,:).*dx(1,:).*dx_w1(2,:)+x(2,:).*dx(2,:).*dx_w1(1,:))+c25*(x(1,:)*2.*dx(2,:).*dx_w1(2,:)+x_w1(1,:).*dx(2,:).^2)+c26*(x(2,:)*2.*dx(2,:).*dx_w1(2,:)+x_w1(2,:).*dx(2,:).^2);
residual(N_dof1+1:2*N_dof1-4,:)=M*ddx_w1+K*x_w1+[non1_w1;non2_w1];
residual(2*N_dof1-3:2*N_dof1-2,:)=0;
residual(2*N_dof1-1,:)=amplification*sum(index_global_w1(:,1).*Harm_parameter_a(:,2));
residual(2*N_dof1,:)=amplification*sum(index_global_w1(:,1).*Harm_parameter_a(:,4));

%% 计算频率的灵敏度  只需要计算频宽的就可以
%生成数组
for i=1:2:index
    temp_index_global_wd((i+1)/2,1)=0;
    temp_index_global_wd((i+1)/2,2:4)=[-1,1,2];
end
size_temp_wd=size(temp_index_global_wd);
index_global_wd=[];
for i=1:size_temp_wd(1,1)
    index_global_wd=[index_global_wd;temp_index_global_wd(i,:)'];
end
%此处vector_wd为一个列向量
%% 计算频率的灵敏度  只需要计算频宽wd的就可以
x_wd=zeros(N_dof,length(Tdata));dx_wd=zeros(N_dof,length(Tdata));ddx_wd=zeros(N_dof,length(Tdata));
for ij=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x_wd(ij,:)=x_wd(ij,:)-index_global_wd(i)*Harm_parameter_a(i,2*ij-1)*Tdata.*sin(vector_w(i)*Tdata)+index_global_wd(i)*Harm_parameter_a(i,2*ij)*Tdata.*cos(vector_w(i)*Tdata);
        dx_wd(ij,:)=dx_wd(ij,:)-(index_global_wd(i)*Harm_parameter_a(i,2*ij-1)*sin(vector_w(i)*Tdata)+vector_w(i)*index_global_wd(i)*Tdata*Harm_parameter_a(i,2*ij-1).*cos(vector_w(i)*Tdata))+...
            (index_global_wd(i)*Harm_parameter_a(i,2*ij)*cos(vector_w(i)*Tdata)-vector_w(i)*index_global_wd(i)*Tdata*Harm_parameter_a(i,2*ij).*sin(vector_w(i)*Tdata));
        ddx_wd(ij,:)=ddx_wd(ij,:)-(2*vector_w(i)*index_global_wd(i)*Harm_parameter_a(i,2*ij-1)*cos(vector_w(i)*Tdata)-vector_w(i)^2*index_global_wd(i)*Tdata*Harm_parameter_a(i,2*ij-1).*sin(vector_w(i)*Tdata))-...
            (2*vector_w(i)*index_global_wd(i)*Harm_parameter_a(i,2*ij)*sin(vector_w(i)*Tdata)+vector_w(i)^2*index_global_wd(i)*Tdata*Harm_parameter_a(i,2*ij).*cos(vector_w(i)*Tdata));
    end
end
non1_wd=3*k13*x(1,:).^2.*x_wd(1,:)+k14*(x(1,:).^2.*x_wd(2,:)+2*x(1,:).*x(2,:).*x_wd(1,:))+...
    k15*(x_wd(1,:).*x(2,:).^2+2*x(1,:).*x(2,:).*x_wd(2,:))+3*k16*x(2,:).^2.*x_wd(2,:)+...
    c11*(x(1,:)*2.*dx(1,:).*dx_wd(1,:)+x_wd(1,:).*dx(1,:).^2)+c12*(x(2,:)*2.*dx(1,:).*dx_wd(1,:)+x_wd(2,:).*dx(1,:).^2)+c13*(x_wd(1,:).*dx(1,:).*dx(2,:)+x(1,:).*dx(1,:).*dx_wd(2,:)+x(1,:).*dx(2,:).*dx_wd(1,:))...
    +c14*(x_wd(2,:).*dx(1,:).*dx(2,:)+x(2,:).*dx(1,:).*dx_wd(2,:)+x(2,:).*dx(2,:).*dx_wd(1,:))+c15*(x(1,:)*2.*dx(2,:).*dx_wd(2,:)+x_wd(1,:).*dx(2,:).^2)+c16*(x(2,:)*2.*dx(2,:).*dx_wd(2,:)+x_wd(2,:).*dx(2,:).^2);
non2_wd=3*k23*x(1,:).^2.*x_wd(1,:)+k24*(x(1,:).^2.*x_wd(2,:)+2*x(1,:).*x(2,:).*x_wd(1,:))+...
    k25*(x_wd(1,:).*x(2,:).^2+2*x(1,:).*x(2,:).*x_wd(2,:))+3*k26*x(2,:).^2.*x_wd(2,:)+...
    c21*(x(1,:)*2.*dx(1,:).*dx_wd(1,:)+x_wd(1,:).*dx(1,:).^2)+c22*(x(2,:)*2.*dx(1,:).*dx_wd(1,:)+x_wd(2,:).*dx(1,:).^2)+c23*(x_wd(1,:).*dx(1,:).*dx(2,:)+x(1,:).*dx(1,:).*dx_wd(2,:)+x(1,:).*dx(2,:).*dx_wd(1,:))...
    +c24*(x_wd(2,:).*dx(1,:).*dx(2,:)+x(2,:).*dx(1,:).*dx_wd(2,:)+x(2,:).*dx(2,:).*dx_wd(1,:))+c25*(x(1,:)*2.*dx(2,:).*dx_wd(2,:)+x_wd(1,:).*dx(2,:).^2)+c26*(x(2,:)*2.*dx(2,:).*dx_wd(2,:)+x_wd(2,:).*dx(2,:).^2);
residual(2*N_dof1+1:3*N_dof1-4,:)=M*ddx_wd+K*x_wd+[non1_wd;non2_wd];
residual(3*N_dof1-3:3*N_dof1-2,:)=0;
residual(3*N_dof1-1,:)=amplification*sum(index_global_wd(:,1).*Harm_parameter_a(:,2));
residual(3*N_dof1,:)=amplification*sum(index_global_wd(:,1).*Harm_parameter_a(:,4));

%% 计算谐波系数的灵敏度
for i=1:2*N_harm*N_dof
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    x_a=zeros(N_dof,length(Tdata));dx_a=zeros(N_dof,length(Tdata));ddx_a=zeros(N_dof,length(Tdata));
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    for k=1:N_dof
        for j=1:N_harm
            x_a(k,:)=x_a(k,:)+sensitivity_parameter_a(j,2*k-1)*cos(vector_w(j)*Tdata)+sensitivity_parameter_a(j,2*k)*sin(vector_w(j)*Tdata);
            dx_a(k,:)=dx_a(k,:)-vector_w(j)*sensitivity_parameter_a(j,2*k-1)*sin(vector_w(j)*Tdata)+vector_w(j)*sensitivity_parameter_a(j,2*k)*cos(vector_w(j)*Tdata);
            ddx_a(k,:)=ddx_a(k,:)-(vector_w(j))^2*sensitivity_parameter_a(j,2*k-1)*cos(vector_w(j)*Tdata)-(vector_w(j))^2*sensitivity_parameter_a(j,2*k)*sin(vector_w(j)*Tdata);
        end
    end
    non1_a=3*k13*x(1,:).^2.*x_a(1,:)+k14*(x(1,:).^2.*x_a(2,:)+2*x(1,:).*x(2,:).*x_a(1,:))+...
        k15*(x_a(1,:).*x(2,:).^2+2*x(1,:).*x(2,:).*x_a(2,:))+3*k16*x(2,:).^2.*x_a(2,:)+...
        c11*(x(1,:)*2.*dx(1,:).*dx_a(1,:)+x_a(1,:).*dx(1,:).^2)+c12*(x(2,:)*2.*dx(1,:).*dx_a(1,:)+x_a(2,:).*dx(1,:).^2)+c13*(x_a(1,:).*dx(1,:).*dx(2,:)+x(1,:).*dx(1,:).*dx_a(2,:)+x(1,:).*dx(2,:).*dx_a(1,:))...
        +c14*(x_a(2,:).*dx(1,:).*dx(2,:)+x(2,:).*dx(1,:).*dx_a(2,:)+x(2,:).*dx(2,:).*dx_a(1,:))+c15*(x(1,:)*2.*dx(2,:).*dx_a(2,:)+x_a(1,:).*dx(2,:).^2)+c16*(x(2,:)*2.*dx(2,:).*dx_a(2,:)+x_a(2,:).*dx(2,:).^2);
    non2_a=3*k23*x(1,:).^2.*x_a(1,:)+k24*(x(1,:).^2.*x_a(2,:)+2*x(1,:).*x(2,:).*x_a(1,:))+...
        k25*(x_a(1,:).*x(2,:).^2+2*x(1,:).*x(2,:).*x_a(2,:))+3*k26*x(2,:).^2.*x_a(2,:)+...
        c21*(x(1,:)*2.*dx(1,:).*dx_a(1,:)+x_a(1,:).*dx(1,:).^2)+c22*(x(2,:)*2.*dx(1,:).*dx_a(1,:)+x_a(2,:).*dx(1,:).^2)+c23*(x_a(1,:).*dx(1,:).*dx(2,:)+x(1,:).*dx(1,:).*dx_a(2,:)+x(1,:).*dx(2,:).*dx_a(1,:))...
        +c24*(x_a(2,:).*dx(1,:).*dx(2,:)+x(2,:).*dx(1,:).*dx_a(2,:)+x(2,:).*dx(2,:).*dx_a(1,:))+c25*(x(1,:)*2.*dx(2,:).*dx_a(2,:)+x_a(1,:).*dx(2,:).^2)+c26*(x(2,:)*2.*dx(2,:).*dx_a(2,:)+x_a(2,:).*dx(2,:).^2);
    residual(N_dof1*(i+N_w0)+1:N_dof1*(i+N_w0+1)-4,:)=M*ddx_a+K*x_a+[non1_a;non2_a];
    if i<N_harm*N_dof
        residual(N_dof1*(i+N_w0+1)-3,:)=amplification*mod(i,2);
        residual(N_dof1*(i+N_w0+1)-2,:)=0;
    else
        residual(N_dof1*(i+N_w0+1)-3,:)=0;
        residual(N_dof1*(i+N_w0+1)-2,:)=amplification*mod(i,2);
    end
    if mod(i,2)==1
        residual(N_dof1*(i+N_w0+1)-1,:)=0;
        residual(N_dof1*(i+N_w0+1),:)=0;
    else
        residual(N_dof1*(i+N_w0+1)-1,:)=amplification*sum(sum(sensitivity_parameter_a(:,1:2).*[vector_w(:,1),vector_w(:,1)]));
        residual(N_dof1*(i+N_w0+1),:)=amplification*sum(sum(sensitivity_parameter_a(:,3:4).*[vector_w(:,1),vector_w(:,1)]));
    end
end

residual=residual';

