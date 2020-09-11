%递推最小二乘（ls）与多新息递推最小二乘（mils）参数估计比较
clear all; close all;

a=[1 -1.5 0.7]'; b=[1 0.5]'; d=1; %对象参数
na=length(a)-1; nb=length(b)-1; %na、nb为A、B阶次

L=400; %仿真长度
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
u=3*randn(L,1); %输入采用白噪声序列
xi=sqrt(0.1)*randn(L,1); %白噪声序列

theta=[a(2:na+1);b]; %对象参数真值

p=10;
thetae_mils1=zeros(na+nb+1,1); %mils算法thetae初值
P_mils=10^6*eye(na+nb+1);      %mils算法数据协方差矩阵P初值
phi=0.0*ones(na+nb+1,p);       %mils算法构造多新息的系统输入输出观测阵phi(p,t)
Y=0.0*ones(p,1);               %mils算法构造多新息系统输出向量y（y,t）
Ip=eye(p);
s=0.0*ones(p,p);
se=0.0*ones(p,p);
for k=1:L
    h=[-yk;uk(1:nb+1)];
    for i=p:-1:2
        phi(:,i)=phi(:,i-1);    %构造多新息观测矩阵
    end
    phi(:,1)=h;
    y(k)=h'*theta+xi(k);        %采集输出数据
    for i=p:-1:2
        Y(i,1)=Y(i-1,1);
    end
    Y(1,:)=y(k);
   
    %多新息递推最小二乘法、递推最小二乘法
    s=(Ip+phi'*P_mils*phi);
    se=inv(s);
    K_mils=P_mils*phi*se;
    thetae_mils(:,k)=thetae_mils1+K_mils*(Y-phi'*thetae_mils1);
    P_mils=P_mils-P_mils*phi*se*phi'*P_mils;
    
    %计算估计误差
    e_mils1=theta-thetae_mils(:,k);
    e_mils=norm(e_mils1,2)/norm(theta,2);    %t时刻mils算法参数估值误差
    e_m(k)=e_mils;
    
    %更新数据
    thetae_mils1=thetae_mils(:,k);
    
    for i=d+nb:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=na:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
end
figure;
plot(1:L,thetae_mils); %line([1,L],[theta,theta]);
xlabel('k');ylabel('MILS参数估计a、b');
legend('a_1','a_2','b_0','b_1'); axis([0 L -2 2]);


uk_ls=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk_ls=zeros(na,1); %输出初值
u_ls=3*randn(L,1); %输入采用白噪声序列
xi_ls=sqrt(0.1)*randn(L,1); %白噪声序列


thetae_1=zeros(na+nb+1,1); %thetae初值
P=10^6*eye(na+nb+1); 
for k=1:L
    phi_ls=[-yk_ls;uk_ls(d:d+nb)]; %此处phi为列向量
    y_ls(k)=phi_ls'*theta+xi_ls(k); %采集输出数据
   
    %递推最小二乘法
    K=P*phi_ls/(1+phi_ls'*P*phi_ls);
    thetae(:,k)=thetae_1+K*(y_ls(k)-phi_ls'*thetae_1);
    P=(eye(na+nb+1)-K*phi_ls')*P;
    
    %更新数据
    thetae_1=thetae(:,k);
    
    for i=d+nb:-1:2
        uk_ls(i)=uk_ls(i-1);
    end
    uk_ls(1)=u_ls(k);
    
    for i=na:-1:2
        yk_ls(i)=yk_ls(i-1);
    end
    yk_ls(1)=y_ls(k);
    
    e_ls1=theta-thetae(:,k);
    e_ls=norm(e_ls1,2)/norm(theta,2);     %t时刻ls算法参数估计误差
    e_l(k)=e_ls;
end

figure;
plot([1:L],thetae);
xlabel('k'); ylabel('LS参数估计a、b');
legend('a_1','a_2','b_0','b_1'); axis([0 L -2 2]);
figure;
plot([1:L],e_m,'--b',[1:L],e_l,'-r');
legend('e mils','e ls');
axis([0 L 0 0.1]);
