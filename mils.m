%多新息递推最小二乘参数估计（mils）
clear all; close all;

a=[1 -1.5 0.7]'; b=[1 0.5]'; d=1; %对象参数
na=length(a)-1; nb=length(b)-1; %na、nb为A、B阶次

L=400; %仿真长度
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
u=3*randn(L,1); %输入采用白噪声序列
xi=sqrt(0.1)*randn(L,1); %白噪声序列

theta=[a(2:na+1);b]; %对象参数真值

p=10;                           %插入新息维数
thetae_mils1=zeros(na+nb+1,1); %mils算法thetae初值
P=10^6*eye(na+nb+1);           %mils算法协方差矩阵P初值
phi=0.0*ones(na+nb+1,p);       %mils算法构造多新息的系统输入输出信息阵phi
Y=0.0*ones(p,1);               %mils算法构造多新息系统输出向量
Ip=eye(p);
s=0.0*ones(p,p);
se=0.0*ones(p,p);
for k=1:L
    h=[-yk;uk(1:nb+1)];
    
    for i=p:-1:2
        phi(:,i)=phi(:,i-1);    %构造多新息观测矩阵，第k步时，phi(p,t)第i列为上一步前一列值
    end
    phi(:,1)=h;                 %phi(p,t)第1列为系统信息向量h的值

    
    y(k)=h'*theta+xi(k);        %采集系统输出真值数据
    
    for i=p:-1:2
        Y(i,1)=Y(i-1,1);        %构造多新息系统输出向量Y（p,t),第k步时，Y(p,t)第i行为上一步        
    end                         %上一行的值
    Y(1,:)=y(k);                %Y(p,t)第一行为该时刻观测到的输出y的值
   
    %多新息递推最小二乘法
    s=(Ip+phi'*P*phi);
    se=inv(s);
    K=P*phi*se;
    thetae(:,k)=thetae_mils1+K*(Y-phi'*thetae_mils1);
    P=P-P*phi*se*phi'*P;
    
    %更新数据
    thetae_mils1=thetae(:,k);
    
    for i=d+nb:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=na:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
end
plot([1:L],thetae); %line([1,L],[theta,theta]);
xlabel('k'); ylabel('MILS参数估计a、b');
legend('a_1','a_2','b_0','b_1'); axis([0 L -2 2]);