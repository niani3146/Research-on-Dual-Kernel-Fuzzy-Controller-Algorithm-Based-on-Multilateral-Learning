% %二自由度机械臂-双边双核模糊-模糊规则1
% %以各和函数输出作为自己的参数更新依据
% clear all;
% close all;
N=2;%选取的边数
n=3;
m=2;

t0=0;tf=30;
T=linspace(t0,tf,100);

para_a=linspace(0,0.6,N);
para_b=linspace(0,0.7,N);
% para_a=zeros(1,N);
% para_b=zeros(1,N);
para_k=zeros(1,N*m*(n^4));
para_w=zeros(1,N*m*(n^4));
para_n=zeros(1,N*m*(n^4));
% para_m=zeros(1,n^4);
para_c=zeros(1,2);          
para_d=zeros(1,N*m);
para_m=zeros(1,8);

x0=[1,1,1,1,1,1,1,1,0,0,0,0,... 
  para_a,para_b,para_k,para_w,para_n,para_c,para_d,para_m];

options = odeset('MaxStep', 0.1); % 设置最大步长
% options = odeset('MaxStep', 0.02); % 设置最大步长
[t,x]=ode45(@mohu_2_fn_20250228,[t0,tf],x0,options);

x_ref1=[sin(t),cos(t)];
x_ref2=[cos(t),-sin(t)];
% x_ref1=zeros(length(t),2); 
% x_ref2=zeros(length(t),2);
% for i=1:length(t)
%     if t(i)>=20 && t(i)<50
%         x_ref1(i,1)=1;x_ref2(i,1)=2;
%     elseif t(i)>50 && t(i)<65
%         x_ref1(i,1)=3;x_ref2(i,1)=1.5;
%     else
%         x_ref1(i,1)=0;x_ref2(i,1)=0;
%     end
% end

%RR-bot系统参数
a0=1;%连杆1长度为1
a1=1;%连杆2长度为1
m1=0.1;%连杆1的质量为0.1kg
m2=0.1;%连杆2的质量为0.2kg
gs=9.8;%重力加速度
y_true=zeros(2,length(t));
for i=1:length(t)
    M=[(a0^2*m1)/3 + a0^2*m2 + (a1^2*m2)/3 + a0*a1*m2*cos(x(i,2)),  (a1*m2*(2*a1 + 3*a0*cos(x(i,2))))/6;
        (a1*m2*(2*a1 + 3*a0*cos(x(i,2))))/6,  (a1^2*m2)/3];
    C=[-(a0*a1*x(i,4)*m2*sin(x(i,2)))/2,   -(a0*a1*m2*sin(x(i,2))*(x(i,3) + x(i,4)))/2;
        (a0*a1*x(i,3)*m2*sin(x(i,2)))/2,   0];
    G=[gs*m2*((a1*cos(x(i,1) + x(i,2)))/2 + a0*cos(x(i,1))) + (a0*gs*m1*cos(x(i,1)))/2;
        (a1*gs*m2*cos(x(i,1) + x(i,2)))/2];

    if abs(det(M))<1e-3
        pseudo_invM=M'*M\M';
        un_model=-pseudo_invM*(C*[x(i,3);x(i,4)]+G);
    else
        un_model=-inv(M)*(C*[x(i,3);x(i,4)]+G);
    end
    y_true(:,i)=[un_model(1);un_model(2)];
end

kr=[-1,0,-2,0,1,0;...
    0,-1,0,-2,0,1]';
xre=[x(:,5),x(:,6),x(:,7),x(:,8),x_ref1(:,1),x_ref2(:,1)];%length(t)*6   增广状态

ke=[50,0,20,0;...
    0,50,0,20]';
er11=x_ref1(:,1)-x(:,5);
er12=x_ref1(:,2)-x(:,7);
er21=x_ref2(:,1)-x(:,6);
er22=x_ref2(:,2)-x(:,8);%er

e11=x(:,5)-x(:,1);
e12=x(:,7)-x(:,3);
e1=[e11;e12];
e21=x(:,6)-x(:,2);
e22=x(:,8)-x(:,4);
e2=[e21;e22];
e=[e11,e21,e12,e22];%e=[e,e']

y_ref1=y_true(1,:)';
y_ref2=y_true(2,:)';

u_tau=zeros(2,length(t));

%模糊控制
% 每个模糊控制中包含两个隶属度函数，经过加权融合的输出结果
xi_12= x(12+2*N+(n^4)*8+1:12+2*N+(n^4)*8+n^4);
xi_22=x(12+2*N+(n^4)*9+1:12+2*N+(n^4)*9+n^4);
xi_32=x(12+2*N+(n^4)*10+1:12+2*N+(n^4)*10+n^4);
xi_42=x(12+2*N+(n^4)*11+1:12+2*N+(n^4)*11+n^4);
% 第一条边的构成

ee_yn_1=x(:,12+2*N+(n^4)*12+2+1:12+2*N+(n^4)*12+2+N*m);
y_hat_m11=zeros(length(t),1);
W_hat_m11=zeros(length(t),n^4);
W_hat_m12=zeros(length(t),n^4);
W_hat_m21=zeros(length(t),n^4);
W_hat_m22=zeros(length(t),n^4);
for j=1:length(t)
            W_hat_m11(j,:)=x(j,12+2*N+1:12+2*N+n^4);
             y_hat_m11(j)= W_hat_m11(j,:)* xi_12';                %高斯 1
end
y_hat_m12=zeros(length(t),1);
for j=1:length(t)
            W_hat_m12(j,:)=x(j,12+2*N+(n^4)*4+1:12+2*N+(n^4)*4+n^4);
             y_hat_m12(j)= W_hat_m12(j,:)* xi_12';              %高斯 2
end  

y_hat_m21=zeros(length(t),1);
for j=1:length(t)
            W_hat_m21(j,:)=x(j,12+2*N+(n^4)+1:12+2*N+(n^4)+n^4);
             y_hat_m21(j)= W_hat_m21(j,:)* xi_22';               %钟型 1
end
y_hat_m22=zeros(length(t),1);
for j=1:length(t)
            W_hat_m22(j,:)=x(j,12+2*N+(n^4)*5+1:12+2*N+(n^4)*5+n^4);
             y_hat_m22(j)= W_hat_m22(j,:)* xi_22';                  %钟型 2
end
y_hat_11=zeros(length(t),1);
y_hat_12=zeros(length(t),1);
y_hat_11(:,1)=y_hat_m11(:,1)+y_hat_m21(:,1);    %第一条模糊输出，关节1
y_hat_12(:,1)=y_hat_m12(:,1)+y_hat_m22(:,1);    %第一条模糊输出，关节2

% 第二条边的构成
y_hat_m31=zeros(length(t),1);
W_hat_m31=zeros(length(t),n^4);
W_hat_m32=zeros(length(t),n^4);
W_hat_m41=zeros(length(t),n^4);
W_hat_m42=zeros(length(t),n^4);
for j=1:length(t)
            W_hat_m31(j,:)=x(j,12+2*N+(n^4)*2+1:12+2*N+(n^4)*2+n^4);
             y_hat_m31(j)= W_hat_m31(j,:)* xi_32';                      %PI型 1
end
y_hat_m32=zeros(length(t),1);
for j=1:length(t)
            W_hat_m32(j,:)=x(j,12+2*N+(n^4)*6+1:12+2*N+(n^4)*6+n^4);
             y_hat_m32(j)= W_hat_m32(j,:)* xi_32';                       %PI型 2
end
y_hat_m41=zeros(length(t),1);
for j=1:length(t)
            W_hat_m41(j,:)=x(j,12+2*N+(n^4)*3+1:12+2*N+(n^4)*3+n^4);
             y_hat_m41(j)= W_hat_m41(j,:)* xi_42';                        %高斯 1
end
y_hat_m42=zeros(length(t),1);
for j=1:length(t)
            W_hat_m42(j,:)=x(j,12+2*N+(n^4)*7+1:12+2*N+(n^4)*7+n^4);
             y_hat_m42(j)= W_hat_m42(j,:)* xi_42';                             %高斯 2
end
y_hat_21=zeros(length(t),1);
y_hat_22=zeros(length(t),1);
y_hat_21(:,1)=y_hat_m31(:,1)+y_hat_m41(:,1);    %第一条模糊输出，关节1
y_hat_22(:,1)=y_hat_m32(:,1)+y_hat_m42(:,1);    %第一条模糊输出，关节2

% y_hat_m1=zeros(length(t),2);
% y_hat_m2=zeros(length(t),2);
y_hat_m1=[y_hat_11';y_hat_21']';
y_hat_m2=[y_hat_12';y_hat_22']';                        %模糊输出，关节1 关节2

%多边输出
y_hat1=zeros(length(t),1);
for j=1:length(t) 
     y_hat1(j,1)=x(j,12+1:12+N)*y_hat_m1(j,1:N)';
end
y_hat2=zeros(length(t),1);
for j=1:length(t) 
     y_hat2(j,1)=x(j,12+N+1:12+N+N)*y_hat_m2(j,1:N)';
end

%控制器设计
for i=1:length(t)
    u_tau(:,i)=ke'*e(i,:)'+kr'*xre(i,:)'-[y_hat1(i);y_hat2(i)];
end

%位置跟踪
Xref4_11=x(:,5);
Xref4_12=x(:,7);
X4_1=x(:,1);
X4_2=x(:,3);
%速度跟踪
Xref4_21=x(:,6);
Xref4_22=x(:,8);
X4_3=x(:,2);
X4_4=x(:,4);
% 跟踪误差
E4_11=e11;
E4_12=e12;
E4_21=e21;
E4_22=e22;
%输入力矩
U4=u_tau;
%非线性部分
yref4_1=y_true(1,:)';
yref4_2=y_true(2,:)';
% 非线性逼近
Y_hat4_1=y_hat1;
Y_hat4_2=y_hat2;
%逼近误差
EE4_1=y_ref1-y_hat1;
EE4_2=y_ref2-y_hat2;

disp(['状态1:',num2str(sum(abs(e11)))])
disp(['状态2:',num2str(sum(abs(e21)))])
disp(['状态3:',num2str(sum(abs(e12)))])
disp(['状态4:',num2str(sum(abs(e22)))])
disp(['状态1:',num2str(sum(abs(e11))/length(t))])
disp(['状态2:',num2str(sum(abs(e21))/length(t))])
disp(['状态3:',num2str(sum(abs(e12))/length(t))])
disp(['状态4:',num2str(sum(abs(e22))/length(t))])
disp(['状态5:',num2str(sum(abs(y_ref1-y_hat1)/length(t)))])
disp(['状态6:',num2str(sum(abs(y_ref2-y_hat2)/length(t)))])

%状态
figure()
subplot(221)
plot(t,x_ref1(:,1),'--r');hold on;
plot(t,x(:,5),'-b')
plot(t,x(:,1),'--m');grid on;
legend('r1','x-ref1','x1')
subplot(222)
plot(t,x_ref2(:,1),'--r');hold on;
plot(t,x(:,6),'-b')
plot(t,x(:,2),'--m');grid on;
legend('r2','x-ref2','x2')
subplot(223)
plot(t,x_ref1(:,2),'--r');hold on;
plot(t,x(:,7),'-b')
plot(t,x(:,3),'--m');grid on;
legend('dr1','dx-ref1','dx1')
subplot(224)
plot(t,x_ref2(:,2),'--r');hold on;
plot(t,x(:,8),'-b')
plot(t,x(:,4),'--m');grid on;
legend('dr2','dx-ref2','dx2')
%误差
figure()
subplot(421);plot(t,er11,'-b');grid on;legend('er1');
subplot(422);plot(t,er21,'-b');grid on;legend('er2');
subplot(423);plot(t,er12,'-b');grid on;legend('der1');
subplot(424);plot(t,er22,'-b');grid on;legend('der2');
subplot(425);plot(t,e11,'-b');grid on;legend('e1');
subplot(426);plot(t,e21,'-b');grid on;legend('e2');
subplot(427);plot(t,e12,'-b');grid on;legend('de1');
subplot(428);plot(t,e22,'-b');grid on;legend('de1');
%控制
figure()
plot(t,u_tau');grid on;legend('u1','u2');

figure()
plot(t,y_ref1,'-b');hold on;
plot(t,y_ref2,'-m');hold on;
legend('y_ref1','y_ref2');
% plot(t,y_ref3,'-r');hold on;
% plot(t,y_ref4,'-g');grid on;

%参数,多边学习权重
figure()
subplot(211)
% plot(t,x(:,13+N:12+2*N));grid on;
plot(t,x(:,12+1:12+N));grid on;
subplot(212)
plot(t,x(:,12+N+1:12+N+N));grid on;
legend('W1','W2');

figure()
%输出1
subplot(221)
plot(t,y_ref1,'-b');hold on;
% plot(t,y_hat_m1,'-g');hold on;
% plot(t,y_hat_m1(2),'-r');hold on;
% plot(t,y_hat_m11,'--r');hold on;
% plot(t,y_hat_m21,'--y');hold on;
% plot(t,y_hat_m31,'--b');hold on;
% plot(t,y_hat_m41,'--g');hold on;     %31、41 
plot(t,y_hat1,'-m');grid on;
% legend('y_ref1','y_hat_m1','y_hat1');
%输出2
subplot(222)
plot(t,y_ref2,'-b');hold on;
% plot(t,y_hat_m2,'-g');hold on;
% plot(t,y_hat_m12,'--r');hold on;
% plot(t,y_hat_m22,'--y');hold on;
% plot(t,y_hat_m32,'--b');hold on;
% plot(t,y_hat_m42,'--g');hold on;
plot(t,y_hat2,'-m');grid on;
% legend('y_ref2','y_hat_m2','y_hat2');
% 逼近误差1
subplot(223)
% plot(t,repmat(y_ref1,1,N)-y_hat_m1);hold on;
plot(t,y_ref1-y_hat1,'--m');grid on;
% legend('E1','E_y{1}');
% 逼近误差2
subplot(224)
% plot(t,repmat(y_ref2,1,N)-y_hat_m2);hold on;
plot(t,y_ref2-y_hat2,'--m');grid on;
% legend('E2','E_y{2}');