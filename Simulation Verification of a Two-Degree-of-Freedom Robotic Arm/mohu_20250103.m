 % %二自由度机械臂-模糊多边-模糊规则1，N=7
% clear all;
% close all;
N=7;%选取的边数
n=3;%隶属度函数 数量

t0=0;tf=30;
T=linspace(t0,tf,100);

para_a=linspace(-0.1,0.2,N);
para_b=linspace(-0.1,0.2,N);
% for i = 1:N
% %     min_val = 0.08 + (i - 1) * 0.015; % 最小值
% %     max_val = 0.08 + (i - 1) * 0.03; % 最大值
%     min_val = 0.0001 + (i - 1) * 0.0001; % 最小值
%     max_val = 0.0001 + (i - 1) * 0.0013; % 最大值
%     para_k((i-1)*n^4+1:(i-1)*n^4+n^4) = linspace(min_val, max_val, n^4)'; % 每列均匀分布
% end
para_k=zeros(1,N*(n^4));
% for i = 1:N
%     min_val = 0.008 + (i - 1) * 0.0008; % 最小值
%     max_val = 0.008 + (i - 1) * 0.0013; % 最大值
%     para_w((i-1)*n^4+1:(i-1)*n^4+n^4) = linspace(min_val, max_val, n^4)'; % 每列均匀分布
% end
% para_w=linspace(-0.01,0,N*n^4);
% para_k=ones(1,N*n^4);
para_w=zeros(1,N*(n^4));
para_n=zeros(1,n^4);
para_m=zeros(1,n^4);
para_c=zeros(1,2*N);

x0=[1,1,1,1,1,1,1,1,0,0,0,0,... 
  para_a,para_b,para_k,para_w,para_n,para_m,para_c];

options = odeset('MaxStep', 0.01); % 设置最大步长
% options = odeset('MaxStep', 0.02); % 设置最大步长
[t,x]=ode45(@mohu_fn_20250103,[t0,tf],x0,options);

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

for j=1:length(t)
        y_ref3(j,:)=x(j,12+2*N+(n^4)*N+(n^4)*N+2*(n^4)+1);
        y_ref4(j,:)=x(j,12+2*N+(n^4)*N+(n^4)*N+2*(n^4)+2);
end

u_tau=zeros(2,length(t));

%模糊控制
xi_1= x(12+2*N+2*(n^4)*N+1:12+2*N+2*(n^4)*N+n^4);

y_hat_m1=zeros(length(t),N);
for j=1:length(t)
       for i=1:N
            W_hat_m1=x(j,12+2*N+(i-1)*(n^4)+1:12+2*N+(i-1)*(n^4)+n^4);
             y_hat_m1(j,i)= W_hat_m1* xi_1';    
       end
end
y_hat_m2=zeros(length(t),N);
for j=1:length(t)
       for i=1:N
            W_hat_m2=x(j,12+2*N+(n^4)*N+(i-1)*(n^4)+1:12+2*N+(n^4)*N+(i-1)*(n^4)+n^4);
             y_hat_m2(j,i)= W_hat_m2* xi_1';    
       end
end
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

disp(['状态1:',num2str(sum(abs(e11)))])
disp(['状态2:',num2str(sum(abs(e21)))])
disp(['状态3:',num2str(sum(abs(e12)))])
disp(['状态4:',num2str(sum(abs(e22)))])
disp(['状态1:',num2str(sum(abs(e11))/length(t))])
disp(['状态2:',num2str(sum(abs(e21))/length(t))])
disp(['状态3:',num2str(sum(abs(e12))/length(t))])
disp(['状态4:',num2str(sum(abs(e22))/length(t))])
disp(['状态5:',num2str(sum(abs(y_ref1-y_hat1)))])
disp(['状态6:',num2str(sum(abs(y_ref2-y_hat2)))])

%位置跟踪
Xref3_1=x(:,5);
X5_1=x(:,1);
X5_2=x(:,3);
%速度跟踪
Xref3_2=x(:,7);
X5_3=x(:,2);
X5_4=x(:,4);

% 跟踪误差
E5_11=e11;
E5_12=e12;
E5_21=e21;
E5_22=e22;
%输入力矩
U5=u_tau;
%非线性参考
yref5_1=y_true(1,:)';
yref5_2=y_true(2,:)';
Y_ref=[yref5_1,yref5_2];
% 非线性逼近
Y_hat5_1=y_hat1;
Y_hat5_2=y_hat2;
%逼近误差
EE5_1=y_ref1-y_hat1;
EE5_2=y_ref2-y_hat2;
% 
% %状态
% figure()
% subplot(221)
% plot(t,x_ref1(:,1),'--r');hold on;
% plot(t,x(:,5),'-b')
% plot(t,x(:,1),'--m');grid on;
% legend('r1','x-ref1','x1')
% subplot(222)
% plot(t,x_ref2(:,1),'--r');hold on;
% plot(t,x(:,6),'-b')
% plot(t,x(:,2),'--m');grid on;
% legend('r2','x-ref2','x2')
% subplot(223)
% plot(t,x_ref1(:,2),'--r');hold on;
% plot(t,x(:,7),'-b')
% plot(t,x(:,3),'--m');grid on;
% legend('dr1','dx-ref1','dx1')
% subplot(224)
% plot(t,x_ref2(:,2),'--r');hold on;
% plot(t,x(:,8),'-b')
% plot(t,x(:,4),'--m');grid on;
% legend('dr2','dx-ref2','dx2')
% %误差
% figure()
% subplot(421);plot(t,er11,'-b');grid on;legend('er1');
% subplot(422);plot(t,er21,'-b');grid on;legend('er2');
% subplot(423);plot(t,er12,'-b');grid on;legend('der1');
% subplot(424);plot(t,er22,'-b');grid on;legend('der2');
% subplot(425);plot(t,e11,'-b');grid on;legend('e1');
% subplot(426);plot(t,e21,'-b');grid on;legend('e2');
% subplot(427);plot(t,e12,'-b');grid on;legend('de1');
% subplot(428);plot(t,e22,'-b');grid on;legend('de1');
% %控制
% figure()
% plot(t,u_tau');grid on;legend('u1','u2');
% 
% figure()
% plot(t,y_ref1,'-b');hold on;
% plot(t,y_ref2,'-m');hold on;
% legend('y_ref1','y_ref2');
% % plot(t,y_ref3,'-r');hold on;
% % plot(t,y_ref4,'-g');grid on;
% 
% %参数,多边学习权重
% figure()
% subplot(211)
% % plot(t,x(:,13+N:12+2*N));grid on;
% plot(t,x(:,12+1:12+N));grid on;
% subplot(212)
% plot(t,x(:,12+N+1:12+N+N));grid on;
% legend('W1','W2');
% % 
% figure()
% %输出1
% subplot(221)
% plot(t,y_ref1,'-b');hold on;
% plot(t,y_hat_m1,'-g');hold on;
% plot(t,y_hat1,'--m');grid on;
% % legend('y_ref1','y_hat_m1','y_hat1');
% %输出2
% subplot(222)
% plot(t,y_ref2,'-b');hold on;
% plot(t,y_hat_m2,'-g');hold on;
% plot(t,y_hat2,'--m');grid on;
% % legend('y_ref2','y_hat_m2','y_hat2');
% % 逼近误差1
% subplot(223)
% plot(t,repmat(y_ref1,1,N)-y_hat_m1);hold on;
% plot(t,y_ref1-y_hat1,'--m');grid on;
% % legend('E1','E_y{1}');
% % 逼近误差2
% subplot(224)
% plot(t,repmat(y_ref2,1,N)-y_hat_m2);hold on;
% plot(t,y_ref2-y_hat2,'--m');grid on;
% % legend('E2','E_y{2}');