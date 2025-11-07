%隶属函数为高斯函数的模糊自适应控制,单边模糊输出
% clear all;
% close all;

options = odeset('MaxStep', 0.001); % 设置最大步长
n=7;
% dx=zeros(6+n+25*n+25+2*n,1);
% dx=zeros(6+25*n+25+2*n,1);
t0=0;
tf=30;
T=linspace(t0,tf,100);%创建了一个包含100个元素的线性空间，这些元素均匀地分布在t0和tf之间
% 
% for i = 1:n
%     min_val = 0.001 + (i - 1) * 0.0001; % 最小值
%     max_val = 0.001 + (i - 1) * 0.0003; % 最大值
%     para_k((i-1)*25+1:(i-1)*25+25) = linspace(min_val, max_val, 25)'; % 每列均匀分布
% end
para_k=zeros(1,25*n);
para_m=zeros(1,25);
para_d=zeros(1,2*n);

% x0=[1,1,1,1,0,0,...
%   para_b,para_k,para_m,para_d];
x0=[1,1,1,1,0,0,...
para_k,para_m,para_d];

[t,x]=ode45(@mohu_fn_test,[t0,tf],x0,options);%求解仿射非线性方程的状态变量
 
A=[0,1;0,0];%倒立摆系统
B=[0;1];%倒立摆系统
Ar=[0,1;-1,-2];%参考系统
Br=[0;1];%参考系统

xi1 = [x(1);x(2)]; % 位置、速度（角度、角速度）信息
xi2 = [x(3);x(4)];
xi = [xi1;xi2]; 
x_ref=[sin(t),cos(t)];

xre=[x(:,3),x(:,4),x_ref(:,1)];

er1=x_ref(:,1)-x(:,3);
er2=x_ref(:,2)-x(:,4);
e1=x(:,3)-x(:,1);   
e2=x(:,4)-x(:,2);   
ee=[e1,e2];

ke=[50;30];
kr=[-1;-2;1];

%真实动态
Wd_vec =[1,-1,0.5]';

y_ref=zeros(length(t),1);
fphi=zeros(length(t),3);
for i=1:length(t)
    fphi(i,:)=[exp(x(i,1)*x(i,2)),sin(x(i,1)),abs(x(i,2))*x(i,2)];
    y_ref(i)=Wd_vec'*fphi(i,:)';
end

%模糊控制 
xi_1= x(189:213);
y_hat_m=zeros(length(t),1);
for j=1:length(t)
        W_hat_m=x(j,6+n+1:6+n+25);
         y_hat_m(j,1)= W_hat_m* xi_1';    
end 

ee_yn=y_ref-y_hat_m; 

u_tau=zeros(length(t),1);
for i=1:length(t)
    u_tau(i)=ke'*ee(i,:)'+kr'*xre(i,:)'-y_hat_m(i)';
end

Xref1=x(:,3);
Xref2=x(:,4);
X21=x(:,1);
X22=x(:,2);
E21=e1;
E22=e2;
U2=u_tau;
EE_1=ee_yn;
Yref=y_ref;
Y_hat_2=y_hat_m;


% disp(['状态1:',num2str(sum(abs(e1)))])
% disp(['状态2:',num2str(sum(abs(e2)))])
% disp(['状态1:',num2str(sum(abs(e1))/length(t))])
% disp(['状态2:',num2str(sum(abs(e2))/length(t))])

%状态
% figure()
% subplot(211)
% plot(t,x_ref(:,1),'--r');hold on;
% plot(t,x(:,3),'-b')
% plot(t,x(:,1),'--m');grid on;
% legend('r1','x-ref1','x1')
% subplot(212)
% plot(t,x_ref(:,2),'--r');hold on;
% plot(t,x(:,4),'-b')
% plot(t,x(:,2),'--m');grid on;
% legend('r2','x-ref2','x2')
% %误差
% figure()
% subplot(211)
% plot(t,er1,'-b');hold on;
% plot(t,er2,'-m');grid on;legend('er1','er2')
% subplot(212)
% plot(t,e1,'-b');hold on;
% plot(t,e2,'-m');grid on;legend('e1','e2')
% %控制
% figure()
% plot(t,u_tau);grid on;legend('u')
% %更新率
% % figure()
% % plot(t,x(:,7:6+n),'-b');
% % hold on;grid on;legend('k1','k2')
% %输出
% figure()
% subplot(211)
% plot(t,y_ref,'--b');hold on;
% plot(t,y_hat_m);hold on;
% % plot(t,y_hat,'--m');
% grid on;
% subplot(212)
% plot(t,repmat(y_ref,1,n)-y_hat_m);hold on;
% grid on;
% figure;
% plot(t, y_hat_m, 'g', 'LineWidth', 2);hold on;
% plot(t, y_ref, 'b', 'LineWidth', 2);hold on;
% legend('y_{hat_m}','y_{ref}')
% xlabel('Time (s)');
% ylabel('y_{ref}');
% title('Reference Output y_{ref}');
% grid on;