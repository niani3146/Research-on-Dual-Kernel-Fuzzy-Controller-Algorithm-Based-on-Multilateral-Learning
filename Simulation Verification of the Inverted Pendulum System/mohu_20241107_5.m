% %隶属函数为高斯函数的模糊自适应控制，多边模糊输出
% clear all;
% close all;
options = odeset('MaxStep', 0.001); % 设置最大步长
n=7;
N=24000;

t0=0;
tf=30;
T=linspace(t0,tf,100);
% para_b=ones(1,n)/(3*n);
para_b=linspace(-0.5,0.5,n);
% para_b=linspace(-0.02,0.02,n);

for i = 1:n
    min_val = 0.0001 + (i - 1) * 0.0001; % 最小值
    max_val = 0.0001 + (i - 1) * 0.0003; % 最大值
    para_k((i-1)*25+1:(i-1)*25+25) = linspace(min_val, max_val, 25)'; % 每列均匀分布
end
% para_k=zeros(1,25*n);
para_m=zeros(1,25);
para_d=zeros(1,2*n);

x0=[1,1,1,1,0,0,...
  para_b,para_k,para_m,para_d];

[t,x]=ode45(@mohu_fn_20241107_5,[t0,tf],x0,options);%求解仿射非线性方程的状态变量
% [t,x] =equidistant_sampling(t1, x1, N);
% % [t,x] =mohu_fn_20241113_5(t1, x1);

A=[0,1;0,0];
B=[0;1];
Ar=[0,1;-1,-2];
Br=[0;1];

xi1 = [x(1);x(2)];
xi2 = [x(3);x(4)];
xi = [xi1;xi2];
x_ref=[sin(t),cos(t)];

xre=[x(:,3),x(:,4),x_ref(:,1)];

er1=x_ref(:,1)-x(:,3);
er2=x_ref(:,2)-x(:,4);
e1=x(:,3)-x(:,1);   
e2=x(:,4)-x(:,2);   
ee=[e1,e2];

% ke=[50;30];
ke=[50;30];
kr=[-1.5;-3;1.5];

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
y_hat_m=zeros(length(t),7);
for j=1:length(t)
       for i=1:7
            W_hat_m=x(j,6+n+(i-1)*25+1:6+n+(i-1)*25+25);
             y_hat_m(j,i)= W_hat_m* xi_1';    
       end
end

y_hat=zeros(length(t),1);
for j=1:length(t) 
     y_hat(j,1)=x(j,7:6+n)*y_hat_m(j,:)';
end

ee_yn1=y_ref-y_hat; 

u_tau1=zeros(length(t),1);
for i=1:length(t)
    u_tau1(i)=ke'*ee(i,:)'+kr'*xre(i,:)'-y_hat(i)';
end


X51=x(:,1);
X52=x(:,2);
E51=e1;
E52=e2;
U5=u_tau1;
EE_5=ee_yn1;
Y_hat_5=y_hat;

disp(['状态1:',num2str(sum(abs(e1)))])
disp(['状态2:',num2str(sum(abs(e2)))])
disp(['状态1:',num2str(sum(abs(e1))/length(t))])
disp(['状态2:',num2str(sum(abs(e2))/length(t))])

% %状态
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
% plot(t,u_tau1);grid on;legend('u')
% %更新率
% figure()
% plot(t,x(:,7:6+n),'-b');
% hold on;grid on;legend('k1','k2')
% %输出
% figure()
% subplot(211)
% plot(t,y_ref,'--b');hold on;
% % plot(t,y_hat_m);hold on;
% plot(t,y_hat,'--m');
% grid on;
% subplot(212)
% % plot(t,repmat(y_ref,1,n)-y_hat_m);hold on;
% plot(t,y_ref-y_hat,'--m');
% grid on;
% figure;
% plot(t, y_hat, 'g', 'LineWidth', 0.5);hold on;
% plot(t, y_ref, 'b', 'LineWidth', 0.5);hold on;
% plot(t, y_hat_m, 'r', 'LineWidth', 0.5);hold on;
% legend('y_{hat}','y_{ref}','y_{hat_m}')
% xlabel('Time (s)');
% ylabel('y_{ref}');
% title('Reference Output y_{ref}');
% grid on;