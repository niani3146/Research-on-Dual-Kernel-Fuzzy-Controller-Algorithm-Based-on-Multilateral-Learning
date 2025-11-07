%PID控制
% clear all;
% close all;
t0=0;tf=30;
x0=[1,1,1,1,0];
options = odeset('MaxStep', 0.001); 
[t,x]=ode45(@PID_test_ode,[t0,tf],x0,options);

A=[0,1;0,0];%倒立摆系统
B=[0;1];%倒立摆系统
Ar=[0,1;-1,-2];%参考系统
Br=[0;1];%参考系统

xref=[sin(t),cos(t)];

e=x(:,3)-x(:,1);
de=x(:,4)-x(:,2);

Kp=[50,30];
Ki=50;
u=Kp*[e,de]'+Ki*x(:,5)';

X11=x(:,1);
X12=x(:,2);
E11=e;
E12=de;
U1=u;

disp(['状态1：',num2str(sum(abs(e))),'；'])
disp(['状态2：',num2str(sum(abs(de))),'；'])
disp(['状态1：',num2str(sum(abs(e))/length(t)),'；'])
disp(['状态2：',num2str(sum(abs(de))/length(t)),'；'])

% figure()
% subplot(2,1,1)
% plot(t,x(:,3),'-b');hold on;grid on;
% plot(t,x(:,1),'--m');legend('x1-l','x1-nl')
% subplot(2,1,2)
% plot(t,x(:,4),'-b');hold on;grid on;
% plot(t,x(:,2),'--m');legend('x2-l','x2-nl')
% figure()
% subplot(211);plot(t,e);grid on;legend('e');
% subplot(212);plot(t,de);grid on;legend('de');
% figure('name','control signal')
% plot(t,u');grid on;legend('u');