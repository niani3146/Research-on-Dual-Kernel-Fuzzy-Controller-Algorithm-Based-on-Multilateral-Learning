%PID
% close all;
% clear all;
t0=0;tf=30;
T=linspace(t0,tf,1000);
x0=[1,1,1,1,....
    1,1,1,1,...
    0,0];
options=odeset('MaxStep',0.01);
[t,x]=ode45(@PID_RRbot_ode,[t0,tf],x0,options);

A=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];%RR-bot系统
B=[0 0;0 0;1 0;0 1];%RR-bot系统
Ar=[0,0,1,0;0,0,0,1;-1,0,-2,0;0,-1,0,-2];%参考系统
Br=[0 0;0 0;1 0;0 1];%参考系统

%RR-+-

% -+++++++++++++++++++++++++++++++++++
                                                                                                                             
a0=1;%连杆1长度为1
a1=1;%连杆2长度为1
m1=0.1;%连杆1的质量为0.1kg
m2=0.1;%连杆2的质量为0.2kg
gs=9.8;%重力加速度

Wd=[2,0,1,0;...
    0,3,0,1]';
fphi=zeros(2,length(t));
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
    fphi(:,i)=[un_model(1);un_model(2)];
end

r1=sin(t);r2=cos(t);
% r1=zeros(size(t));%列向量
% r2=zeros(size(t));
% for i=1:length(t)
%     if t(i)>=20 && t(i)<50
%         r1(i)=1;r2(i)=2;
% %     elseif t(i)>50 && t(i)<65
% %         r1(i)=3;r2(i)=1.5;
%     else
%         r1(i)=0;r2(i)=0;
%     end
% end

Kp=[500,500];
Kd=[100,100];
Ki=[200,200];

e11=x(:,5)-x(:,1);
e12=x(:,7)-x(:,3);
e21=x(:,6)-x(:,2);
e22=x(:,8)-x(:,4);
e=[e11,e21,e12,e22];

U=zeros(2,length(t));
for i=1:size(t,1)
    U(1,i)=Kp(1)*e11(i)'+Kd(1)*e12(i)+Ki(1)*x(i,9);
    U(2,i)=Kp(2)*e21(i)'+Kd(2)*e22(i)+Ki(2)*x(i,10);
end

%位置跟踪
Xref1_11=x(:,5);
Xref1_12=x(:,7);
X1_1=x(:,1);
X1_2=x(:,3);
%速度跟踪
Xref1_21=x(:,6);
Xref1_22=x(:,8);
X1_3=x(:,2);
X1_4=x(:,4);
% 跟踪误差
E1_11=e11;
E1_12=e12;
E1_21=e21;
E1_22=e22;
%输入力矩
U1=U;

% %非线性部分
% yref2_1=y_true(1,:)';
% yref2_2=y_true(2,:)';
% %逼近误差
% EE2_1=y_ref1-y_hat1;
% EE2_2=y_ref2-y_hat2;

% 
% disp(['状态1:',num2str(sum(abs(e(:,1)))),';'])
% disp(['状态2:',num2str(sum(abs(e(:,2)))),';'])
% disp(['状态3:',num2str(sum(abs(e(:,3)))),';'])
% disp(['状态4:',num2str(sum(abs(e(:,4)))),'.'])
% disp(['状态1:',num2str(sum(abs(e(:,1)))/length(t))])
% disp(['状态2:',num2str(sum(abs(e(:,2)))/length(t))])
% disp(['状态3:',num2str(sum(abs(e(:,3)))/length(t))])
% disp(['状态4:',num2str(sum(abs(e(:,4)))/length(t))])
% 
% 
% 
% %状态
% figure()
% subplot(221)
% plot(t,x(:,1),'-r')
% grid on;hold on;
% plot(t,x(:,5),'--b')
% legend('x1-dyn','x1-ref')
% xlabel('time(t/s)')
% ylabel('state1 value')
% subplot(222)
% plot(t,x(:,2),'-r')
% grid on;hold on;
% plot(t,x(:,6),'--b')
% legend('x2-dyn','x2-ref')
% xlabel('time(t/s)')
% ylabel('state2 value')
% subplot(223)
% plot(t,x(:,3),'-r')
% grid on;hold on;
% plot(t,x(:,7),'--b')
% legend('x3-dyn','x3-ref')
% xlabel('time(t/s)')
% ylabel('state3 value')
% subplot(224)
% plot(t,x(:,4),'-r')
% grid on;hold on;
% plot(t,x(:,8),'--b')
% legend('x4-dyn','x4-ref')
% xlabel('time(t/s)')
% ylabel('state4 value')
% %误差
% figure()
% subplot(221);plot(t,e11);grid on;legend('e11');
% subplot(222);plot(t,e21);grid on;legend('e21');
% subplot(223);plot(t,e12);grid on;legend('e12');
% subplot(224);plot(t,e22);grid on;legend('e22');
% %控制
% figure()
% plot(t,U');grid on
% legend('U1','U2')
% xlabel('time(t/s)')
% ylabel('control signal')