n=2;%选取的边数
ee_yn_pre=zeros(n,1);

t0=0;tf=30;
T=linspace(t0,tf,100);

para_y=linspace(-1,1,n);
para_k=ones(1,n)/n;

para_d=zeros(1,2*n+n);
x0=[1,1,1,1,0,0,...
    para_y,para_k,para_d];
options = odeset('MaxStep', 0.001); % 设置最大步长
[t,x]=ode45(@fn_MRAC_ode,[t0,tf],x0,options);

x_ref=[sin(t),cos(t)];

Wd_vec=[1,-1,0.5]';

kr=[-1 -2 1]';
xre=[x(:,3),x(:,4),x_ref(:,1)];
ke=[50;30];
er1=x_ref(:,1)-x(:,3);
er2=x_ref(:,2)-x(:,4);
e1=x(:,3)-x(:,1);
e2=x(:,4)-x(:,2);
ee=[e1,e2];

y_ref=zeros(length(t),1);
fphi=zeros(length(t),3);
for i=1:length(t)
    fphi(i,:)=[exp(x(i,1)*x(i,2)),sin(x(i,1)),abs(x(i,2))*x(i,2)];
    y_ref(i)=Wd_vec'*fphi(i,:)';
end

y_hat=zeros(length(t),1);
u_tau=zeros(length(t),1);
for i=1:length(t)
    y_hat(i)=x(i,7+n:6+2*n)*x(i,7:6+n)';
    u_tau(i)=ke'*ee(i,:)'+kr'*xre(i,:)'-y_hat(i);
end
e_y=y_ref-y_hat;  

disp(['state 1:',num2str(sum(abs(e1)))])
disp(['state 2:',num2str(sum(abs(e2)))])
disp(['state 1:',num2str(sum(abs(e1))/length(t))])
disp(['state 2:',num2str(sum(abs(e2))/length(t))])

X31=x(:,1);
X32=x(:,2);
E31=e1;
E32=e2;
U3=u_tau;
EE_3=e_y;
Y_hat_3=y_hat;

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

% figure()
% subplot(211)
% plot(t,er1,'-b');hold on;
% plot(t,er2,'-m');grid on;legend('er1','er1')
% subplot(212)
% plot(t,e1,'-b');hold on;
% plot(t,e2,'-m');grid on;legend('e1','e2')

% figure()
% plot(t,u_tau);grid on;legend('u')
% 
% figure()
% plot(t,x(:,7+n:6+2*n));

% hold on;grid on;legend('k1','k2','k3','k4','k5','k6','k7')
