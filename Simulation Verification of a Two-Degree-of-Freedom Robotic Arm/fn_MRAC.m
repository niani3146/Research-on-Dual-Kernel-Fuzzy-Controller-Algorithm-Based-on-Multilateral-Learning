n=2;

t0=0;tf=30;
T=linspace(t0,tf,100);
para_y=linspace(-1,1,n);
para_k=ones(1,n)/n;
para_i=zeros(1,2*n);
para_m=zeros(1,2);
x0=[1,1,1,1,1,1,1,1,0,0,0,0,...
    repmat([para_y,para_k],1,2),para_i,para_m];
options = odeset('MaxStep', 0.02);
[t,x]=ode45(@fn_MRAC_ode,[t0,tf],x0,options);

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

a0=1;
a1=1;
m1=0.1;
m2=0.1;
gs=9.8;
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
xre=[x(:,5),x(:,6),x(:,7),x(:,8),x_ref1(:,1),x_ref2(:,1)];

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
        y_ref3(j,:)=x(j,2*6+4*n+2*n+1);
        y_ref4(j,:)=x(j,2*6+4*n+2*n+2);
end

y_hat1=zeros(length(t),1);
y_hat2=zeros(length(t),1);
u_tau=zeros(2,length(t));
for i=1:length(t)
    y_hat1(i)=x(i,13+n:12+2*n)*x(i,13:12+n)';
    y_hat2(i)=x(i,13+3*n:12+4*n)*x(i,13+2*n:12+3*n)';
    u_tau(:,i)=ke'*e(i,:)'+kr'*xre(i,:)'-[y_hat1(i);y_hat2(i)];
end

Xref3_11=x(:,5);
Xref3_12=x(:,7);
X3_1=x(:,1);
X3_2=x(:,3);

Xref3_21=x(:,6);
Xref3_22=x(:,8);
X3_3=x(:,2);
X3_4=x(:,4);

E3_11=e11;
E3_12=e12;
E3_21=e21;
E3_22=e22;

U3=u_tau;

yref3_1=y_true(1,:)';
yref3_2=y_true(2,:)';
Y_ref=[yref3_1,yref3_2];

Y_hat3_1=y_hat1;
Y_hat3_2=y_hat2;

EE3_1=y_ref1-y_hat1;
EE3_2=y_ref2-y_hat2;

disp(['state 1:',num2str(sum(abs(e11)))])
disp(['state 2:',num2str(sum(abs(e21)))])
disp(['state 3:',num2str(sum(abs(e12)))])
disp(['state 4:',num2str(sum(abs(e22)))])
disp(['state 1:',num2str(sum(abs(e11))/length(t))])
disp(['state 2:',num2str(sum(abs(e21))/length(t))])
disp(['state 3:',num2str(sum(abs(e12))/length(t))])
disp(['state 4:',num2str(sum(abs(e22))/length(t))])

% % figure()
% % subplot(221)
% % plot(t,x_ref1(:,1),'--r');hold on;
% % plot(t,x(:,5),'-b')
% % plot(t,x(:,1),'--m');grid on;
% % legend('r1','x-ref1','x1')
% % subplot(222)
% % plot(t,x_ref2(:,1),'--r');hold on;
% % plot(t,x(:,6),'-b')
% % plot(t,x(:,2),'--m');grid on;
% % legend('r2','x-ref2','x2')
% % subplot(223)
% % plot(t,x_ref1(:,2),'--r');hold on;
% % plot(t,x(:,7),'-b')
% % plot(t,x(:,3),'--m');grid on;
% % legend('dr1','dx-ref1','dx1')
% % subplot(224)
% % plot(t,x_ref2(:,2),'--r');hold on;
% % plot(t,x(:,8),'-b')
% % plot(t,x(:,4),'--m');grid on;
% % legend('dr2','dx-ref2','dx2')

% % figure()
% % subplot(421);plot(t,er11,'-b');grid on;legend('er1');
% % subplot(422);plot(t,er21,'-b');grid on;legend('er2');
% % subplot(423);plot(t,er12,'-b');grid on;legend('der1');
% % subplot(424);plot(t,er22,'-b');grid on;legend('der2');
% % subplot(425);plot(t,e11,'-b');grid on;legend('e1');
% % subplot(426);plot(t,e21,'-b');grid on;legend('e2');
% % subplot(427);plot(t,e12,'-b');grid on;legend('de1');
% % subplot(428);plot(t,e22,'-b');grid on;legend('de1');

% % figure()
% % plot(t,u_tau');grid on;legend('u1','u2');
% % 
% % figure()
% % plot(t,y_ref1,'-b');hold on;
% % plot(t,y_ref2,'-m');hold on;
% % plot(t,y_ref3,'-r');hold on;
% % plot(t,y_ref4,'-g');grid on;

% % figure()
% % subplot(211)
% % plot(t,x(:,13+n:12+2*n));grid on;
% % subplot(212)
% % plot(t,x(:,13+3*n:12+4*n));grid on;
% % figure()
% % subplot(221)
% % plot(t,y_ref1,'-b');hold on;
% % plot(t,x(:,13:12+n));
% % plot(t,y_hat1,'--m');grid on;
% % subplot(222)
% % plot(t,y_ref2,'-b');hold on;
% % plot(t,x(:,13+2*n:12+3*n));
% % plot(t,y_hat2,'--m');grid on;
% % subplot(223)
% % plot(t,repmat(y_ref1,1,n)-x(:,13:12+n));hold on;
% % plot(t,y_ref1-y_hat1,'--m');grid on;
% % subplot(224)
% % plot(t,repmat(y_ref2,1,n)-x(:,13+2*n:12+3*n));hold on;

% % plot(t,y_ref2-y_hat2,'--m');grid on;
