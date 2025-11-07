N=2;
n=3;
m=2;

t0=0;tf=30;
T=linspace(t0,tf,100);

para_a=linspace(0.25,0.25,N);
para_k=zeros(1,N*m*(n^2));
para_w=zeros(1,N*m*(n^2));
para_c=zeros(1,1);

x0=[1,1,1,1,0,0,... 
  para_a,para_k,para_w,para_c];

options = odeset('MaxStep', 0.001); 
[t,x]=ode45(@mohu_2_fn_20250304,[t0,tf],x0,options);

x_ref1=[sin(t),cos(t)];

xi1 = [x(1);x(2)];
xi2 = [x(3);x(4)];
xi = [xi1;xi2];
x_ref=[sin(t),cos(t)];

xre=[x(:,3),x(:,4),x_ref(:,1)];

er1=x_ref(:,1)-x(:,3);
er2=x_ref(:,2)-x(:,4);
e1=x(:,3)-x(:,1);   
e2=x(:,4)-x(:,2);   
e=[e1,e2];

% ke=[50;30];
ke=[50;30];
kr=[-1.5;-3;1.5];

Wd_vec =[1,-1,0.5]';

y_ref=zeros(length(t),1);
fphi=zeros(length(t),3);
for i=1:length(t)
    fphi(i,:)=[exp(x(i,1)*x(i,2)),sin(x(i,1)),abs(x(i,2))*x(i,2)];
    y_ref(i)=Wd_vec'*fphi(i,:)';
end

u_tau=zeros(2,length(t));

xi_12= x(6+N+(n^2)*4+1:6+N+(n^2)*4+n^2);
xi_22=x(6+N+(n^2)*5+1:6+N+(n^2)*5+n^2);
xi_32=x(6+N+(n^2)*6+1:6+N+(n^2)*6+n^2);
xi_42=x(6+N+(n^2)*7+1:6+N+(n^2)*7+n^2);

y_hat_m11=zeros(length(t),1);
y_hat_m21=zeros(length(t),1);
W_hat_m11=zeros(length(t),n^2);
W_hat_m21=zeros(length(t),n^2);

for j=1:length(t)
            W_hat_m11(j,:)=x(j,6+N+1:6+N+n^2);
             y_hat_m11(j)= W_hat_m11(j,:)* xi_12';            
end

for j=1:length(t)
            W_hat_m21(j,:)=x(j,6+N+(n^2)+1:6+N+(n^2)+n^2);
             y_hat_m21(j)= W_hat_m21(j,:)* xi_22';             
end

y_hat_11=zeros(length(t),1);

y_hat_11(:,1)=y_hat_m11(:,1)+y_hat_m21(:,1);  

y_hat_m31=zeros(length(t),1);
y_hat_m41=zeros(length(t),1);
W_hat_m31=zeros(length(t),n^2);
W_hat_m41=zeros(length(t),n^2);

for j=1:length(t)
            W_hat_m31(j,:)=x(j,6+N+(n^2)*2+1:6+N+(n^2)*2+n^2);
             y_hat_m31(j)= W_hat_m31(j,:)* xi_32';                
end

for j=1:length(t)
            W_hat_m41(j,:)=x(j,6+N+(n^2)*3+1:6+N+(n^2)*3+n^2);
             y_hat_m41(j)= W_hat_m41(j,:)* xi_42';              
end
y_hat_21=zeros(length(t),1);
y_hat_21(:,1)=y_hat_m31(:,1)+y_hat_m41(:,1);  

y_hat_m1=[y_hat_11,y_hat_21];

y_hat1=zeros(length(t),1);
for j=1:length(t) 
     y_hat1(j,1)=x(j,6+1:6+N)*y_hat_m1(j,:)';
end

for i=1:length(t)
    u_tau(:,i)=ke'*e(i,:)'+kr'*xre(i,:)'-y_hat1(i);
end

Xref4_1=x(:,3);
X41=x(:,1);
X42=x(:,2);

Xref4_2=x(:,4);
X43=x(:,3);
X44=x(:,4);

E41=e1;
E42=e2;

U4=u_tau;

yref_4=y_ref;

Y_hat_4=y_hat1;

EE4_1=y_ref-y_hat1;


% disp(['state 1:',num2str(sum(abs(e1)))])
% disp(['state 2:',num2str(sum(abs(e2)))])
% disp(['state 1:',num2str(sum(abs(e1))/length(t))])
% disp(['state 2:',num2str(sum(abs(e2))/length(t))])
% disp(['state 5:',num2str(sum(abs(y_ref-y_hat1)/length(t)))])
% 
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
% plot(t,er2,'-m');grid on;legend('er1','er2')
% subplot(212)
% plot(t,e1,'-b');hold on;
% plot(t,e2,'-m');grid on;legend('e1','e2')

% figure()
% plot(t,u_tau);grid on;legend('u')

% figure()
% plot(t,x(:,7:6+N),'-b');
% hold on;grid on;legend('k1','k2')

% figure()
% subplot(211)
% plot(t,y_ref,'--b');hold on;
% % plot(t,y_hat_m);hold on;
% plot(t,y_hat1,'--m');
% grid on;
% subplot(212)
% % plot(t,repmat(y_ref,1,n)-y_hat_m);hold on;
% plot(t,y_ref-y_hat1,'--m');
% grid on;
% figure;
% plot(t, y_hat1, 'r', 'LineWidth', 0.5);hold on;
% plot(t, y_ref, 'b', 'LineWidth', 0.5);hold on;
% plot(t, y_hat_m1, 'g', 'LineWidth', 0.5);hold on;
% legend('y_{hat}','y_{ref}','y_{hat_m}')
% xlabel('Time (s)');
% ylabel('y_{ref}');
% title('Reference Output y_{ref}');
% grid on;
% 
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
% % legend('y_ref1','y_ref2');
% % % plot(t,y_ref3,'-r');hold on;
% % % plot(t,y_ref4,'-g');grid on;

% % figure()
% % subplot(211)
% % % plot(t,x(:,13+N:12+2*N));grid on;
% % plot(t,x(:,12+1:12+N));grid on;
% % subplot(212)
% % plot(t,x(:,12+N+1:12+N+N));grid on;
% % legend('W1','W2');
% % 
% % figure()
% % subplot(221)
% % plot(t,y_ref1,'-b');hold on;
% % % plot(t,y_hat_m1,'-g');hold on;
% % % plot(t,y_hat_m1(2),'-r');hold on;
% % % plot(t,y_hat_m11,'--r');hold on;
% % % plot(t,y_hat_m21,'--y');hold on;
% % % plot(t,y_hat_m31,'--b');hold on;
% % % plot(t,y_hat_m41,'--g');hold on;     
% % plot(t,y_hat1,'-m');grid on;
% % % legend('y_ref1','y_hat_m1','y_hat1');

% % subplot(222)
% % plot(t,y_ref2,'-b');hold on;
% % % plot(t,y_hat_m2,'-g');hold on;
% % % plot(t,y_hat_m12,'--r');hold on;
% % % plot(t,y_hat_m22,'--y');hold on;
% % % plot(t,y_hat_m32,'--b');hold on;
% % % plot(t,y_hat_m42,'--g');hold on;
% % plot(t,y_hat2,'-m');grid on;
% % % legend('y_ref2','y_hat_m2','y_hat2');

% % subplot(223)
% % % plot(t,repmat(y_ref1,1,N)-y_hat_m1);hold on;
% % plot(t,y_ref1-y_hat1,'--m');grid on;
% % % legend('E1','E_y{1}');

% % subplot(224)
% % % plot(t,repmat(y_ref2,1,N)-y_hat_m2);hold on;
% % plot(t,y_ref2-y_hat2,'--m');grid on;

% % % legend('E2','E_y{2}');
