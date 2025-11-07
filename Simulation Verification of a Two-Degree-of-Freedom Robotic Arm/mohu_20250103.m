% clear all;
% close all;
N=7;
n=3;

t0=0;tf=30;
T=linspace(t0,tf,100);

para_a=linspace(-0.1,0.2,N);
para_b=linspace(-0.1,0.2,N);
% for i = 1:N
% %     min_val = 0.08 + (i - 1) * 0.015;
% %     max_val = 0.08 + (i - 1) * 0.03; 
%     min_val = 0.0001 + (i - 1) * 0.0001;
%     max_val = 0.0001 + (i - 1) * 0.0013;
%     para_k((i-1)*n^4+1:(i-1)*n^4+n^4) = linspace(min_val, max_val, n^4)';
% end
para_k=zeros(1,N*(n^4));
% for i = 1:N
%     min_val = 0.008 + (i - 1) * 0.0008; 
%     max_val = 0.008 + (i - 1) * 0.0013; 
%     para_w((i-1)*n^4+1:(i-1)*n^4+n^4) = linspace(min_val, max_val, n^4)'; 
% end
% para_w=linspace(-0.01,0,N*n^4);
% para_k=ones(1,N*n^4);
para_w=zeros(1,N*(n^4));
para_n=zeros(1,n^4);
para_m=zeros(1,n^4);
para_c=zeros(1,2*N);

x0=[1,1,1,1,1,1,1,1,0,0,0,0,... 
  para_a,para_b,para_k,para_w,para_n,para_m,para_c];

options = odeset('MaxStep', 0.01);
% options = odeset('MaxStep', 0.02); 
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
        y_ref3(j,:)=x(j,12+2*N+(n^4)*N+(n^4)*N+2*(n^4)+1);
        y_ref4(j,:)=x(j,12+2*N+(n^4)*N+(n^4)*N+2*(n^4)+2);
end

u_tau=zeros(2,length(t));

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

y_hat1=zeros(length(t),1);
for j=1:length(t) 
     y_hat1(j,1)=x(j,12+1:12+N)*y_hat_m1(j,1:N)';
end
y_hat2=zeros(length(t),1);
for j=1:length(t) 
     y_hat2(j,1)=x(j,12+N+1:12+N+N)*y_hat_m2(j,1:N)';
end

for i=1:length(t)
    u_tau(:,i)=ke'*e(i,:)'+kr'*xre(i,:)'-[y_hat1(i);y_hat2(i)];
end

disp(['state 1:',num2str(sum(abs(e11)))])
disp(['state 2:',num2str(sum(abs(e21)))])
disp(['state 3:',num2str(sum(abs(e12)))])
disp(['state 4:',num2str(sum(abs(e22)))])
disp(['state 1:',num2str(sum(abs(e11))/length(t))])
disp(['state 2:',num2str(sum(abs(e21))/length(t))])
disp(['state 3:',num2str(sum(abs(e12))/length(t))])
disp(['state 4:',num2str(sum(abs(e22))/length(t))])
disp(['state 5:',num2str(sum(abs(y_ref1-y_hat1)))])
disp(['state 6:',num2str(sum(abs(y_ref2-y_hat2)))])

Xref3_1=x(:,5);
X5_1=x(:,1);
X5_2=x(:,3);

Xref3_2=x(:,7);
X5_3=x(:,2);
X5_4=x(:,4);

E5_11=e11;
E5_12=e12;
E5_21=e21;
E5_22=e22;

U5=u_tau;

yref5_1=y_true(1,:)';
yref5_2=y_true(2,:)';
Y_ref=[yref5_1,yref5_2];

Y_hat5_1=y_hat1;
Y_hat5_2=y_hat2;

EE5_1=y_ref1-y_hat1;
EE5_2=y_ref2-y_hat2;
