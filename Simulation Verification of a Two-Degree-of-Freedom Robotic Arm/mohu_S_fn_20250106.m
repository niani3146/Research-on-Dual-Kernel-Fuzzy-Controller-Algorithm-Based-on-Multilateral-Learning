function dx = mohu_S_fn_20250106(t,x)
N=1;
n=5;

dx=zeros(2*6+n^4*N+n^4*N+n^4+n^4+2,1);

A=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
B=[0 0;0 0;1 0;0 1];
Ar=[0,0,1,0;0,0,0,1;-1,0,-2,0;0,-1,0,-2];
Br=[0 0;0 0;1 0;0 1];

r1=sin(t);r2=cos(t);

kr=[-1,0,-2,0,1,0;...
    0,-1,0,-2,0,1]';
xre=[x(5);x(6);x(7);x(8);r1;r2];

ke=[50,0,20,0;...
    0,50,0,20]';

e11=x(5)-x(1);
e12=x(7)-x(3);
e1=[e11;e12];
e21=x(6)-x(2);
e22=x(8)-x(4);
e2=[e21;e22];
e=[e11;e21;e12;e22];%e=[e,e']

a0=1;
a1=1;
m1=0.1;
m2=0.1;
gs=9.8;
M=[(a0^2*m1)/3 + a0^2*m2 + (a1^2*m2)/3 + a0*a1*m2*cos(x(2)),  (a1*m2*(2*a1 + 3*a0*cos(x(2))))/6;
    (a1*m2*(2*a1 + 3*a0*cos(x(2))))/6,  (a1^2*m2)/3];
C=[-(a0*a1*x(4)*m2*sin(x(2)))/2,   -(a0*a1*m2*sin(x(2))*(x(3) + x(4)))/2;
    (a0*a1*x(3)*m2*sin(x(2)))/2,   0];
G=[gs*m2*((a1*cos(x(1) + x(2)))/2 + a0*cos(x(1))) + (a0*gs*m1*cos(x(1)))/2;
    (a1*gs*m2*cos(x(1) + x(2)))/2];

if abs(det(M))<1e-3  
    pseudo_invM=M'*M\M'; 
    un_model=-pseudo_invM*(C*[x(3);x(4)]+G);
else
    un_model=-inv(M)*(C*[x(3);x(4)]+G);
end

y_true=un_model;

Al=A-B*ke';
Q=10*eye(4);
P=lyap(Al,Q);

zeta=0.7;
omega=100;
AA=[1,0;0,1];
BB=[0;1];
keke=[50,20]';
All=AA-BB*keke';
en1=BB'*e1;
en2=BB'*e2;

dx(9)=x(10);
dx(10)=-2*zeta*omega*x(10)+omega^2*(en1-x(9));
dx(11)=x(12);
dx(12)=-2*zeta*omega*x(12)+omega^2*(en2-x(11));

sigma1 =2;
centers1 = [-15, -7.5, 0, 7.4, 15];

membership.q1 = zeros(1, n);
membership.q2 = zeros(1, n);
membership.dq1 = zeros(1, n);
membership.dq2 = zeros(1, n);
for i = 1:n
    membership.q1(i) = exp(-((x(1) - centers1(i))^2) /(2*sigma1^2));
    membership.q2(i) = exp(-((x(2) - centers1(i))^2) /(2*sigma1^2));      
    membership.dq1(i) = exp(-((x(3) - centers1(i))^2) /(2*sigma1^2));
    membership.dq2(i) = exp(-((x(4) - centers1(i))^2) /(2*sigma1^2));
end

FS2=zeros(n^4,1);
FS1=0;
for L1=1:n
    for L2=1:n
        for L3=1:n
            for L4=1:n
              FS2(n^3*(L1-1)+n^2*(L2-1)+n*(L3-1)+L4)=membership.q1(L1)*membership.dq1(L2)'*membership.q2(L3)*membership.dq2(L4);
              FS1=FS1+membership.q1(L1)*membership.dq1(L2)'*membership.q2(L3)*membership.dq2(L4);
            end
        end
    end
end
xi_1 = FS2 / (FS1+0.0005);

xi_2=reshape(xi_1,n^4,1);                     

dx(12+(n^4)*N+(n^4)*N+1:12+(n^4)*N+(n^4)*N+n^4)=xi_2;

W_hat1=x(12+1:12+n^4);
y_hat_m1= W_hat1'* xi_2;

W_hat2=x(12+(n^4)+1:12+(n^4)+n^4);   
y_hat_m2= W_hat2'* xi_2;   

y_ref1=-(x(10)-All(2,1)*e11-All(2,2)*e12)+y_hat_m1;
y_ref2=-(x(12)-All(2,1)*e21-All(2,2)*e22)+y_hat_m2;

ee_yn1=y_ref1-y_hat_m1;
ee_yn2=y_ref2-y_hat_m2; 


u=ke'*e+kr'*xre-[y_hat_m1;y_hat_m2];

eta1 =100;         
eta2 =80; 
m_1=100;
m_2=60;  

dW_hat1=zeros(1,n^4);
dW_hat1 = eta1 * (ee_yn1'* xi_2') +dW_hat1 ;

for j=1:n^4
    if isnan(dW_hat1(j))
        dW_hat1(j)=0;
    end
end
for j=1:n^4
    if dW_hat1(j) >= m_1
        dW_hat1(j) = m_1;
    elseif dW_hat1(j) <= -m_1
        dW_hat1(j) = -m_1;
    end
end

dW_hat2=zeros(N,n^4);
dW_hat2(:) = eta2 * (ee_yn2'* xi_2') +dW_hat2;

for j=1:n^4
    if isnan(dW_hat2(j))
        dW_hat2(j)=0;
    end
end
for j=1:n^4
    if dW_hat2(j) >= m_2
        dW_hat2(j) = m_2;
    elseif dW_hat2(j) <= -m_2
        dW_hat2(j) = -m_2;
    end
end

dx(12+1:12+n^4)=dW_hat1;
dx(12+n^4+1:12+n^4+n^4)=dW_hat2;

dx(1)=x(3);
dx(2)=x(4);
dx(3)=y_true(1)+u(1);
dx(4)=y_true(2)+u(2);
dx(5)=x(7);
dx(6)=x(8);
dx(7)=-x(5)-2*x(7)+r1;
dx(8)=-x(6)-2*x(8)+r2;

end
