function dx = mohu_2_fn_20250228(t,x)
N=2;
n=3; 
m=2; 

% dx=zeros(2*6+2*N+(n^4)*N+(n^4)*N+(n^4)+(n^4)+2*N,1);
dx=zeros(2*6+2*N+(n^4)*N*m+(n^4)*N*m+(n^4)*N*m+2+N*m+8,1);
A=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
B=[0 0;0 0;1 0;0 1];
Ar=[0,0,1,0;0,0,0,1;-1,0,-2,0;0,-1,0,-2];
Br=[0 0;0 0;1 0;0 1];

r1=sin(t);r2=cos(t);
% if t>=20 && t<50
%     r1=1;r2=2;
% % elseif t>50 && t<65
% %     r1=3;r2=1.5;
% else
%     r1=0;r2=0;
% end

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

%fref=f
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

% sigma1 =0.6; 
% centers1 = [-6, 0, 6];
sigma1 =3;
centers1 = [-45, 0,45];
membership.q11 = zeros(1, n);
membership.q12 = zeros(1, n);
membership.dq11 = zeros(1, n);
membership.dq12 = zeros(1, n);
for i = 1:n
    membership.q11(i) = exp(-((x(1) - centers1(i))^2) /(2*sigma1^2));
    membership.q12(i) = exp(-((x(2) - centers1(i))^2) /(2*sigma1^2));       
    membership.dq11(i) = exp(-((x(3) - centers1(i))^2) /(2*sigma1^2));
    membership.dq12(i) = exp(-((x(4) - centers1(i))^2) /(2*sigma1^2));
end

sigma2 =8;
centers2 = [-25, 0, 25];
% centers2 = [-40,-20, 0, 20,40];
membership.q21=zeros(1, n);
membership.q22=zeros(1, n);
membership.dq21=zeros(1, n);
membership.dq22=zeros(1, n);
for i = 1:n
    membership.q21(i) = 1 /(1+sigma2*(x(1)-centers2(i))^2);
    membership.q22(i) = 1 /(1+sigma2*(x(2)-centers2(i))^2);      
    membership.dq21(i) = 1 /(1+sigma2*(x(3)-centers2(i))^2);
    membership.dq22(i) = 1 /(1+sigma2*(x(4)-centers2(i))^2);
end

params1=[1 6 14 20];
params2=[-10 -5 8 10];
params3=[10 18 23 30];
% params4=[5 8 12 1];
% params5=[5 8 12 1];
% params=[params1',params2',params3',params4',params5']';
params=[params1',params2',params3']';
for i = 1:n
    membership.q31(i) = pimf(x(1), params(i,:));
    membership.q32(i) = pimf(x(2), params(i,:));
    membership.dq31(i) = pimf(x(3), params(i,:));
    membership.dq32(i) = pimf(x(4), params(i,:));
end

sigma4 =8; 
centers4 = [-45, 0, 45];
membership.q41 = zeros(1, n);
membership.q42 = zeros(1, n);
membership.dq41 = zeros(1, n);
membership.dq42 = zeros(1, n);
for i = 1:n
    membership.q41(i) = exp(-((x(1) - centers4(i))^2) /(2*sigma4^2));
    membership.q42(i) = exp(-((x(2) - centers4(i))^2) /(2*sigma4^2));     
    membership.dq41(i) = exp(-((x(3) - centers4(i))^2) /(2*sigma4^2));
    membership.dq42(i) = exp(-((x(4) - centers4(i))^2) /(2*sigma4^2));
end

FS12=zeros(n^4,1);
FS11=0;
for L1=1:n
    for L2=1:n
        for L3=1:n
            for L4=1:n
              FS12(n^3*(L1-1)+n^2*(L2-1)+n*(L3-1)+L4)=membership.q11(L1)*membership.dq11(L2)'*membership.q12(L3)*membership.dq12(L4);
              FS11=FS11+membership.q11(L1)*membership.dq11(L2)'*membership.q12(L3)*membership.dq12(L4);
            end
        end
    end
end
xi_11 = FS12 / (FS11+0.0005);
xi_12=reshape(xi_11,n^4,1);   

FS22=zeros(n^4,1);
FS21=0;
for L1=1:n
    for L2=1:n
        for L3=1:n
            for L4=1:n
              FS22(n^3*(L1-1)+n^2*(L2-1)+n*(L3-1)+L4)=membership.q21(L1)*membership.dq21(L2)'*membership.q22(L3)*membership.dq22(L4);
              FS21=FS21+membership.q21(L1)*membership.dq21(L2)'*membership.q22(L3)*membership.dq22(L4);
            end
        end
    end
end
xi_21 = FS22 / (FS21+0.0005);
xi_22=reshape(xi_21,n^4,1);   

FS32=zeros(n^4,1);
FS31=0;
for L1=1:n
    for L2=1:n
        for L3=1:n
            for L4=1:n
              FS32(n^3*(L1-1)+n^2*(L2-1)+n*(L3-1)+L4)=membership.q31(L1)*membership.dq31(L2)'*membership.q32(L3)*membership.dq32(L4);
              FS31=FS31+membership.q31(L1)*membership.dq31(L2)'*membership.q32(L3)*membership.dq32(L4);
            end
        end
    end
end
xi_31 = FS32 / (FS31+0.0005);
xi_32=reshape(xi_31,n^4,1);   

FS42=zeros(n^4,1);
FS41=0;
for L1=1:n
    for L2=1:n
        for L3=1:n
            for L4=1:n
              FS42(n^3*(L1-1)+n^2*(L2-1)+n*(L3-1)+L4)=membership.q41(L1)*membership.dq41(L2)'*membership.q42(L3)*membership.dq42(L4);
              FS41=FS41+membership.q41(L1)*membership.dq41(L2)'*membership.q42(L3)*membership.dq42(L4);
            end
        end
    end
end
xi_41 = FS42 / (FS41+0.0005);
xi_42=reshape(xi_41,n^4,1);   

% xi_2=reshape(xi_1,n^4,1);               

dx(12+2*N+(n^4)*8+1:12+2*N+(n^4)*8+n^4)=xi_12; 
dx(12+2*N+(n^4)*9+1:12+2*N+(n^4)*9+n^4)=xi_22; 
dx(12+2*N+(n^4)*10+1:12+2*N+(n^4)*10+n^4)=xi_32;
dx(12+2*N+(n^4)*11+1:12+2*N+(n^4)*11+n^4)=xi_42;

% 第一条边的构成
W_hat_m11=x(12+2*N+1:12+2*N+n^4);
y_hat_m11= W_hat_m11(:)'* xi_12;                

W_hat_m12(:)=x(12+2*N+(n^4)*4+1:12+2*N+(n^4)*4+n^4);
y_hat_m12= W_hat_m12(:)'* xi_12;            
W_hat_m21(:)=x(12+2*N+(n^4)+1:12+2*N+(n^4)+n^4);
y_hat_m21= W_hat_m21(:)'* xi_22;                

W_hat_m22(:)=x(12+2*N+(n^4)*5+1:12+2*N+(n^4)*5+n^4);
y_hat_m22= W_hat_m22(:)'* xi_22;
% 
% y_hat_11=zeros(length(t),1);
% y_hat_12=zeros(length(t),1);


W_hat_m31(:)=x(12+2*N+(n^4)*2+1:12+2*N+(n^4)*2+n^4);
y_hat_m31= W_hat_m31(:)'* xi_32;                        

W_hat_m32(:)=x(12+2*N+(n^4)*6+1:12+2*N+(n^4)*6+n^4);
y_hat_m32= W_hat_m32(:)'* xi_32;                         

W_hat_m41(:)=x(12+2*N+(n^4)*3+1:12+2*N+(n^4)*3+n^4);
y_hat_m41= W_hat_m41(:)'* xi_42;                            

W_hat_m42(:)=x(12+2*N+(n^4)*7+1:12+2*N+(n^4)*7+n^4);
y_hat_m42= W_hat_m42(:)'* xi_42;                            

%单边双核输出 融合
y_hat_11=y_hat_m11+y_hat_m21;               
y_hat_12=y_hat_m12+y_hat_m22;                 

y_hat_21=y_hat_m31+y_hat_m41;               
y_hat_22=y_hat_m32+y_hat_m42;                 

y_hat_1=[y_hat_m11,y_hat_m21,y_hat_m31,y_hat_m41];         
y_hat_2=[y_hat_m12,y_hat_m22,y_hat_m32,y_hat_m42];          
xi_1=[xi_12,xi_22,xi_32,xi_42];

y_hat_m1=[y_hat_11,y_hat_21];             
y_hat_m2=[y_hat_12,y_hat_22];

y_hat1=x(12+1:12+N)'*y_hat_m1'; 
y_hat2=x(12+N+1:12+N+N)'*y_hat_m2'; 

y_ref1=-(x(10)-All(2,1)*e11-All(2,2)*e12)+y_hat1;
y_ref2=-(x(12)-All(2,1)*e21-All(2,2)*e22)+y_hat2;

dx(12+2*N+(n^4)*N*m+(n^4)*N*m+(n^4)*N*m+1)=y_ref1;
dx(12+2*N+(n^4)*N*m+(n^4)*N*m+(n^4)*N*m+2)=y_ref2;         

ee_yn_1=y_ref1-y_hat_1;       
ee_yn_2=y_ref2-y_hat_2;   

e_y1=y_ref1-y_hat1; 
e_y2=y_ref2-y_hat2; 

u=ke'*e+kr'*xre-[y_hat1;y_hat2];
eta1 =[40,40,40,50]; eta2 =[50,50,60,50];
m_1=[32,32,32,42];    m_2=[40,40,50,40];

dW_hat_1=zeros(N*m,n^4);
for i=1:N*m
   dW_hat_1(i,:) = eta1(i) * (ee_yn_1(i)* xi_1(:,i)')  ;
end
% dW_hat_1=[dW_hat_11',dW_hat_12']';
for i=1:N*m
    for j=1:(n^4)
        if isnan(dW_hat_1(i,j))
            dW_hat_1(i,j)=0;
        end
    end
end

for i=1:N*m
    for j=1:(n^4)
        if dW_hat_1(i,j) >= m_1(i)
            dW_hat_1(i,j) = m_1(i);
        elseif dW_hat_1(i,j) <= -m_1(i)
            dW_hat_1(i,j) = -m_1(i);
        end
    end
end

dW_hat_2=zeros(N*m,n^4);
for i=1:N*m
   dW_hat_2(i,:) = eta2(i) * (ee_yn_2(i)* xi_1(:,i)')  ;
end

for i=1:N*m
    for j=1:n^4
        if isnan(dW_hat_2(i,j))
            dW_hat_2(i,j)=0;
        end
    end
end

for i=1:N*m
    for j=1:n^4
        if dW_hat_2(i,j) >= m_2(i)
            dW_hat_2(i,j) = m_2(i);
        elseif dW_hat_2(i,j) <= -m_2(i)
            dW_hat_2(i,j) = -m_2(i);
        end
    end
end
% dx(rowi:coli)=dW_hat11;
for i=1:N*m
         rowi = 12+2*N+(i-1)*(n^4)+1;coli =12+2*N+(i-1)*(n^4)+n^4;
         dx(rowi:coli)=dW_hat_1(i,:);
end
for i=1:N*m
         rowi = 12+2*N+(n^4)*N*m+(i-1)*(n^4)+1;coli =12+2*N+(n^4)*N*m+(i-1)*(n^4)+n^4;
         dx(rowi:coli)=dW_hat_2(i,:);
end

% alpha1=[1,1];
% alpha2=[1,1];
% dW_hat1=zeros(N,1);
% for i=1:N
%    dW_hat1(i,:) = alpha1(i) * (ee_m1(i)* y_hat_m1(:,i)')  ;
% end
% dW_hat2=zeros(N,1);
% for i=1:N
%    dW_hat2(i,:) = alpha2(i) * (ee_m2(i)* y_hat_m2(:,i)')  ;
% end
% dx(12+2*N+(n^4)*N*m+(n^4)*N*m+(n^4)*N*m+2+4+1:12+2*N+(n^4)*N*m+(n^4)*N*m+(n^4)*N*m+2+4+2)=dW_hat1;
% dx(12+2*N+(n^4)*N*m+(n^4)*N*m+(n^4)*N*m+2+4+3:12+2*N+(n^4)*N*m+(n^4)*N*m+(n^4)*N*m+2+4+4)=dW_hat2;

gamma1= 0.8;
gamma2= 0.4;

up_1=gamma1*e_y1*y_hat_m1(1:N);
up_2=gamma2*e_y2*y_hat_m2(1:N);
% up_1=gamma1*e_y1*y_1;
% up_2=gamma2*e_y2*y_2;

dx(12+1:12+N)=sig_sat_line(up_1,0.1);
dx(12+N+1:12+N+N)=sig_sat_line(up_2,0.1);
% dx(12+1)=sig_sat_line(up_1,0.1);
% dx(12+N+1)=sig_sat_line(up_2,0.1);

dx(1)=x(3);
dx(2)=x(4);
dx(3)=y_true(1)+u(1);
dx(4)=y_true(2)+u(2);
dx(5)=x(7);
dx(6)=x(8);
dx(7)=-x(5)-2*x(7)+r1;
dx(8)=-x(6)-2*x(8)+r2;
end
