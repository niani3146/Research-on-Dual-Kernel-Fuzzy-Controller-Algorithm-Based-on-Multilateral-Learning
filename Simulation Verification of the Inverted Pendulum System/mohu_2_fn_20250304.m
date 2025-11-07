function dx = mohu_2_fn_20250304(t,x)
N=2;
n=3; 
m=2;

% dx=zeros(2*6+2*N+(n^4)*N+(n^4)*N+(n^4)+(n^4)+2*N,1);
dx=zeros(6+N+(n^2)*N*m+(n^2)*N*m+1,1);
A=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
B=[0 0;0 0;1 0;0 1];
Ar=[0,0,1,0;0,0,0,1;-1,0,-2,0;0,-1,0,-2];
Br=[0 0;0 0;1 0;0 1];

r=sin(t);
kr=[-1,-2,1];

ke=[50,20];

Wd=[1;-1;0.5];
fphi=[exp(x(1)*x(2)),sin(x(1)),abs(x(2))*x(2)]';

x_ref=[sin(t),cos(t)];

e11=x(3)-x(1);
e12=x(4)-x(2); 
e1=[e11;e12];
e=e1;

xre=[x(3);x(4);x_ref(:,1)];

Al=A-B*ke';
Q=10*eye(4);
P=lyap(Al,Q);

zeta=0.7;
omega=100;
AA=[0,1;0,0];
keke=[50,20]';
All=AA-BB*keke';
en=BB'*e1;

dx(5)=x(6);          
dx(6)=-2*zeta*omega*x(6)+omega^2*(en-x(5));

sigma1 =3;
centers1 = [-45, 0,45];
membership.q11 = zeros(1, n);
membership.q12 = zeros(1, n);
for i = 1:n
    membership.q11(i) = exp(-((x(1) - centers1(i))^2) /(2*sigma1^2));
    membership.q12(i) = exp(-((x(2) - centers1(i))^2) /(2*sigma1^2));      
end

sigma2 =8;
centers2 = [-30, 0,30];

membership.q21=zeros(1, n);
membership.q22=zeros(1, n);
for i = 1:n
    membership.q21(i) = 1 /(1+sigma2*(x(1)-centers2(i))^2);
    membership.q22(i) = 1 /(1+sigma2*(x(2)-centers2(i))^2);      
end

params1=[1 6 14 20];
params2=[-10 -5 8 10];
params3=[10 18 23 30];
params=[params1',params2',params3']';
for i = 1:n
    membership.q31(i) = pimf(x(1), params(i,:));
    membership.q32(i) = pimf(x(2), params(i,:));
end

sigma4 =4; 
centers4 = [-30, 0, 30];
membership.q41 = zeros(1, n);
membership.q42 = zeros(1, n);
for i = 1:n
    membership.q41(i) = exp(-((x(1) - centers4(i))^2) /(2*sigma4^2));
    membership.q42(i) = exp(-((x(2) - centers4(i))^2) /(2*sigma4^2));     
end

FS12=zeros(n,n);
FS11=0;
for i=1:n
    for j=1:n
        FS12(i,j)=membership.q11(i)*membership.q12(j)';
        FS11=FS11+membership.q11(i)*membership.q12(j)';
    end
end
xi_11 = FS12 / (FS11+0.0005);
xi_12=reshape(xi_11,n^2,1);   

FS22=zeros(n,n);
FS21=0;
for i=1:n
    for j=1:n
        FS22(i,j)=membership.q21(i)*membership.q22(j)';
        FS21=FS21+membership.q21(i)*membership.q22(j)';
    end
end
xi_21 = FS22 / (FS21+0.0005);
xi_22=reshape(xi_21,n^2,1);   

FS32=zeros(n,n);
FS31=0;
for i=1:n
    for j=1:n
        FS32(i,j)=membership.q31(i)*membership.q32(j)';
        FS31=FS31+membership.q31(i)*membership.q32(j)';
    end
end
xi_31 = FS32 / (FS31+0.0005);
xi_32=reshape(xi_31,n^2,1);   

FS42=zeros(n,n);
FS41=0;
for i=1:n
    for j=1:n
        FS42(i,j)=membership.q41(i)*membership.q42(j)';
        FS41=FS41+membership.q41(i)*membership.q42(j)';
    end
end
xi_41 = FS42 / (FS41+0.0005);
xi_42=reshape(xi_41,n^2,1);   

dx(6+N+(n^2)*4+1:6+N+(n^2)*4+n^2)=xi_12; 
dx(6+N+(n^2)*5+1:6+N+(n^2)*5+n^2)=xi_22;           
dx(6+N+(n^2)*6+1:6+N+(n^2)*6+n^2)=xi_32;
dx(6+N+(n^2)*7+1:6+N+(n^2)*7+n^2)=xi_42;

W_hat_m11=x(6+N+1:6+N+n^2);
y_hat_m11= W_hat_m11(:)'* xi_12;              

W_hat_m21(:)=x(6+N+(n^2)+1:6+N+(n^2)+n^2);
y_hat_m21= W_hat_m21(:)'* xi_22;                   

W_hat_m31(:)=x(6+N+(n^2)*2+1:6+N+(n^2)*2+n^2);
y_hat_m31= W_hat_m31(:)'* xi_32;                        
W_hat_m41(:)=x(6+N+(n^2)*3+1:6+N+(n^2)*3+n^2);
y_hat_m41= W_hat_m41(:)'* xi_42;                            

y_hat_11=y_hat_m11+y_hat_m21;              
y_hat_21=y_hat_m31+y_hat_m41;                 

y_hat_1=[y_hat_m11,y_hat_m21,y_hat_m31,y_hat_m41];             

xi_1=[xi_12,xi_22,xi_32,xi_42];

y_hat_m1=[y_hat_11,y_hat_21];             

y_hat1=x(6+1:6+N)'*y_hat_m1'; 

y_ref1=-(x(6)-All(2,1)*e11-All(2,2)*e12)+y_hat1;

dx(6+N+(n^2)*N*m+(n^2)*N*m+1)=y_ref1;

ee_yn_1=y_ref1-y_hat_1;       

e_y1=y_ref1-y_hat1; 

u=ke*e+kr*xre-y_hat1;

eta1 =[60,60,60,60];
m_1=[42,42,42,42]; 

dW_hat_1=zeros(N*m,n^2);
for i=1:N*m
   dW_hat_1(i,:) = eta1(i) * (ee_yn_1(i)* xi_1(:,i)')  ;
end

for i=1:N*m
    for j=1:(n^2)
        if isnan(dW_hat_1(i,j))
            dW_hat_1(i,j)=0;
        end
    end
end

for i=1:N*m
    for j=1:(n^2)
        if dW_hat_1(i,j) >= m_1(i)
            dW_hat_1(i,j) = m_1(i);
        elseif dW_hat_1(i,j) <= -m_1(i)
            dW_hat_1(i,j) = -m_1(i);
        end
    end
end

for i=1:N*m
         rowi = 6+N+(i-1)*(n^2)+1;coli =6+N+(i-1)*(n^2)+n^2;
         dx(rowi:coli)=dW_hat_1(i,:);
end

gamma1=7;

up_1=gamma1*e_y1*y_hat_m1(1:N);

dx(6+1:6+N)=sig_sat_line(up_1,1);

dx(1)=x(2);
dx(2)=Wd'*fphi+u;
dx(3)=x(4);
dx(4)=-x(3)-2*x(4)+r;
end
