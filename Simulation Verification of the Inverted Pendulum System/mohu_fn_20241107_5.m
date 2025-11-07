function dx = mohu_fn_20241107_5(t,x)
n=7;
dx=zeros(6+n+25*n+25+2*n,1);

A=[0,1;0,0];
B=[0;1];
% Ar=[0,1;-1,-2];
% Br=[0;1];

r=sin(t);

ke=[50;30];
kr=[-1;-2;1];

Wd=[1;-1;0.5];
fphi=[exp(x(1)*x(2)),sin(x(1)),abs(x(2))*x(2)]';

x_ref=[sin(t),cos(t)];

xi1 = [x(1);x(2)]; 
xi2 = [x(3);x(4)];
xi = [xi1;xi2]; 

e1=x(3)-x(1);
e2=x(4)-x(2);
e=[e1;e2];

xre=[x(3);x(4);x_ref(:,1)];

sigma1 =0.4; 
sigma2 =0.3; 

centers1 = [-5, -2.5, 0, 2.5, 5];

membership.q1 = zeros(1, 5);
membership.q2 = zeros(1, 5);
for i = 1:5
    membership.q1(i) = exp(-((x(1) - centers1(i))^2) /(2*sigma1^2));
    membership.q2(i) = exp(-((x(2) - centers1(i))^2) /(2*sigma2^2));
end

FS2=zeros(5,5);
FS1=0;
for i=1:5
    for j=1:5
        FS2(i,j)=membership.q1(i)*membership.q2(j)';
        FS1=FS1+membership.q1(i)*membership.q2(j)';
    end
end
xi_1 = FS2 / (FS1+0.0005);
xi_2=reshape(xi_1,25,1);
for i=1:5
         rowi = 6+n+25*n+(i-1)*5+1;coli =6+n+25*n+(i-1)*5+5;
         dx(rowi:coli)=xi_1(i,:);
end
Al=A-B*ke';
Q=10*eye(2);
P=lyap(Al,Q);

y_hat_m=zeros(1,7);
for i=1:7
     W_hat1(i,:)=x(6+n+(i-1)*25+1:6+n+(i-1)*25+25);
     y_hat_m(i)= W_hat1(i,:)* xi_2;     
end

zeta=0.7;
omega=100;
en=B'*e;
dx(5)=x(6);          
dx(6)=-2*zeta*omega*x(6)+omega^2*(en-x(5));

y_hat=x(7:6+n)'*y_hat_m(:,1:7)';  
y_ref=-(x(6)-Al(2,1)*e1-Al(2,2)*e2)+y_hat;

ee_yn=y_ref-y_hat_m;  
dx(214:220)=ee_yn;
e_y=y_ref-y_hat;   
dx(221:227)=e_y;

eta =40;
dW_hat1=zeros(7,25);
for i=1:7
   dW_hat1(i,:) = eta * (ee_yn(i)'* xi_2') +dW_hat1(i,:) ;
end
for i=1:7
    for j=1:25
        if isnan(dW_hat1(i,j))
            dW_hat1(i,j)=0;
        end
    end
end

for i=1:7
    for j=1:25
        if dW_hat1(i,j) >= 1
            dW_hat1(i,j) = 1;
        elseif dW_hat1(i,j) <= -1
            dW_hat1(i,j) = -1;
        end
    end
end

for i=1:7
         rowi = 6+n+(i-1)*25+1;coli =6+n+(i-1)*25+25;
         dx(rowi:coli)=dW_hat1(i,:);
end

u=ke'*e+kr'*xre-y_hat;
gamma= 10;
up_1=gamma*e_y*y_hat_m(1:7);
dx(7:6+n)=sig_sat_line(up_1,1);

dx(1)=x(2);
dx(2)=Wd'*fphi+u;
dx(3)=x(4);
dx(4)=-x(3)-2*x(4)+r;

end
