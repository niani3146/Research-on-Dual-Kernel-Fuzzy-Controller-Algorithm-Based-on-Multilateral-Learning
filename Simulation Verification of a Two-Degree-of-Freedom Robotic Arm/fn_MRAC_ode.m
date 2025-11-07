function dx = fn_MRAC_ode(t,x)
n=2;
dx=zeros(2*6+4*n+2*n+2,1);

A=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];%RR-bot系统
B=[0 0;0 0;1 0;0 1];%RR-bot系统
Ar=[0,0,1,0;0,0,0,1;-1,0,-2,0;0,-1,0,-2];%参考系统
Br=[0 0;0 0;1 0;0 1];%参考系统

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

%RR-bot系统参数
a0=1;%连杆1长度为1
a1=1;%连杆2长度为1
m1=0.1;%连杆1的质量为0.1kg
m2=0.1;%连杆2的质量为0.2kg
gs=9.8;%重力加速度
M=[(a0^2*m1)/3 + a0^2*m2 + (a1^2*m2)/3 + a0*a1*m2*cos(x(2)),  (a1*m2*(2*a1 + 3*a0*cos(x(2))))/6;
    (a1*m2*(2*a1 + 3*a0*cos(x(2))))/6,  (a1^2*m2)/3];
C=[-(a0*a1*x(4)*m2*sin(x(2)))/2,   -(a0*a1*m2*sin(x(2))*(x(3) + x(4)))/2;
    (a0*a1*x(3)*m2*sin(x(2)))/2,   0];
G=[gs*m2*((a1*cos(x(1) + x(2)))/2 + a0*cos(x(1))) + (a0*gs*m1*cos(x(1)))/2;
    (a1*gs*m2*cos(x(1) + x(2)))/2];

%计算系统不确定项f
if abs(det(M))<1e-3   %计算矩阵 M 的行列式，并取其绝对值，值小于 1e-3，认为矩阵 M接近奇异
    pseudo_invM=M'*M\M';   %矩阵M是 奇异矩阵，计算伪逆
    un_model=-pseudo_invM*(C*[x(3);x(4)]+G);
else
    un_model=-inv(M)*(C*[x(3);x(4)]+G);
end

%fref=f
y_true=un_model;

Al=A-B*ke';
Q=10*eye(4);
P=lyap(Al,Q);

%二阶线性滤波
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
dx(12)=-2*zeta*omega*x(12)+omega^2*(en2-x(11));%二阶滤波的误差及其导数的输出

% 自适应输出
y_hat1=x(13+n:12+2*n)'*x(13:12+n);
y_hat2=x(13+3*n:12+4*n)'*x(13+2*n:12+3*n);
y_ref1=-(x(10)-All(2,1)*e11-All(2,2)*e12)+y_hat1; 
y_ref2=-(x(12)-All(2,1)*e21-All(2,2)*e22)+y_hat2;
dx(2*6+4*n+2*n+1)=y_ref1;
dx(2*6+4*n+2*n+2)=y_ref2;

ee_yn1=y_ref1-x(13:12+n);
ee_yn2=y_ref2-x(13+2*n:12+3*n);

dx(13+4*n:12+5*n)=ee_yn1;
dx(13+5*n:12+6*n)=ee_yn2;

e_y1=y_ref1-y_hat1;
e_y2=y_ref2-y_hat2;

u=ke'*e+kr'*xre-[y_hat1;y_hat2];

gamma1=2;
gamma2=2;
k_gain1=5;%预测误差比例增益
k_gain2=5;%预测误差比例增益
dx(13:12+n)=k_gain1*ee_yn1;
dx(13+2*n:12+3*n)=k_gain2*ee_yn2;

up_1=gamma1*e_y1*x(13:12+n);
up_2=gamma2*e_y2*x(13+2*n:12+3*n);
dx(13+n:12+2*n)=sig_sat_line(up_1,0.1);
dx(13+3*n:12+4*n)=sig_sat_line(up_2,0.1);

dx(1)=x(3);
dx(2)=x(4);
dx(3)=y_true(1)+u(1);
dx(4)=y_true(2)+u(2);
dx(5)=x(7);
dx(6)=x(8);
dx(7)=-x(5)-2*x(7)+r1;
dx(8)=-x(6)-2*x(8)+r2;
end