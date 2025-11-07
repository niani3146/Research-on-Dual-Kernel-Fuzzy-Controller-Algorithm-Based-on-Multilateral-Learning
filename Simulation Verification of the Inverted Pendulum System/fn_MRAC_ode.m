function dx = fn_MRAC_ode(t,x)
n=2;
dx=zeros(6+2*n+2*n+n,1);


A=[0,1;0,0];%倒立摆系统
B=[0;1];%倒立摆系统
Ar=[0,1;-1,-2];%参考系统
Br=[0;1];%参考系统

%激励信号
r=sin(t);
% % r=1;
% if t>=20 && t<25
%     r=1;
% elseif t>50 && t<65
%     r=3;
% else
%     r=0;
% end

kr=[-1;-2;1];
xre=[x(3);x(4);r];

ke=[50;30];
e=x(3)-x(1);
de=x(4)-x(2);
e1=[e;de];

Wd=[1;-1;0.5];
fphi=[exp(x(1)*x(2)),sin(x(1)),abs(x(2))*x(2)]';
% fphi=[sin(x(1)),cos(x(1)),tan(x(1))]';
% fphi=[sin(t),cos(t),1]';

Al=A-B*ke';
Q=10*eye(2);
P=lyap(Al,Q);

%二阶线性滤波
zeta=0.7;
omega=100;
en=B'*e1;
dx(5)=x(6);
dx(6)=-2*zeta*omega*x(6)+omega^2*(en-x(5));%二阶线性滤波求导数

y_hat=x(7+n:6+2*n)'*x(7:6+n);    %多边学习输出，系数*自适应输出           ！！！！！！！！
y_ref=-(x(6)-Al(2,1)*e-Al(2,2)*de)+y_hat;     %~f+fi=f，不确定项参考值

ee_yn=y_ref-x(7:6+n);          %自适应的逼近误差
% dx(7+4*n:6+5*n)=ee_yn;

zeta=0.7;
omega=100;
dx(12)=-2*zeta*omega*x(12)+omega^2*(en-x(11));
dx(7+2*n:6+3*n)=x(7+3*n:6+4*n);
dx(7+3*n:6+4*n)=-2*zeta*omega*x(7+3*n:6+4*n)+omega^2*(ee_yn-x(7+2*n:6+3*n));%二阶线性滤波求导数

e_y=y_ref-y_hat;         %多边学习的逼近误差
% dx(7+1)=e_y;

u=ke'*e1+kr'*xre-y_hat;           %控制器输出

gamma=5;%学习率
k_gain=5;%预测误差比例增益
dx(7:6+n)=k_gain*ee_yn;

up_1=gamma*e_y*x(7:6+n);  %多边学习的权重更新
dx(7+n:6+2*n)=sig_sat_line(up_1,5);%多边学习的权重更新值存储在dx中

dx(1)=x(2);
dx(2)=Wd'*fphi+u;
dx(3)=x(4);
dx(4)=-x(3)-2*x(4)+r;   %返回dx的计算结果
end