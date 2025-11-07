 function dx = mohu_fn_20241107_5(t,x)
% function [dx, y_ref] = mohu_fn_20241107_5(t, x)
n=7;
dx=zeros(6+n+25*n+25+2*n,1);

A=[0,1;0,0];%倒立摆系统
B=[0;1];%倒立摆系统
% Ar=[0,1;-1,-2];%参考系统
% Br=[0;1];%参考系统
%激励信号
r=sin(t);
% if t>=10 && t<15
%     r=1;
% elseif t>20 && t<25
%     r=-1.5;
% else
%     r=0;
% end

ke=[50;30];
kr=[-1;-2;1];

Wd=[1;-1;0.5];
fphi=[exp(x(1)*x(2)),sin(x(1)),abs(x(2))*x(2)]';

% 模糊控制
x_ref=[sin(t),cos(t)];

xi1 = [x(1);x(2)]; % 位置、速度（角度、角速度）信息
xi2 = [x(3);x(4)];
xi = [xi1;xi2]; 

e1=x(3)-x(1);
e2=x(4)-x(2); %e2为e1的导数
e=[e1;e2];

xre=[x(3);x(4);x_ref(:,1)];

sigma1 =0.4; % 高斯函数的标准差,增大标准差，提高灵敏度
sigma2 =0.3; 

centers1 = [-5, -2.5, 0, 2.5, 5];% 隶属度中心
% centers2 = [-3, -1.5, 0, 1.5, 3];
% centers1 = [-pi, -pi/2, 0, pi/2, pi]; 

% 计算高斯隶属度
membership.q1 = zeros(1, 5);
membership.q2 = zeros(1, 5);
for i = 1:5
    membership.q1(i) = exp(-((x(1) - centers1(i))^2) /(2*sigma1^2));
    membership.q2(i) = exp(-((x(2) - centers1(i))^2) /(2*sigma2^2));
%     membership.q1(i) = exp(-((x(1) - centers1(i))^2) /( (2*sigma1)^2));
%     membership.q2(i) = exp(-((x(2) - centers1(i))^2) /( (2*sigma2)^2));
end
% 归一化激活强度向量
FS2=zeros(5,5);
FS1=0;
for i=1:5
    for j=1:5
        FS2(i,j)=membership.q1(i)*membership.q2(j)';
        FS1=FS1+membership.q1(i)*membership.q2(j)';
    end
end
xi_1 = FS2 / (FS1+0.0005);
% xi_1 = FS2 / FS1;
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

%二阶线性滤波
zeta=0.7;
omega=100;
en=B'*e;
dx(5)=x(6);          %x(5)、x(6)：跟踪误差e及其导数
dx(6)=-2*zeta*omega*x(6)+omega^2*(en-x(5));

y_hat=x(7:6+n)'*y_hat_m(:,1:7)';  %x(7+n:6+2*n)'表示参数k,多边学习输出
y_ref=-(x(6)-Al(2,1)*e1-Al(2,2)*e2)+y_hat;%根据 二次滤波得到不确定项误差，模型不确定性参考误差：f(t)+y

ee_yn=y_ref-y_hat_m;   %模糊控制逼近误差
dx(214:220)=ee_yn;
e_y=y_ref-y_hat;     %多边学习逼近误差
dx(221:227)=e_y;

% eta = 45;
eta =40; % 学习率
dW_hat1=zeros(7,25);
for i=1:7
   dW_hat1(i,:) = eta * (ee_yn(i)'* xi_2') +dW_hat1(i,:) ;
%    gamma1=20;
%    dW_hat1(i,:) = eta * ee_yn(i)'* xi_2'+gamma1*e'*P*B*xi_2';
%    dW_hat1(i,:) = gamma*e'*P*B*xi_1;
end
for i=1:7
    for j=1:25
        if isnan(dW_hat1(i,j))
            dW_hat1(i,j)=0;
        end
    end
end
% 权重更新率限幅
% for i=1:7
%     for j=1:25
%         if dW_hat1(i,j) >= 1
%             dW_hat1(i,j) = 1;
%         elseif dW_hat1(i,j) <= -1
%             dW_hat1(i,j) = -1;
%         end
%     end
% end
for i=1:7
    for j=1:25
        if dW_hat1(i,j) >= 1
            dW_hat1(i,j) = 1;
        elseif dW_hat1(i,j) <= -1
            dW_hat1(i,j) = -1;
        end
    end
end
%     dW_hat1 = min(max(dW_hat1, -1), 1); % 限幅

for i=1:7
         rowi = 6+n+(i-1)*25+1;coli =6+n+(i-1)*25+25;
         dx(rowi:coli)=dW_hat1(i,:);
end

u=ke'*e+kr'*xre-y_hat;%控制器设计

% gamma= 2.5;
% up_1=gamma*e_y*y_hat_m(1:7);
% dx(7:6+n)=sig_sat_line(up_1,5);
gamma= 10;
up_1=gamma*e_y*y_hat_m(1:7);
dx(7:6+n)=sig_sat_line(up_1,1);

dx(1)=x(2);
dx(2)=Wd'*fphi+u;
dx(3)=x(4);
dx(4)=-x(3)-2*x(4)+r;
end