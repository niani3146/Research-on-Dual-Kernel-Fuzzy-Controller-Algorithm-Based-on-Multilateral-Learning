function dx = PID_RRbot_ode(t,x)
dx=zeros(10,1);
A=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
B=[0 0;0 0;1 0;0 1];
Ar=[0,0,1,0;0,0,0,1;-1,0,-2,0;0,-1,0,-2];
Br=[0 0;0 0;1 0;0 1];


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
fphi=[un_model(1);un_model(2)];
xref=[sin(t),cos(t),cos(t),-sin(t)];

e1=x(5)-x(1);
e2=x(6)-x(2);

Kp=[500,500];Kd=[100,100];Ki=[200,200];
u1=Kp(1)*(x(5)-x(1))+Kd(1)*(x(7)-x(3))+Ki(1)*x(9);
u2=Kp(2)*(x(6)-x(2))+Kd(2)*(x(8)-x(4))+Ki(2)*x(10);

dx(1)=x(3);
dx(2)=x(4);
dx(3)=fphi(1)+u1;
dx(4)=fphi(2)+u2;
dx(5)=x(7);
dx(6)=x(8);
dx(7)=-x(5)-2*x(7)+xref(1);
dx(8)=-x(6)-2*x(8)+xref(2);
dx(9)=e1;
dx(10)=e2;

end
