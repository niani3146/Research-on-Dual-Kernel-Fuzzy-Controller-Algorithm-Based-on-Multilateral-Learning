function dx = PID_test_ode(t,x)
dx=zeros(5,1);

A=[0,1;0,0];
B=[0;1];
Ar=[0,1;-1,-2];
Br=[0;1];
Wd=[1,-1,0.5]';
fphi=[exp(x(1)*x(2)),sin(x(1)),abs(x(2))*x(2)]';
dyn_true=Wd'*fphi;

e=x(3)-x(1);
de=x(4)-x(2);

xref=[sin(t);cos(t)];

Kp=[50,30];
Ki=50;
u=Kp*[e;de]+Ki*x(5);

dx(1)=x(2);
dx(2)=dyn_true+u;
dx(3)=x(4);
dx(4)=-x(3)-2*x(4)+xref(1);
dx(5)=e;

end
