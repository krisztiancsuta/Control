function dxdt = nonlinear_model(x,u,psi_dot_des)
% x = [e1,de1,e2,de2,vx]
% u = [steering angle, force]
% d = disturbance
Caf = 222685.8 / 2; 
Car = 136242.8 / 2;
lf = 1236e-3;
lr = 2789e-3 - lf; 
m = 2300; 
Iz = 2873; 
b = 50;

vx = max(x(5), 0.1);

a22 = -(2*Caf+2*Car)/(m*vx);
a23 = (2*Caf+2*Car)/m;
a24 = (-2*lf*Caf+2*lr*Car)/(m*vx);
a42 = -(2*lf*Caf-2*lr*Car)/(Iz*vx);
a43 = (2*lf*Caf-2*lr*Car)/Iz;
a44 = -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*vx);
a55 = -b/m;

bc1_2 = 2*Caf/m;
bc1_4 = 2*lf*Caf/Iz;

bc2_5 = 1/m;

d2 = -((2*lf*Caf-2*lr*Car)/(m*vx))-vx;
d4 = -(2*Caf*lf^2 + 2*Car*lr^2) / (Iz * vx);

dxdt = zeros(5,1);
dxdt(1) = x(2); %e1
dxdt(2) = a22*x(2) + a23*x(3) + a24*x(4) + bc1_2*u(1)+ d2*psi_dot_des;%de1
dxdt(3) = x(4);%e2
dxdt(4) = a42*x(2) + a43*x(3) + a44*x(4) + bc1_4*u(1) + d4*psi_dot_des;%de2
dxdt(5) = a55*x(5) + bc2_5*u(2);%vx


end
