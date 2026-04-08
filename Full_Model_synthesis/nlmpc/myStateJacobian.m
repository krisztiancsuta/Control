function [A,Bmv] = myStateJacobian(x,u,psi_dot_des)
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

%d2 = -((2*lf*Caf-2*lr*Car)/(m*vx))-vx;
d4 = -(2*Caf*lf^2 + 2*Car*lr^2) / (Iz * vx);

% Calculate Non-linear Partial Derivatives with respect to vx (x5)
% Since a22, a24, a42, a44, d2, and d4 depend on 1/vx, their derivative
% with respect to vx introduces a -1/vx^2 factor.
d2_const = -(2*lf*Caf - 2*lr*Car) / m;

%d(de1)/dvx
dde1_dvx = -(a22 / vx) * x(2) - (a24 / vx) * x(4) + ( -d2_const / (vx^2) - 1 ) * psi_dot_des;
%d(de2)/dvx
dde2_dvx = -(a42 / vx) * x(2) - (a44 / vx) * x(4) - (d4 / vx) * psi_dot_des;

% (If vx was bounded to 0.1, the gradient with respect to x(5) is strictly 0)
if x(5) < 0.1
    dde1_dvx = 0;
    dde2_dvx = 0;
end

% Assemble Jacobian Matrix A (df/dx) - 5x5
A = [0,   1,   0,   0,         0;
     0, a22, a23, a24,   dde1_dvx;
     0,   0,   0,   1,         0;
     0, a42, a43, a44,   dde2_dvx;
     0,   0,   0,   0,       a55];

% Assemble Jacobian Matrix Bmv (df/du) - 5x2
Bmv = [0,           0;
       bc1_2,       0;
       0,           0;
       bc1_4,       0;
       0,       bc2_5];


end


