% Predictive Control within One Optimization Window
% Y= Fx(ki) + ΦΔU,

% Optimization
% Rs'=[1 1. . . 1] * r(ki), setpoint signals for Np period
% J = (Rs− Y )'* (Rs− Y ) + ΔU'RΔU, R is a diagonal matrix
% R = [ rw  0 ;... rw is a tuning parameter
%       0   rw;...]
% first term is linked to the objective of minimizin g the errors between
% the predicted output and the set-point signal while the second term reflects
% the consideration given to the size of ΔU when the objective function J is 
% made to be as small as possible.

% Goal is to minimize errors between output and set-point signal
% xm(k + 1) = axm(k) + bu(k)
% y(k) = xm(k)
a = 0.8;
b = 0.1;
c = 1;
Np = 10; % Prediction horizon
Nc = 4;
% Augmented state-space model
A = [a   0;...
     c*a 1];
B = [b;
     c*b];
C = [0 1]; 

%% Y= Fx(ki) + ΦΔU,

[n, m] = size(B); % n: állapotok száma, m: bemenetek száma
[q, ~] = size(C); % q: kimenetek száma
% F
F = zeros(Np*q, n);
for i = 1:Np
    F((i-1)*q + 1 : i*q, :) = C * A^i;
end

% Φ
% 2. Phi mátrix felépítése (Np*q x Nc*m)
Phi = zeros(Np*q, Nc*m);
for i = 1:Np
    for j = 1:Nc
        if i >= j
            Phi((i-1)*q + 1 : i*q, (j-1)*m + 1 : j*m) = C * A^(i-j) * B;
        end
    end
end


% Φ'*Φ,
Phi_Phi = Phi'*Phi;
% Φ'*F 
Phi_F = Phi'*F;
% Φ'*Rs
Rs = ones(Np, 1);
Phi_Rs = Phi'*Rs;

% At time ki = 10, the state vector x(ki) = [0.1 0.2]T .

% ΔU = inv((Φ'Φ))(Φ'*Rs− Φ'*F*x(ki)) 
% Mivel rw = 0 
dU = inv(Phi_Phi)*(Phi_Rs-Phi_F*[0.1; 0.2])

% Mivel rw = 10
rw = 10;
R = diag([rw rw rw rw]);
dU = inv(Phi_Phi + R)*(Phi_Rs-Phi_F*[0.1; 0.2])

%% Todo make it as a function