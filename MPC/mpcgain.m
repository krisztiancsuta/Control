function [Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e,F,Phi]=mpcgain(Ap,Bp,Cp,Nc,Np);
% This function computes the following matrices:
% Phi is a toeplitz matrix
% where Phi_Phi,Phi_F,Phi_R are for helping calculations for delta U
% A_e, B_e,C_e are the augmented states space model
% n1 is the number of states in Ap, rows of Ap
% m1 is the number of outputs, rows of Cp
% n_in the number of inputs, columns of Bp
[m1,n1]=size(Cp);
[~,n_in]=size(Bp);

% The augmented model for this plant is given by
% A_e = [Ap    0m1;...  B_e = [Bp;... C_e = [0m1 1]
%       Cp*Ap eye];           Cp*Bp]; 
A_e=eye(n1+m1,n1+m1); % Creating correct matrix dimension
A_e(1:n1,1:n1)=Ap; % Fill in Ap
A_e(n1+1:n1+m1,1:n1)=Cp*Ap; % Fill in Cp*Ap
B_e=zeros(n1+m1,n_in); % Creating correct matrix dimension
B_e(1:n1,:)=Bp; % Fill in Bp
B_e(n1+1:n1+m1,:)=Cp*Bp; %Fill in Cp*Bp
C_e=zeros(m1,n1+m1); % Creating correct matrix dimension
C_e(:,n1+1:n1+m1)=eye(m1,m1);

% Np is the lenght of the prediction horizon
% Nc is the length of control horizon
n=n1+m1; % Number of states in the augmented state space model

% Pre-allocating F and h for MIMO (dimensions scaled by m1)
F = zeros(Np*m1, n);
h = zeros(Np*m1, n);

% Initialize the first block row
h(1:m1,:)=C_e;
F(1:m1,:)=C_e*A_e;

% Iteratively calculate the free response matrices for the horizon
for kk=2:Np
    h((kk-1)*m1+1:kk*m1,:) = h((kk-2)*m1+1:(kk-1)*m1,:) * A_e;
    F((kk-1)*m1+1:kk*m1,:) = F((kk-2)*m1+1:(kk-1)*m1,:) * A_e;
end

% v represents the impulse response blocks
v=h*B_e; 

% Declare the dimension of Phi for MIMO (Np*m1 rows and Nc*n_in columns)
Phi=zeros(Np*m1, Nc*n_in); 

% Construct the Block-Toeplitz matrix Phi
for i=1:Nc
    % Each column block is shifted down by m1 rows for every control step
    Phi((i-1)*m1+1:Np*m1, (i-1)*n_in+1:i*n_in) = v(1:(Np-i+1)*m1, :);
end

% BarRs represents the reference vector over the horizon
% For MIMO, we assume a setpoint for each output
BarRs = repmat(eye(m1), Np, 1);

Phi_Phi= Phi'*Phi;
Phi_F= Phi'*F;
Phi_R=Phi'*BarRs;