function [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e, F, Phi] = mpcgain(Ap, Bp, Cp, Nc, Np, Q)
% Bemenetek:
% Q: Kimeneti hibát súlyozó mátrix (m1 x m1)
% R: Beavatkozás változását (\Delta U) súlyozó mátrix (n_in x n_in)

[m1, n1] = size(Cp);
[~, n_in] = size(Bp);

% --- Augmentált állapottér-modell (Delta U formuláció) ---
A_e = eye(n1+m1, n1+m1); 
A_e(1:n1, 1:n1) = Ap; 
A_e(n1+1:n1+m1, 1:n1) = Cp*Ap; 

B_e = zeros(n1+m1, n_in); 
B_e(1:n1, :) = Bp; 
B_e(n1+1:n1+m1, :) = Cp*Bp; 

C_e = zeros(m1, n1+m1); 
C_e(:, n1+1:n1+m1) = eye(m1, m1);

% --- Predikciós mátrixok (F és h) számolása ---
n = n1 + m1; 
F = zeros(Np*m1, n);
h = zeros(Np*m1, n);

h(1:m1, :) = C_e;
F(1:m1, :) = C_e * A_e;

for kk = 2:Np
    h((kk-1)*m1+1:kk*m1, :) = h((kk-2)*m1+1:(kk-1)*m1, :) * A_e;
    F((kk-1)*m1+1:kk*m1, :) = F((kk-2)*m1+1:(kk-1)*m1, :) * A_e;
end

% --- Toeplitz mátrix (Phi) építése ---
v = h * B_e; 
Phi = zeros(Np*m1, Nc*n_in); 

for i = 1:Nc
    Phi((i-1)*m1+1:Np*m1, (i-1)*n_in+1:i*n_in) = v(1:(Np-i+1)*m1, :);
end

% --- Referencia követés és Költségfüggvény mátrixok ---
BarRs = repmat(eye(m1), Np, 1);

% Kronecker-szorzat a horizont menti súlyozáshoz
Q_extended = kron(eye(Np), Q);

% Phi_Phi immár tartalmazza az R mátrixot is, így ez a teljes Hesse-mátrix (H)
Phi_Phi = Phi' * Q_extended * Phi; 
Phi_F   = Phi' * Q_extended * F;
Phi_R   = Phi' * Q_extended * BarRs;

end