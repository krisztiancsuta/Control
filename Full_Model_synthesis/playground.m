close all;
clear; clc;

%% PARAMÉTEREK DEFINIÁLÁSA
Caf = 222685.8 / 2; 
Car = 136242.8 / 2; 
lf = 1236e-3; 
lr = 2789e-3 - lf; 
m = 2300; 
Iz = 2873; 
b = 50;
Nc = 10;
Np = 30;
% Súlyozás: [y, vx, filtered_ax, filtered_ay]
Q = diag([10, 1, 5, 5]); 
Rw = diag([1, 0.00001]);

% Diszkretizálási paraméterek
Ts = 0.033; 
vx_vector = 1 : 0.1 : 30; 
num_models = length(vx_vector);

%% Alap mátrixok (Vehicle Dynamics)
Av_const = [0 1 0                      0 0;...
      0 0 (2*Caf+2*Car)/m        0 0;...
      0 0 0                      1 0;...
      0 0 (2*lf*Caf-2*lr*Car)/Iz 0 0;...
      0 0 0                      0 -b/m];
      
Av_var = [0 0                         0 0                             0;...
      0 -(2*Caf+2*Car)/m          0 (-2*lf*Caf+2*lr*Car)/m            0;...
      0 0                         0 0                                 0;...
      0 -(2*lf*Caf-2*lr*Car)/Iz   0 -(2*lf*lf*Caf+2*lr*lr*Car)/Iz     0;...
      0 0                         0 0                                 0];

Bv1 = [0; 2*Caf/m; 0; 2*lf*Caf/Iz; 0]; % Steering
Bv2 = [0; 0; 0; 0; 1/m];               % Pedal force
Bd1 = [0; -(2*lf*Caf-2*lr*Car)/m; 0; -(2*lf*lf*Caf+2*lr*lr*Car)/Iz; 0];
Bd2 = [0; -1; 0; 0; 0];

Cv = [1 0 0 0 0;...
      0 0 0 0 1];

%% Filter transfer function (ISO 2631)
num_wf = [0.02633, 0.0238, 2.286, 0.2335, 0.02902];
den_wf = [1, 2.527, 4.584, 2.993, 1.373];
[Af, Bf, Cf, Df] = tf2ss(num_wf, den_wf);

%% MODELLEK GENERÁLÁSA (LPV + Dual Filter)
sys_c_models = cell(num_models, 1);
sys_d_models = cell(num_models, 1);
Phi_Phi_models = cell(num_models, 1);
Phi_F_models   = cell(num_models, 1);
Phi_R_models   = cell(num_models, 1);
F_models       = cell(num_models, 1);
Phi_models     = cell(num_models, 1);

for i = 1:num_models
    vx = vx_vector(i);
    Av_current = Av_const + (1/vx) * Av_var;
    Bd_current = (1/vx) * Bd1 + vx * Bd2;
    Bv_current = [Bv1, Bv2, Bd_current];
    
    C_ax = Av_current(5, :); D_ax = Bv_current(5, :);
    C_ay = Av_current(2, :); D_ay = Bv_current(2, :);

    Ac_current = [ Av_current   zeros(5,4)   zeros(5,4);...
                   Bf*C_ax      Af           zeros(4,4);...
                   Bf*C_ay      zeros(4,4)   Af ];
    Bc_current = [ Bv_current; Bf*D_ax; Bf*D_ay ];
    Cc_current = [ Cv           zeros(2,4)   zeros(2,4);...
                   Df*C_ax      Cf           zeros(1,4);...
                   Df*C_ay      zeros(1,4)   Cf ];
    Dc_current = [ zeros(2,3); Df*D_ax; Df*D_ay ];

    sys_c_models{i} = ss(Ac_current, Bc_current, Cc_current, Dc_current);
    sys_d_models{i} = c2d(sys_c_models{i}, Ts, 'zoh');
    
    [Phi_Phi,Phi_F,Phi_R,~,~,~,F,Phi] = mpcgain(sys_d_models{i}.A, sys_d_models{i}.B(:,1:2), sys_d_models{i}.C, Nc, Np, Q);
    
    Phi_Phi_models{i} = Phi_Phi;
    Phi_F_models{i}   = Phi_F;
    Phi_R_models{i}   = Phi_R;
    F_models{i}       = F;
    Phi_models{i}     = Phi;
end

%% SZIMULÁCIÓ INICIALIZÁLÁSA
t_vec = 0:Ts:35; 
N_sim = length(t_vec);
n_in = 2; m_out = 4; n_total = 13; 

xm = zeros(n_total, 1); xm(5) = 0.1; % 1 m/s kezdeti sebesség
u = [0; 0];
Xf = [zeros(n_total, 1); zeros(m_out, 1)]; 

% Tárolók inicializálása - EZ A RÉSZ KRITIKUS
Y = zeros(N_sim, m_out);
U = zeros(N_sim, n_in);
delta_x = zeros(N_sim, n_total); % <--- Itt jön létre a delta_x
a_raw_x = zeros(N_sim, 1);
a_raw_y = zeros(N_sim, 1);

% Korlátok és referencia
R_mat = kron(eye(Nc), Rw);
ref_v = 8 * ones(1, N_sim);
r_dist = generate8likePath(t_vec, 8/30);

u_max = [0.5; 10000]; u_min = -u_max;
du_max = [0.1; 2000]; du_min = -du_max;
Ymax = repmat([1; 30; 20; 20], Np, 1); 
Ymin = repmat([-1; 0; -20; -20], Np, 1);

% MPC Mátrixok a korlátokhoz
C1 = repmat(eye(n_in), Nc, 1);
C2 = kron(tril(ones(Nc)), eye(n_in));
M1 = [-C2; C2]; M2 = [-eye(Nc*n_in); eye(Nc*n_in)];
lambda0 = zeros(2*(Nc*n_in) + 2*(Nc*n_in) + 2*(Np*m_out), 1);

%% SZIMULÁCIÓS CIKLUS
for kk = 1:N_sim
    v_idx = max(1, min(num_models, round((xm(5) - 1) / 0.1) + 1));
    
    Phi_Phi_k = Phi_Phi_models{v_idx};
    Phi_F_k = Phi_F_models{v_idx};
    Phi_R_k = Phi_R_models{v_idx};
    Phi_k = Phi_models{v_idx};
    F_k = F_models{v_idx};
    Cc_k = sys_d_models{v_idx}.C;

    M = [M1; M2; [-Phi_k; Phi_k]];
    gamma = [-repmat(u_min,Nc,1) + C1*u; repmat(u_max,Nc,1) - C1*u; ...
             -repmat(du_min,Nc,1); repmat(du_max,Nc,1); ...
             -Ymin + F_k*Xf; Ymax - F_k*Xf];
    
    ref_current = [0; ref_v(kk); 0; 0]; 
    H = Phi_Phi_k + R_mat;
    f = -(Phi_R_k * ref_current - Phi_F_k * Xf);
    
    [deltaU, ~, lambda0, ~] = QPhild(H, f, M, gamma, lambda0);
    delta_u = deltaU(1:n_in);
    
    % Mentés előtt mentsük el az állapotot
    xm_old = xm;
    
    % Update
    u = u + delta_u;
    xm = sys_d_models{v_idx}.A * xm + sys_d_models{v_idx}.B(:,1:2) * u + sys_d_models{v_idx}.B(:,3) * r_dist(kk);
    y_vec = Cc_k * xm;
    
    % Nyers gyorsulások
    a_raw_x(kk) = sys_c_models{v_idx}.A(5,1:5) * xm(1:5) + sys_c_models{v_idx}.B(5,1:2) * u;
    a_raw_y(kk) = sys_c_models{v_idx}.A(2,1:5) * xm(1:5) + sys_c_models{v_idx}.B(2,1:3) * [u; r_dist(kk)];
    
    % ADATOK ELMENTÉSE - EZT HIÁNYOLTA A MATLAB
    delta_x(kk, :) = (xm - xm_old)'; 
    Y(kk, :) = y_vec';
    U(kk, :) = u';
    
    % Xf frissítése
    Xf = [(xm - xm_old); y_vec];
end

%% VIZUALIZÁCIÓ
t = t_vec(1:N_sim);

% 1. Ábra: Gyorsulás
figure('Name', 'Gyorsulás és Kényelmi Analízis', 'Position', [100, 100, 1000, 700]);
subplot(2, 1, 1);
plot(t, a_raw_x, 'Color', [0.8 0.8 0.8]); hold on;
plot(t, Y(:, 3), 'r', 'LineWidth', 2); title('Hosszirányú gyorsulás (ax)'); grid on;
subplot(2, 1, 2);
plot(t, a_raw_y, 'Color', [0.8 0.8 0.8]); hold on;
plot(t, Y(:, 4), 'b', 'LineWidth', 2); title('Oldalirányú gyorsulás (ay)'); grid on;

% 2. Ábra: Járműdinamika
figure('Name', 'MPC Járműdinamika Analízis', 'Position', [150, 150, 1200, 900]);
subplot(3, 2, 1); plot(t, Y(:, 1), 'b'); title('Laterális eltérés (y)'); grid on;
subplot(3, 2, 2); plot(t, Y(:, 2), 'g'); hold on; plot(t, ref_v, 'r--'); title('Sebesség (vx)'); grid on;
subplot(3, 2, 3); stairs(t, U(:, 1), 'k'); title('Kormányszög'); grid on;
subplot(3, 2, 4); stairs(t, U(:, 2), 'm'); title('Pedálerő'); grid on;
subplot(3, 2, 5); plot(t, r_dist, 'c'); title('Út görbülete'); grid on;
subplot(3, 2, 6); 
plot(t, delta_x(:, 1), 'b', t, delta_x(:, 3), 'r', t, delta_x(:, 5), 'k');
title('Állapotváltozók változása (\Delta x)'); grid on;
legend('\Delta y', '\Delta \psi', '\Delta v_x');

if exist('drawpath', 'file')
    drawpath(t_vec(1:N_sim), Y(:, 1), ref_v(1), r_dist);
end