close all;
%% PARAMÉTEREK DEFINIÁLÁSA
Caf = 222685.8 / 2; % Cornering stiffness az első kerékhez
Car = 136242.8 / 2; % Cornering stiffness a hátsó kerékhez
lf = 1236e-3; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 2789e-3 - lf; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 2300; % Az autó tömege
Iz = 2873; % Az autó tehetetlenségi nyomatéka
b = 50;

Nc = 10;
Np = 30;
Q = diag([10,1,40]);
Rw = diag([1, 0.00001]);

% Diszkretizálási paraméterek
Ts = 0.033; 

%% 2D MUNKAPONTI RÁCS (Grid) DEFINÍCIÓ
vx_vector = 1 : 0.5 : 30;          % Sebesség tengely [m/s]
psi_des_dot_vector = -1.5 : 0.1 : 1.5; % Kanyarodási sebesség tengely [rad/s]
num_vx  = length(vx_vector);
num_psi = length(psi_des_dot_vector);

%% Model 
% Av_const = constant matrix
% Av_var = variable matrix
Av_const =   [0 1 0                      0 0;...
              0 0 (2*Caf+2*Car)/m        0 0;...
              0 0 0                      1 0;...
              0 0 (2*lf*Caf-2*lr*Car)/Iz 0 0;...
              0 0 0                      0 -b/m];
      
Av_var = [0 0                         0 0                             0;...
      0 -(2*Caf+2*Car)/m          0 (-2*lf*Caf+2*lr*Car)/m            0;...
      0 0                         0 0                                 0;...
      0 -(2*lf*Caf-2*lr*Car)/Iz   0 -(2*lf*lf*Caf+2*lr*lr*Car)/Iz     0;...
      0 0                         0 0                                 0];

a1 = (2*Caf+2*Car)/m;
a2 = (2*lf*Caf-2*lr*Car)/m;
b1 = (2*lf*Caf-2*lr*Car)/Iz;
b2 = (2*lf*lf*Caf+2*lr*lr*Car)/Iz;


Bv1 = [0; 2*Caf/m; 0; 2*lf*Caf/Iz; 0]; % Steering
Bv2 = [0; 0; 0; 0; 1/m];% Pedal force

Bd1 = [0; -(2*lf*Caf-2*lr*Car)/m; 0; -(2*lf*lf*Caf+2*lr*lr*Car)/Iz; 0];
Bd2 = [0; -1; 0; 0; 0];

% Kimeneti mátrix
Cv = [1 0 0 0 0;...
      0 0 0 0 1];
D = zeros(size(Cv, 1), 3); % 3 bemenet van: Bc1, Bc2, Bd

%% Filter transfer function

num_wf = [0.02633, 0.0238, 2.286, 0.2335, 0.02902];
den_wf = [1, 2.527, 4.584, 2.993, 1.373];
[Af, Bf, Cf, Df] = tf2ss(num_wf, den_wf);

%% OFF-LINE MÁTRIX-GENERÁLÁS (2D Pre-calculation)
sys_c_models   = cell(num_vx, num_psi);
sys_d_models   = cell(num_vx, num_psi);
Phi_Phi_models = cell(num_vx, num_psi);
Phi_F_models   = cell(num_vx, num_psi);
Phi_R_models   = cell(num_vx, num_psi);
F_models       = cell(num_vx, num_psi);
Phi_models     = cell(num_vx, num_psi);

e1_dot_0 = 0.0;
e2_dot_0 = 0.0;


fprintf('2D rács: %d vx × %d psi_dot = %d modell generálása...\n', ...
        num_vx, num_psi, num_vx * num_psi);

for i = 1:num_vx
    vx = vx_vector(i);

    Av_current_base = Av_const + (1/vx) * Av_var;
    Bd_current = (1/vx) * Bd1 + vx * Bd2;
    Da1 = Bv1(2, :);
    Da2 = Bd_current(2, :);

    for j = 1:num_psi
        psi_des_dot_0 = psi_des_dot_vector(j);

        Av_current = Av_current_base;
        A25 = (a1/vx^2)*e1_dot_0 + (a2/vx^2)*e2_dot_0 + (a2/vx^2 - 1)*psi_des_dot_0;
        A45 = (b1/vx^2)*e1_dot_0 + (b2/vx^2)*e2_dot_0 + (b2/vx^2)*psi_des_dot_0;
        Av_current(2,5) = A25;
        Av_current(4,5) = A45;

        Bv_current = [Bv1, Bv2, Bd_current];
        Ca = Av_current(2, 1:5);

        Ac_current = [Av_current zeros(5,4); Bf*Ca Af];
        Bc_current = [Bv_current; Bf*Da1 zeros(4,1) Bf*Da2];
        Cc_current = [Cv zeros(2,4); Df*Ca Cf];
        Dc_current = [zeros(2,3); Df*Da1 0 Df*Da2];

        sys_c_models{i,j} = ss(Ac_current, Bc_current, Cc_current, Dc_current);
        sys_d_models{i,j} = c2d(sys_c_models{i,j}, Ts, 'zoh');

        [Phi_Phi,Phi_F,Phi_R,~,~,~,F,Phi] = mpcgain( ...
            sys_d_models{i,j}.A, sys_d_models{i,j}.B(:,1:2), ...
            sys_d_models{i,j}.C, Nc, Np, Q);

        Phi_Phi_models{i,j} = Phi_Phi;
        Phi_F_models{i,j}   = Phi_F;
        Phi_R_models{i,j}   = Phi_R;
        F_models{i,j}       = F;
        Phi_models{i,j}     = Phi;
    end

    if mod(i, 10) == 0
        fprintf('  ... %d / %d vx kész\n', i, num_vx);
    end
end
fprintf('Modellek generálva: %d × %d = %d munkapont.\n', num_vx, num_psi, num_vx*num_psi);

%% Receding Horizon control
t_vec = 0:Ts:35; 
N_sim = length(t_vec);
n_in = 2; % Bemenetek szama
m1 = 3; % Kimenetek szama [y; vx; filtered_accel]
n = 9; % allapotok szama 5 vehicle + 4 filter states

xm = zeros(n, 1); 
xm(5) = 0.1;
y = zeros(m1,1);% Outputs: [y; vx; filtered_accel]
u = [0; 0];
Xf = [zeros(n, 1); y]; % A kiterjesztett állapot kezdeti értéke: [dx; y]

Y = zeros(N_sim, m1);
U = zeros(N_sim, n_in);
delta_x = zeros(N_sim, n);
DU = zeros(N_sim, n_in);

%% Referencia jel definialasa
R = kron(eye(Nc), Rw);

frekvencia = 0.05;
ref_v = 10 + 5 * sign(sin(2 * pi * frekvencia * t_vec));

r = generate8likePath(t_vec,10/10);
%% Constraints
u_max = [0.5;...% Max steering angle
         10000];% Max pedal force
u_min = -u_max;
du_max = [0.08;... % Max rate of steering angle change
          2000];    % Max rate of pedal force change
du_min = -du_max; 

Umax = repmat(u_max, Nc, 1);
Umin = repmat(u_min, Nc, 1);
dUmax = repmat(du_max, Nc, 1);
dUmin = repmat(du_min, Nc, 1);

% Kimeneti korlátok definiálása (ha az út széle pl. +/- 1 méter)
y_max = [1;... % Max letaral deviation [m]
         1000;...
         1000]; % Max speed [m/s]
y_min = [-1;... % Min letaral deviation [m]
         -1000;...
         -1000];
Ymax = repmat(y_max, Np, 1);
Ymin = repmat(y_min, Np, 1);

I_n_in = eye(n_in);
C1 = repmat(I_n_in, Nc, 1);
C2 = kron(tril(ones(Nc)), I_n_in);

% The 3 contstraints enaquility can be described as:
M1 = [-C2;C2];
M2 = [-eye(Nc*n_in); eye(Nc*n_in)];

n_u_constraints = 2 * (Nc * n_in);    % Umin és Umax (2 bemenet * Nc időlépés * 2 irány)
n_du_constraints = 2 * (Nc * n_in);   % dUmin és dUmax
n_y_constraints = 2 * (Np * m1);      % Ymin és Ymax (2 kimenet * Np időlépés * 2 irány)

total_constraints = n_u_constraints + n_du_constraints + n_y_constraints;
lambda0 = zeros(total_constraints, 1);

a_raw = zeros(N_sim, 1);
grid_idx_log = zeros(N_sim, 2); % [v_idx, psi_idx] trajectory on the 2D grid

%% SZIMULÁCIÓS CIKLUS (On-line Gain Scheduling)
for kk = 1:N_sim

    v_current = xm(5);
    psi_des_current = r(kk);

    v_step   = vx_vector(2) - vx_vector(1);
    psi_step = psi_des_dot_vector(2) - psi_des_dot_vector(1);

    v_idx   = round((v_current - vx_vector(1)) / v_step) + 1;
    v_idx   = max(1, min(num_vx, v_idx));

    psi_idx = round((psi_des_current - psi_des_dot_vector(1)) / psi_step) + 1;
    psi_idx = max(1, min(num_psi, psi_idx));

    grid_idx_log(kk, :) = [v_idx, psi_idx];

    Phi_Phi_k    = Phi_Phi_models{v_idx, psi_idx};
    Cc_current_k = sys_d_models{v_idx, psi_idx}.C;
    Phi_F_k      = Phi_F_models{v_idx, psi_idx};
    Phi_R_k      = Phi_R_models{v_idx, psi_idx};
    F_k          = F_models{v_idx, psi_idx};
    Phi_k        = Phi_models{v_idx, psi_idx};

    % Constraints
    M3 = [-Phi_k;Phi_k];
  
    M = [M1; M2;M3];

    N1 = [-Umin + C1*u;... 
           Umax - C1*u];
    N2 = [-dUmin;dUmax];
    N3 = [-Ymin + F_k*Xf;... 
           Ymax - F_k*Xf];

    
    gamma = [N1;N2;N3];
    
    % Az aktuális referencia vektor frissítése
    
    ref_current = [0;ref_v(kk);0];
    H = (Phi_Phi_k + R);
    f = -(Phi_R_k *ref_current - Phi_F_k * Xf);
    
    [deltaU, iterations, lambda0, ~] = QPhild(H, f, M, gamma, lambda0);

    %deltaU = (Phi_Phi_k + R) \ (Phi_R_k * ref_current - Phi_F_k * Xf);
    
    % Csak a legelső lépés levágása (Receding Horizon: 2 elem kell)
    delta_u_k = deltaU(1:n_in);
    
    % --- III. BEAVATKOZÓJEL FRISSÍTÉSE ÉS KORLÁTOZÁSA ---
    u_old = u;
    u = u_old + delta_u_k;
    
    %u = max(min(u, u_max), u_min);
    
    delta_u_k_applied = u - u_old; 
    
    xm_old = xm;

    Ad_k   = sys_d_models{v_idx, psi_idx}.A;
    Bd_u_k = sys_d_models{v_idx, psi_idx}.B(:, 1:2);
    Bd_d_k = sys_d_models{v_idx, psi_idx}.B(:, 3);
  
    xm = Ad_k * xm + Bd_u_k * u + Bd_d_k *r(kk);
    y = Cc_current_k * xm;
    
    a_raw(kk) = sys_c_models{v_idx, psi_idx}.A(2,:) * xm ...
              + sys_c_models{v_idx, psi_idx}.B(2,1:2) * u ...
              + sys_c_models{v_idx, psi_idx}.B(2,3) * r(kk);

    
    % --- V. KITERJESZTETT ÁLLAPOT FRISSÍTÉSE A KÖVETKEZŐ LÉPÉSHEZ ---
    % Xf = [ \Delta x ; y ]
    Xf = [(xm - xm_old); y];
    
    % --- VI. ADATOK MENTÉSE ---
    % Mivel az xm, y, u oszlopvektorok, transzponáljuk (') őket sorként
    Y(kk, :) = y';
    U(kk, :) = u';
    delta_x(kk, :) = (xm - xm_old)';
    DU(kk, :) = delta_u_k_applied';
end

fprintf('Szimuláció sikeresen véget ért.\n');

%% Útvonal vizualizáció: referencia (integrált yaw-rate profil) vs. tényleges pálya
% drawpath: X_des,Y_des a r(t) yaw-rate és Vx alapján; kék görbe a laterális eltérés (Y(:,1)) alapján
drawpath(t_vec(1:N_sim), Y(:, 1), ref_v(1), r(1:N_sim));

%% Ábrázolás (Plotting Results)
t = t_vec(1:N_sim); % Idővektor a plotokhoz

% --- ÚJ ÁBRA: GYORSULÁS ANALÍZIS ---
figure('Name', 'Gyorsulás és Kényelmi Analízis', 'Position', [200, 200, 800, 500]);
plot(t, a_raw, 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'Nyers gyorsulás (Raw)');
hold on;
plot(t, Y(:, 3), 'r', 'LineWidth', 2, 'DisplayName', 'Szűrt gyorsulás (ISO/Filtered)');
grid on;
xlabel('Idő [s]');
ylabel('Gyorsulás [m/s^2]');
title('Laterális gyorsulás: Fizikai vs. Szűrt jel');
legend('Location', 'northeast');

% --- EREDETI ÁBRÁK FRISSÍTÉSE ---
figure('Name', 'MPC Járműdinamika Analízis', 'Position', [100, 100, 1200, 900]);

% Oldalirányú pozíció
subplot(3, 2, 1);
plot(t, Y(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(t, zeros(size(t)), 'r--', 'LineWidth', 1.5);
title('Laterális eltérés (y)'); grid on;

% Sebesség
subplot(3, 2, 2);
plot(t, Y(:, 2), 'g-', 'LineWidth', 2); hold on;
plot(t, ref_v, 'r--', 'LineWidth', 1.5);
title('Hosszirányú sebesség (v_x)'); grid on;

% Kormányszög
subplot(3, 2, 3);
stairs(t, U(:, 1), 'k', 'LineWidth', 1.5); hold on;
yline(u_max(1), 'r:', 'LineWidth', 1.5);
yline(u_min(1), 'r:', 'LineWidth', 1.5);
title('Kormányszög (\delta_f)'); grid on;

% Pedálerő
subplot(3, 2, 4);
stairs(t, U(:, 2), 'm', 'LineWidth', 1.5); hold on;
yline(u_max(2), 'r:', 'LineWidth', 1.5);
yline(u_min(2), 'r:', 'LineWidth', 1.5);
title('Hajtó/Fékező erő (F_x)'); grid on;

% Zavarójel (Út görbülete/Yaw-rate ref)
subplot(3, 2, 5);
plot(t, r(1:N_sim), 'c-', 'LineWidth', 1.5);
title('Zavaró jel (r) - Útprofil'); grid on;

% Állapotváltozók növekményei
subplot(3, 2, 6);
plot(t, delta_x(:, 1), 'b', t, delta_x(:, 3), 'r', t, delta_x(:, 5), 'k');
title('Állapotváltozók változása (\Delta x)'); grid on;
legend('\Delta y', '\Delta \psi', '\Delta v_x');

sgtitle('LPV-MPC Járműirányítási Szimuláció Eredményei');

%% 2D MUNKAPONTI RÁCS VIZUALIZÁCIÓ
v_traj   = vx_vector(grid_idx_log(:,1));
psi_traj = psi_des_dot_vector(grid_idx_log(:,2));

figure('Name', '2D Gain-Scheduling Munkaponti Palya', 'Position', [300, 100, 950, 750]);

subplot(2,2,[1 3]);
visit_count = zeros(num_psi, num_vx);
for kk = 1:N_sim
    vi = grid_idx_log(kk,1);
    pi = grid_idx_log(kk,2);
    visit_count(pi, vi) = visit_count(pi, vi) + 1;
end
imagesc(vx_vector, psi_des_dot_vector, log10(visit_count + 1));
set(gca, 'YDir', 'normal');
colormap(gca, parula);
cb = colorbar; cb.Label.String = 'log_{10}(latogatasok + 1)';
hold on;
plot(v_traj, psi_traj, 'r-', 'LineWidth', 1.2);
xlabel('v_x [m/s]'); ylabel('d\psi_{des}/dt [rad/s]');
title('Munkapont-latogatottsag es palya (piros)');

subplot(2,2,2);
plot(t, v_traj, 'b', 'LineWidth', 1.5);
xlabel('Ido [s]'); ylabel('v_x [m/s]');
title('Aktualis sebesseg munkapont');
grid on;

subplot(2,2,4);
plot(t, psi_traj, 'r', 'LineWidth', 1.5);
xlabel('Ido [s]'); ylabel('d\psi_{des}/dt [rad/s]');
title('Aktualis kanyarodasi sebesseg munkapont');
grid on;
