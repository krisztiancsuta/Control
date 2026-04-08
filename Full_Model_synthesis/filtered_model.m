close all;
%% PARAMÉTEREK DEFINIÁLÁSA
Caf = 222685.8 / 2; % Cornering stiffness az első kerékhez
Car = 136242.8 / 2; % Cornering stiffness a hátsó kerékhez
lf = 1236e-3; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 2789e-3 - lf; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 2300; % Az autó tömege
Iz = 2873; % Az autó tehetetlenségi nyomatéka
b = 50;

Nc = 4;
Np = 30;
Q = [10 0 0;...
     0 1 0;
     0 0 1];
Rw = [60 0 0;
      0 0.00001 0;...
      0 0 1];

% Diszkretizálási paraméterek
Ts = 0.033; 

vx_vector = 1 : 1 : 30; 
num_models = length(vx_vector);

%%
A0 = [0 1 0                      0 0;...
      0 0 (2*Caf+2*Car)/m        0 0;...
      0 0 0                      1 0;...
      0 0 (2*lf*Caf-2*lr*Car)/Iz 0 0;...
      0 0 0                      0 -b/m];
      
A1 = [0 0                         0 0                                 0;...
      0 -(2*Caf+2*Car)/m          0 (-2*lf*Caf+2*lr*Car)/m            0;...
      0 0                         0 0                                 0;...
      0 -(2*lf*Caf-2*lr*Car)/Iz   0 -(2*lf*lf*Caf+2*lr*lr*Car)/Iz     0;...
      0 0                         0 0                                 0];

Bc1 = [0; 2*Caf/m; 0; 2*lf*Caf/Iz; 0];
Bc2 = [0; 0; 0; 0; 1/m];

Bd1 = [0; -(2*lf*Caf-2*lr*Car)/m; 0; -(2*lf*lf*Caf+2*lr*lr*Car)/Iz; 0];
Bd2 = [0; -1; 0; 0; 0];

% Kimeneti mátrix
C = [1 0 0 0 0;...
     0 0 0 0 1];
D = zeros(size(C, 1), 3); % 3 bemenet van: Bc1, Bc2, Bd

%% MODELLEK GENERÁLÁSA
sys_c_models = cell(num_models, 1); % Folytonos modellek
sys_d_models = cell(num_models, 1); % Diszkrét modellek
Phi_Phi_models = cell(num_models, 1);
Phi_F_models   = cell(num_models, 1);
Phi_R_models   = cell(num_models, 1);
F_models       = cell(num_models, 1);
Phi_models     = cell(num_models, 1);

for i = 1:num_models
    % Aktuális sebesség a ciklusban
    vx = vx_vector(i);
    
    % Mátrixok összeállítása az aktuális sebességre
    A_current = A0 + (1/vx) * A1;
    Bd_current = (1/vx) * Bd1 + vx * Bd2;
    
    % Teljes bemeneti mátrix összerakása [Kormányzás, Erő, Zavarás]
    B_current = [Bc1, Bc2, Bd_current];
    
    % Folytonos állapottér-modell létrehozása
    sys_c_models{i} = ss(A_current, B_current, C, D);
    
    % Diszkretizálás Zero-Order Hold (ZOH) módszerrel
    sys_d_models{i} = c2d(sys_c_models{i}, Ts, 'zoh');
    
    % mpc gian szamolasa minden sebessegre.
    [Phi_Phi,Phi_F,Phi_R,A_e,B_e,~,F,Phi]=mpcgain(sys_d_models{i}.A,sys_d_models{i}.B(:,1:2),sys_d_models{i}.C,Nc,Np,Q);
    % MPC mátrixok elmentése a cellatömbökbe
    Phi_Phi_models{i} = Phi_Phi;
    Phi_F_models{i}   = Phi_F;
    Phi_R_models{i}   = Phi_R;
    F_models{i}       = F;
    Phi_models{i}     = Phi;
       
end
 fprintf('Modellek generálva.\n');

%% Receding Horizon control
t_vec = 0:Ts:35; 
N_sim = length(t_vec);
n_in = 2; % Bemenetek szama
m1 = 2; % Kimenetek szama 
n = 5; % allapotok szama

xm = [0; 0; 0; 0; 0.1]; 
y = C * xm;
u = [0; 0];
Xf = [zeros(n, 1); y]; % A kiterjesztett állapot kezdeti értéke: [dx; y]

Y = zeros(N_sim, m1);
U = zeros(N_sim, n_in);
delta_x = zeros(N_sim, n);
DU = zeros(N_sim, n_in);

%% Referencia jel definialasa
R = kron(eye(Nc), Rw);

frekvencia = 0.05;
ref_v = 10 + 2 * sign(sin(2 * pi * frekvencia * t_vec));

r = generate8likePath(t_vec,10/10);
%% Constraints
u_max = [0.5;...% Max steering angle
         10000];% Max pedal force
u_min = -u_max;
du_max = [0.02;... % Max rate of steering angle change
          2000];    % Max rate of pedal force change
du_min = -du_max; 

Umax = repmat(u_max, Nc, 1);
Umin = repmat(u_min, Nc, 1);
dUmax = repmat(du_max, Nc, 1);
dUmin = repmat(du_min, Nc, 1);

% Kimeneti korlátok definiálása (ha az út széle pl. +/- 1 méter)
y_max = [1;... % Max letaral deviation [m]
         20]; % Max speed [m/s]
y_min = [-1;... % Min letaral deviation [m]
         0.1];
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
%% SZIMULÁCIÓS CIKLUS
for kk = 1:N_sim
    v_current = y(2)
    v_idx = round(v_current);
    v_idx = max(1, min(30, v_idx)); % Ne tudjon túlcsordulni az index
    
    % Mátrixok kinyerése a cellatömbökből a v_idx alapján
    Phi_Phi_k = Phi_Phi_models{v_idx};
    Phi_F_k   = Phi_F_models{v_idx};
    Phi_R_k   = Phi_R_models{v_idx};
    F_k = F_models{v_idx};
    Phi_k = Phi_models{v_idx};

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
    ref_current = [0; ref_v(kk)];

    H = (Phi_Phi_k + R);
    f = -(Phi_R_k * ref_current - Phi_F_k * Xf);
    
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

    Ad_k   = sys_d_models{v_idx}.A;
    Bd_u_k = sys_d_models{v_idx}.B(:, 1:2); % Szabályozott bemenetek
    Bd_d_k = sys_d_models{v_idx}.B(:, 3);   % Zavaró bemenet 
  
    xm = Ad_k * xm + Bd_u_k * u + Bd_d_k *r(kk);
    y = C * xm
    
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

figure('Name', 'MPC Járműdinamika Analízis', 'Position', [100, 100, 1200, 800]);

% --- 1. ábra: Kimenetek (Y) és Referenciakövetés ---
% Oldalirányú pozíció (y)
subplot(3, 2, 1);
plot(t, Y(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(t, zeros(size(t)), 'r--', 'LineWidth', 1.5);
xlabel('Idő [s]'); ylabel('y [m]');
title('Oldalirányú pozíció követése');
legend('Valós y', 'Referencia (0)');

% Hosszirányú sebesség (v_x)
subplot(3, 2, 2);
plot(t, Y(:, 2), 'g-', 'LineWidth', 2); hold on;
plot(t, ref_v, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Idő [s]'); ylabel('v_x [m/s]');
title('Hosszirányú sebesség követése');
legend('Valós v_x', 'Referencia (10)');

% --- 2. ábra: Bemenetek (U) és Korlátok ---
% Kormányszög beavatkozás (\delta_f)
subplot(3, 2, 3);
stairs(t, U(:, 1), 'k', 'LineWidth', 1.5); hold on;
plot(t, u_max(1) * ones(size(t)), 'r:', 'LineWidth', 1.5);
plot(t, u_min(1) * ones(size(t)), 'r:', 'LineWidth', 1.5);
grid on;
xlabel('Idő [s]'); ylabel('\delta_f [rad]');
title('Kormányszög (u_1) és korlátai');
legend('u_1', 'u_{max/min}');

% Hosszirányú erő beavatkozás (F_x)
subplot(3, 2, 4);
stairs(t, U(:, 2), 'm', 'LineWidth', 1.5); hold on;
plot(t, u_max(2) * ones(size(t)), 'r:', 'LineWidth', 1.5);
plot(t, u_min(2) * ones(size(t)), 'r:', 'LineWidth', 1.5);
grid on;
xlabel('Idő [s]'); ylabel('F_x [N]');
title('Hosszirányú erő (u_2) és korlátai');
legend('u_2', 'u_{max/min}');

% --- 3. ábra: Zavarójel profilja ---
subplot(3, 2, 5);
plot(t, r(1:N_sim), 'c-', 'LineWidth', 1.5);
grid on;
xlabel('Idő [s]'); ylabel('Zavarás amplitúdó');
title('Alkalmazott zavarójel (Path profil)');
legend('r(t)');

% --- 4. ábra: Állapotváltozók növekményei (\Delta x) ---
subplot(3, 2, 6);
plot(t, delta_x(:, 1), 'b', 'LineWidth', 1); hold on;
plot(t, delta_x(:, 3), 'r', 'LineWidth', 1);
plot(t, delta_x(:, 5), 'k', 'LineWidth', 1);
grid on;
xlabel('Idő [s]'); ylabel('\Delta x');
title('Főbb állapotváltozók növekményei');
legend('\Delta y', '\Delta \psi', '\Delta v_x');

sgtitle('LPV-MPC Járműirányítási Szimuláció Eredményei');
