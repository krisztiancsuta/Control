close all;
%% PARAMÉTEREK DEFINIÁLÁSA
Caf = 222685.8 / 2; % Cornering stiffness az első kerékhez
Car = 136242.8 / 2; % Cornering stiffness a hátsó kerékhez
lf = 1236e-3; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 2789e-3 - lf; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 2300; % Az autó tömege
Iz = 2873; % Az autó tehetetlenségi nyomatéka
b = 50;

Nc = 20;
Np = 30;
Q = [10 0;...
     0 1];
Rw = [1 0;
      0 0.00001];

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

xm = [0; 0; 0; 0; 5]; 
y = C * xm;
u = [0; 0];
Xf = [zeros(n, 1); y]; % A kiterjesztett állapot kezdeti értéke: [dx; y]

Y = zeros(N_sim, m1);
U = zeros(N_sim, n_in);
delta_x = zeros(N_sim, n);
DU = zeros(N_sim, n_in);

%% Referencia jel definialasa
R = kron(eye(Nc), Rw);
% Referencia jel
u_max = [ 0.5;  10000]; % Kormányszög max [rad], Erő max [N]
u_min = [-0.5; -10000]; % Kormányszög min [rad], Erő min [N]

frekvencia = 0.02;
ref_v = 10 + 3 * sign(sin(2 * pi * frekvencia * t_vec));

r = generate8likePath(t_vec,1);
%% SZIMULÁCIÓS CIKLUS
for kk = 1:N_sim
    % --- I. AZ AKTUÁLIS LPV MODELL KIVÁLASZTÁSA ---
    v_current = xm(5); % Aktuális hosszirányú sebesség
    
    % Sebesség indexszé alakítása (1 és 30 közé szorítva)
    v_idx = round(v_current);
    v_idx = max(1, min(30, v_idx)); % Ne tudjon túlcsordulni az index
    
    % Mátrixok kinyerése a cellatömbökből a v_idx alapján
    Phi_Phi_k = Phi_Phi_models{v_idx};
    Phi_F_k   = Phi_F_models{v_idx};
    Phi_R_k   = Phi_R_models{v_idx};
    
    Ad_k   = sys_d_models{v_idx}.A;
    Bd_u_k = sys_d_models{v_idx}.B(:, 1:2); % Szabályozott bemenetek
    Bd_d_k = sys_d_models{v_idx}.B(:, 3);   % Zavaró bemenet
    
    ref_current = [0;ref_v(kk)];
    deltaU = (Phi_Phi_k + R) \ (Phi_R_k * ref_current - Phi_F_k * Xf);
    
    % Csak a legelső lépés levágása (Receding Horizon: 2 elem kell)
    delta_u_k = deltaU(1:n_in);
    
    % --- III. BEAVATKOZÓJEL FRISSÍTÉSE ÉS KORLÁTOZÁSA ---
    u_old = u;
    u = u_old + delta_u_k;
    
    u = max(min(u, u_max), u_min);
    
    delta_u_k_applied = u - u_old; 
    
    xm_old = xm;

    xm = Ad_k * xm + Bd_u_k * u + Bd_d_k * r(kk);
    y = C * xm;
    
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

%% Ábrázolás (Plotting Results)
t = t_vec(1:N_sim); % Idővektor a plotokhoz

figure('Name', 'MPC Járműdinamika Analízis', 'Position', [100, 100, 1200, 800]);

% --- 1. ábra: Kimenetek (Y) és Referenciakövetés ---
% Oldalirányú pozíció (y)
subplot(3, 2, 1);
plot(t, Y(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(t, 0 , 'r--', 'LineWidth', 1.5);
grid on;
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
