%% Longitudinális Járműirányítás MPC-vel: FFWD vs. No-FFWD
clear; clc; close all;

%% 1. Fizikai paraméterek és Modell
air_density = 1.225;
Cd = 0.3;
Am = 2.2;
v0 = 25;
b_aero = air_density * Cd * Am ;
m  = 1000;
b = 50;

% Folytonos állapot-tér modell (v' = -b/m*v + 1/m*u)
A = -b/m;
B = 1/m;
C = 1; 
c_ss = ss(A, B, C, 0);

% Diszkretizálás
Ts = 33e-3; % 33 ms mintavételi idő
d_ss = c2d(c_ss, Ts);
Ap = d_ss.A;
Bp = d_ss.B;
Cp = d_ss.C;

%% 2. MPC tervezés
Nc = 50; % Kontroll horizont
Np = 100; % Predikciós horizont
rw = 0.0001; % Súlyozás a beavatkozó jel változására

% MPC mátrixok generálása (Feltételezve, hogy az mpcgain elérhető)
[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e, F, Phi] = mpcgain(Ap, Bp, Cp, Nc, Np);

%% 3. Szimulációs beállítások
t_vec = 0:Ts:35.7; 
N_sim = length(t_vec);
[m1, n1] = size(Cp);
[~, n_in] = size(B_e);

% Referencia jel és FFWD komponens
r = zeros(size(t_vec));
u_ff_all = zeros(size(t_vec));

r(t_vec >= 10) = 25; 
r(t_vec >= 20) = 5; 

% FFWD számítása a referencia sebesség alapján (stacionárius egyensúlyhoz)
u_ff_all(t_vec >= 10) = 0.5 * Cd * air_density * Am * 25^2;
u_ff_all(t_vec >= 20) = 0.5 * Cd * air_density * Am * 5^2;

% Korlátok
u_max = 10000; u_min = -10000;
du_max = 500; du_min = -500;

% Megszorítás mátrixok (QP formátumhoz)
I_n_in = eye(n_in);
C1 = repmat(I_n_in, Nc, 1);
C2 = kron(tril(ones(Nc)), I_n_in);
M1 = [-C2; C2];
M2 = [-eye(Nc*n_in); eye(Nc*n_in)];
M = [M1; M2];
H = (Phi_Phi + rw*eye(Nc*n_in));

%% 4. Futtatás: Szimulációs loop függvényben az újrahasznosíthatóságért
sim_with_ffwd = run_simulation(true, u_ff_all, N_sim, n_in, n1, m1, r, H, M, C1, Phi_R, Phi_F, Ap, Bp, Cp, u_min, u_max, du_min, du_max, Nc);
sim_no_ffwd   = run_simulation(false, u_ff_all, N_sim, n_in, n1, m1, r, H, M, C1, Phi_R, Phi_F, Ap, Bp, Cp, u_min, u_max, du_min, du_max, Nc);

%% 5. Ábrázolás
figure('Color', 'w', 'Position', [100, 100, 900, 600]);

% Sebesség összehasonlítása
subplot(2,1,1);
plot(t_vec, r, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Referencia'); hold on;
plot(t_vec, sim_with_ffwd.Y, 'b', 'LineWidth', 2, 'DisplayName', 'MPC + FFWD');
plot(t_vec, sim_no_ffwd.Y, 'r', 'LineWidth', 1.5, 'DisplayName', 'Csak MPC');
grid on; ylabel('Sebesség [m/s]');
title('MPC alapú sebességkövetés összehasonlítása');
legend('Location', 'best');

% Beavatkozó jel (Pedál) összehasonlítása
subplot(2,1,2);
plot(t_vec, sim_with_ffwd.U / 10000, 'b', 'LineWidth', 2, 'DisplayName', 'MPC + FFWD'); hold on;
plot(t_vec, sim_no_ffwd.U / 10000, 'r', 'LineWidth', 1.5, 'DisplayName', 'Csak MPC');
grid on; ylabel('Pedálpozíció (normált)');
xlabel('Idő [s]');
title('Beavatkozó jelek');
legend('Location', 'best');


%% --- Módosított Segédfüggvény ---
function res = run_simulation(use_ffwd, phys, N_sim, n_in, n1, m1, r, H, M, C1, Phi_R, Phi_F, Ap, Bp, Cp, u_min, u_max, du_min, du_max, Nc)
    xm = 0;      % Aktuális sebesség (állapot)
    u_mpc = 0;   % Az MPC által számolt kumulált beavatkozó jel
    Xf = zeros(n1 + m1, 1); % Kiterjesztett állapotvektor [delta_x; y]
    
    U_total_out = zeros(N_sim, 1);
    Y_out = zeros(N_sim, 1);
    
    for kk = 1:N_sim
        % 1. QP KORLÁTOK FRISSÍTÉSE
        % A korlátokat az MPC által vezérelt u_mpc-re vetítjük
        gamma = [- (ones(Nc,1)*u_min) + C1*u_mpc; (ones(Nc,1)*u_max) - C1*u_mpc; -ones(Nc,1)*du_min; ones(Nc,1)*du_max];
        
        % 2. QP MEGOLDÁS (Inkrementális delta_u)
        f = -(Phi_R * r(kk) - Phi_F * Xf);
        deltaU = QPhild(H, f, M, gamma);
        delta_u_k = deltaU(1:n_in);
        
        % 3. MPC JEL FRISSÍTÉSE
        u_mpc = u_mpc + delta_u_k;
        
        % 4. DINAMIKUS ELŐRECSATOLÁS (FFD) SZÁMÍTÁSA
        % Itt használjuk az aktuális sebességet (xm)
        u_ff = 0;
        if use_ffwd
            v_act = xm; 
            u_ff = 0.5 * 0.8085 * v_act^2;
        end
        
        % 5. TELJES BEAVATKOZÓ JEL (MPC + FFD)
        u_total = u_mpc + u_ff;
        u_total = max(min(u_total, u_max), u_min); % Fizikai korlátozás
        
        % 6. RENDSZER (PLANT) SZIMULÁCIÓ
        xm_old = xm;
        % A plant a teljes erőt (u_total) kapja meg
        xm = Ap * xm + Bp * u_total; 
        y = Cp * xm;
        
        % 7. MPC ÁLLAPOT-UPDATE (Fontos a visszacsatoláshoz)
        % Az MPC-nek tudnia kell, mi történt a kimenettel
        Xf = [(xm - xm_old); y];
        
        % Mentés
        Y_out(kk) = y;
        U_total_out(kk) = u_total;
    end
    res.Y = Y_out;
    res.U = U_total_out;
end