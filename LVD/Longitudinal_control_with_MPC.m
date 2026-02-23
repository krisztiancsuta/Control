%% Longitudinális Járműirányítás MPC-vel: FFWD vs. No-FFWD
clear; clc; close all;

%% 1. Fizikai paraméterek és Modell
air_density = 1.225;
Cd = 0.3;
Am = 2.2;
v0 = 25;
b_aero = air_density * Cd * Am * v0;
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
Nc = 30; % Kontroll horizont
Np = 300; % Predikciós horizont
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
du_max = 2000; du_min = -2000;

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

%% --- Segédfüggvény a szimulációhoz ---
function res = run_simulation(use_ffwd, u_ff_vec, N_sim, n_in, n1, m1, r, H, M, C1, Phi_R, Phi_F, Ap, Bp, Cp, u_min, u_max, du_min, du_max, Nc)
    xm = 0; u = 0;
    Xf = zeros(n1 + m1, 1);
    
    U_out = zeros(N_sim, 1);
    Y_out = zeros(N_sim, 1);
    
    Umax = ones(Nc, 1) * u_max;
    Umin = ones(Nc, 1) * u_min;
    dUmax = ones(Nc, 1) * du_max;
    dUmin = ones(Nc, 1) * du_min;

    for kk = 1:N_sim
        % Korlátok frissítése az aktuális u(k-1) alapján
        gamma = [-Umin + C1*u; Umax - C1*u; -dUmin; dUmax];
        
        % QP célfüggvény gradiens
        f = -(Phi_R * r(kk) - Phi_F * Xf);

        % QP megoldás (QPhild függvény hívása)
        deltaU = QPhild(H, f, M, gamma);
        delta_u_k = deltaU(1:n_in);
        
        % Beavatkozó jel kiszámítása
        ff_term = 0;
        if use_ffwd, ff_term = u_ff_vec(kk); end
        
        u = u + delta_u_k + ff_term;
        u = max(min(u, u_max), u_min); % Szaturáció
        
        % Plant szimuláció
        xm_old = xm;
        xm = Ap * xm + Bp * u;
        y = Cp * xm;
        
        % Állapot frissítés a következő lépéshez
        Xf = [(xm - xm_old); y];
        
        % Mentés
        Y_out(kk) = y;
        U_out(kk) = u;
    end
    res.Y = Y_out;
    res.U = U_out;
end