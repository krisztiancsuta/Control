%% JÁRMŰDINAMIKAI LQR SZIMULÁCIÓ - MANUÁLIS INTEGRÁCIÓVAL
% Készítette: Senior MATLAB Fejlesztő
% Leírás: Laterális hiba követése LQR-rel és Feedforward kompenzációval

clear; clc; close all;

%% 1. Fizikai paraméterek
R = 20;             % Útkanyarulat sugara [m]
Caf = 222685.8 / 2; % Első tengely oldalkúszási merevsége
Car = 136242.8 / 2; % Hátsó tengely oldalkúszási merevsége
lf = 1.85;         % Tömegközéppont távolsága az első tengelytől [m]
lr = 2.43;    % Tömegközéppont távolsága a hátsó tengelytől [m]
m = 1465;           % Jármű tömege [kg]
Vx = 10;            % Haladási sebesség [m/s]
Iz = 1600;          % Tehetetlenségi nyomaték [kg*m^2]
Ts = 0.033;
%% 2. Állapottér modell (Hiba-dinamika)
% Állapotok: [e1; e1_dot; e2; e2_dot]
A = [0 1                            0                      0;...
     0 -(2*Caf+2*Car)/(m*Vx)        (2*Caf+2*Car)/m        (-2*lf*Caf+2*lr*Car)/(m*Vx);...
     0 0                            0                      1;...
     0 -(2*lf*Caf-2*lr*Car)/(Iz*Vx) (2*lf*Caf-2*lr*Car)/Iz -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*Vx)];

B1 = [0; 2*Caf/m; 0; 2*lf*Caf/Iz]; % Bemenet: kormányszög (delta)

% Zavarás mátrix (útkanyarulat miatti igényelt legyezési sebesség hatása)
B2 = [0; -(2*lf*Caf-2*lr*Car)/(m*Vx)-Vx; 0; -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*Vx)];

C = [1 0 0 0]; % Csak az e1 (laterális hiba) figyelése


ss_c = ss(A,[B1 B2],C,0)
ss_d = c2d(ss_c,Ts);

%% 3. LQR Tervezés (Kiterjesztett rendszer integrátorral)
A_ext = [ss_d.A, zeros(4, 1); -C*Ts, 1];
B_ext = [ss_d.B(:,1); 0];

Q = diag([4 0.1 1 0.1 60]); % Súlyozás: e1, de1, e2, de2, integrált hiba
r_weight = 0.1;

K = dlqr(A_ext,B_ext, Q, r_weight);
K_pd = K(1:4); % Állapot-visszacsatolási erősítés
K_i  = K(5);   % Integrátor erősítése

%% 4. Szimulációs beállítások
dt = Ts;
t = 0:dt:30;
n = length(t);
max_yaw_rate = Vx / R;

% Referencia és zavarás generálása (külső függvény hívása)
% Megjegyzés: A függvénynek léteznie kell az elérésben!
[u_vphides, u_delta_ff] = generate8likePathExt(t, max_yaw_rate, m, lf, lr, Caf, Car, Vx, K(3));


% Kezdeti feltételek [e1, e1_dot, e2, e2_dot, int_e1]
x0 = [0; 0; 0; 0; 0]; 

% Adattárolók inicializálása
states_ff = zeros(5, n);
states_no_ff = zeros(5, n);
u_ctrl_ff = zeros(1, n);
u_ctrl_no_ff = zeros(1, n);

states_ff(:, 1) = x0;
states_no_ff(:, 1) = x0;

Acl = [(ss_d.A-ss_d.B(:,1)*K(1:4)) -ss_d.B(:,1)*K(5);...
       -ss_d.C*Ts       1];

%% 5. Szimulációs hurok (Manuális Euler-integráció)
for k = 1:n-1
    % --- 1. ESET: FEEDBACK + FEEDFORWARD ---
    u_ctrl_ff(k) = -(K * states_ff(:, k))+ u_delta_ff(k);
    
    % DISZKRÉT léptetés: x[k+1] = Ad*x[k] + Bd*u[k]
    % Nem használunk "+ dx*dt"-t, mert a mátrixok már tartalmazzák az időlépést!
    states_ff(:, k+1) = Acl * states_ff(:, k) + ...
                        [ss_d.B(:,1); 0] * u_delta_ff(k) + ...
                        [ss_d.B(:,2); 0] * u_vphides(k);

    % --- 2. ESET: CSAK FEEDBACK ---
    u_ctrl_no_ff(k) = -(K * states_no_ff(:, k));
    
    states_no_ff(:, k+1) = Acl * states_no_ff(:, k) + ...
                           [ss_d.B(:,2); 0] * u_vphides(k);
end

% Utolsó vezérlő jelek számítása a grafikonhoz
u_ctrl_ff(n) = -K * states_ff(:, n) + u_delta_ff(n);
u_ctrl_no_ff(n) = -K * states_no_ff(:, n);

drawpath(t, states_ff', states_no_ff', Vx, u_vphides);
%% 6. Vizualizáció
figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], 'Color', 'w');

% Zavarás (Yaw Rate)
subplot(4,1,1);
plot(t, u_vphides, 'k', 'LineWidth', 1.5);
grid on; ylabel('\psi_{des} [rad/s]');
title('Bemeneti zavarás (Kanyarodás)');

% e1 hiba
subplot(4,1,2);
plot(t, states_no_ff(1,:), 'r--', 'LineWidth', 1.2); hold on;
plot(t, states_ff(1,:), 'b', 'LineWidth', 1.5);
grid on; ylabel('e_1 [m]');
legend('Csak Feedback', 'Feedback + Feedforward');
title('Laterális pozícióhiba');

% e2 hiba
subplot(4,1,3);
plot(t, states_no_ff(3,:), 'r--', 'LineWidth', 1.2); hold on;
plot(t, states_ff(3,:), 'b', 'LineWidth', 1.5);
grid on; ylabel('e_2 [rad]');
title('Orientációs hiba');

% Kormányszög (u)
subplot(4,1,4);
plot(t, u_ctrl_no_ff, 'r--', 'LineWidth', 1.2); hold on;
plot(t, u_ctrl_ff, 'g', 'LineWidth', 1.5);
grid on; xlabel('Idő [s]'); ylabel('\delta [rad]');
title('Beavatkozó jel (Kormányszög)');
legend('\delta_{FB}', '\delta_{FB+FF}');

sgtitle(sprintf('LQR + Feedforward összehasonlítás (Vx = %d m/s)', Vx));