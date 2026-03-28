%% MATLAB Szimuláció: LQR+FFWD vs. Csak LQR (For-loop)
clear; clc;

%% 1. Rendszer paraméterek
air_density = 1.225;
Cd = 0.3;
Am = 2.2;
m  = 1000;
b = 50;

% Alaprendszer
A = -b/m;
B = 1/m;
C = 1;

% Kiterjesztett rendszer (Integrátorral)
Aaug = [A 0; -C 0];
Baug = [B; 0];

% LQR Tervezés
Q = diag([100, 50]);
R = 0.001;
K = lqr(Aaug, Baug, Q, R);

%% 2. Szimulációs beállítások
dt = 0.033;
t = 0:dt:50;
n = length(t);

v_ref_val = 30;
v_ref = (t >= 1) * v_ref_val; 
% Előrecsatolt jel (FFWD) - csak az első esethez

% Tárolók az 1. esethez (LQR + FFWD)
x_ff = zeros(2, n);
y_ff = zeros(1, n);
u_total_ff = zeros(1, n);

% Tárolók a 2. esethez (Csak LQR)
x_no = zeros(2, n);
y_no = zeros(1, n);
u_total_no = zeros(1, n);

%% 3. Szimulációs hurok
for k = 1:n-1
    % --- 1. ESET: LQR + FEEDFORWARD ---
    v_act = x_ff(1,k);
    u_ff_dynamic = 0.5 * Cd * air_density * Am * v_act^2;

    u_total_ff(k) = -K(1)*x_ff(1,k) - K(2)*x_ff(2,k) + u_ff_dynamic;
    u_total_ff(k) = max(0, min(10000, u_total_ff(k))); % Hard limit példa

    dx_ff = Aaug * x_ff(:,k) + Baug * u_total_ff(k) + [0; 1] * v_ref(k);
    x_ff(:, k+1) = x_ff(:, k) + dx_ff * dt;
    y_ff(k) = x_ff(1,k);

    % --- 2. ESET: CSAK LQR ---
    u_total_no(k) = -K(1)*x_no(1,k) - K(2)*x_no(2,k); % Nincs + u_ff!
    dx_no = Aaug * x_no(:,k) + Baug * u_total_no(k) + [0; 1] * v_ref(k);
    x_no(:, k+1) = x_no(:, k) + dx_no * dt;
    y_no(k) = x_no(1,k);
    u_total_no(k) = max(0, min(10000, u_total_no(k))); % Hard limit példa
end

% Utolsó értékek mentése
y_ff(n) = x_ff(1,n);
y_no(n) = x_no(1,n);

%% 4. Összehasonlító ábrázolás
figure('Color', 'w', 'Position', [100 100 800 600]);

% Sebességkövetés
subplot(2,1,1);
plot(t, v_ref, '--r', 'LineWidth', 1.2, 'DisplayName', 'Referencia');
hold on;
plot(t, y_ff, 'Color', [0 0.447 0.741], 'LineWidth', 2, 'DisplayName', 'LQR + Feedforward');
plot(t, y_no, 'Color', [0.85 0.325 0.098], 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', 'Csak LQR');
grid on;
ylabel('Sebesség [m/s]');
title('Jármű sebességkövetése: Előrecsatolás hatása');
legend('Location', 'southeast');

% Beavatkozó jel (Pedálpozíció)
subplot(2,1,2);
plot(t, u_total_ff / 10000, 'Color', [0.466 0.674 0.188], 'LineWidth', 2, 'DisplayName', 'LQR + Feedforward');
hold on;
plot(t, u_total_no / 10000, 'Color', [0.635 0.078 0.184], 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Csak LQR');
grid on;
ylabel('Pedálpozíció [-]');
xlabel('Idő [s]');
title('Beavatkozó jelek összehasonlítása');
legend;