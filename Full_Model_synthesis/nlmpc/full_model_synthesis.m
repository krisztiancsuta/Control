%% Lane Following NLMPC Szimuláció - Változó útvonal (Nyolcas pálya)

% 1. Alapbeállítások (maradnak a korábbiak)
nx = 5; ny = 2; nu = 2;
Ts = 0.033;
nlobj = nlmpc(nx, ny, nu);
nlobj.Ts = Ts;
nlobj.PredictionHorizon = 10;
nlobj.ControlHorizon = 2;

nlobj.Model.StateFcn = "nonlinear_model";
nlobj.Jacobian.StateFcn = "myStateJacobian";
nlobj.Model.OutputFcn = "myOutputFunction";
nlobj.Jacobian.OutputFcn = "myOutputJacobian";
nlobj.Model.NumberOfParameters = 1;

% 2. Korlátok és skálázás (a dokumentáció alapján beállítva)
nlobj.MV(1).Min = -0.8; nlobj.MV(1).Max = 0.8; 
nlobj.MV(2).Min = -10000; nlobj.MV(2).Max = 10000;
nlobj.OV(1).ScaleFactor = 0.5; nlobj.OV(2).ScaleFactor = 15;
nlobj.MV(1).ScaleFactor = 2.26; nlobj.MV(2).ScaleFactor = 10000;
nlobj.Weights.OutputVariables = [1, 0.2];
nlobj.Weights.ManipulatedVariablesRate = [0.1, 0.3];

% 3. ÚTVONAL GENERÁLÁSA
T_sim = 30; % Hosszabb szimuláció, hogy látszódjon a kanyarodás
N_steps = round(T_sim / Ts);
tHistory = (0:N_steps) * Ts;

% Desired yaw rate generálása a nyolcas pálya függvényeddel
% max_yaw_rate = 0.1 rad/s (kb. 5.7 fok/s)
max_yaw_rate = 0.3;
psi_dot_ref_series = generate8likePath(tHistory, max_yaw_rate);

% 4. SZIMULÁCIÓ INDÍTÁSA
x0 = [0; 0; 0; 0; 1.1]; % Induljunk 0 hibával és 10 m/s sebességgel
u0 = [0; 0];
yref = [0, 25];

xHistory = zeros(nx, N_steps + 1);
xHistory(:,1) = x0;
uHistory = zeros(nu, N_steps);

options = nlmpcmoveopt;
xk = x0;
uk = u0;

disp('Szimuláció futtatása változó yaw rate mellett...');
for k = 1:N_steps
    % Aktuális referencia yaw rate lekérése az időpillanathoz
    current_psi_dot_des = psi_dot_ref_series(k);
    
    % Átadjuk az MPC-nek az aktuális paramétert
    options.Parameters = {current_psi_dot_des};
    
    % MPC számítás
    [uk, options, info] = nlmpcmove(nlobj, xk, uk, yref, [], options);    
    uHistory(:,k) = uk;
    
    % Rendszer léptetése ode45-tel a változó paraméterrel
    [~, x_temp] = ode45(@(t,x) nonlinear_model(x, uk, current_psi_dot_des), [0 Ts], xk);
    xk = x_temp(end, :)'; 
    
    % Fizikai korlát
    xk(5) = max(xk(5), 0.1); 
    xHistory(:,k+1) = xk;
end
disp('Szimuláció kész!');

%% --- EREDMÉNYEK KIRAJZOLÁSA ---

figure('Name', 'NLMPC Nyolcas pálya követési eredmények', 'NumberTitle', 'off');

% 1. Keresztirányú pozícióhiba (e1)
subplot(3,2,1);
plot(tHistory, xHistory(1,:), 'b', 'LineWidth', 1.5);
yline(0, 'r--', 'Referencia');
title('Keresztirányú pozícióhiba (e1)');
xlabel('Idő [s]'); ylabel('e1 [m]');
grid on;

% 2. Hosszirányú sebesség (vx)
subplot(3,2,2);
plot(tHistory, xHistory(5,:), 'b', 'LineWidth', 1.5);
yline(yref(2), 'r--', 'Referencia');
title('Hosszirányú sebesség (vx)');
xlabel('Idő [s]'); ylabel('vx [m/s]');
grid on;

% 3. Kormányzási szög beavatkozás (u1 - delta)
% Itt látható, hogyan kormányoz az MPC a kanyarokban
subplot(3,2,3);
stairs(tHistory(1:end-1), uHistory(1,:), 'k', 'LineWidth', 1.5);
hold on;
% Összehasonlításképp kitehetjük a kért yaw rate-et kicsinyítve, 
% hogy lássuk az összefüggést (opcionális)
title('Kormányzási szög bemenet (\delta)');
xlabel('Idő [s]'); ylabel('\delta [rad]');
grid on;

% 4. Gyorsító/Lassító erő beavatkozás (u2 - F)
subplot(3,2,4);
stairs(tHistory(1:end-1), uHistory(2,:), 'k', 'LineWidth', 1.5);
title('Hajtó/Fékező erő (F)');
xlabel('Idő [s]'); ylabel('F [N]');
grid on;

% 5. Yaw Rate követés (Referencia vs. Aktuális)
subplot(3,2,5:6); % Alsó sor teljes szélességben
actual_yaw_rate = xHistory(4,1:end-1) + psi_dot_ref_series(1:end-1);
plot(tHistory, psi_dot_ref_series, 'r--', 'LineWidth', 1.5); hold on;
plot(tHistory(1:end-1), actual_yaw_rate, 'b', 'LineWidth', 1);
title('Yaw Rate: Előírt (Pálya) vs. Megvalósult');
xlabel('Idő [s]'); ylabel('\psi_{dot} [rad/s]');
legend('Kívánt (Pálya)', 'Autó válasza');
grid on;