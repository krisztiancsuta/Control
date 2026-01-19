%% DYNAMIC MODEL IN TERMS OF ERROR WITH RESPECT TO ROAD
R = 1000;  % radious of the road
Caf = 80000; % Cornering stiffnes az első kerékhez
Car = 80000; % Cornering stiffnes a hátsó kerékhez
lf = 1.1; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 1.58; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 1573; % Az autó tömege
Vx = 30; % Az autó haladási sebessége a saját koordinátarendszerében
Iz = 2873; % Az autó tehetetlenségi nyomatéka

% Az állapotvektorunk e1, e1', e2, e2' == 
% laterális pozíciócióhiba a sávközéptől,
% laterális pozíciócióhiba sebessége,
% legyezési szöghiba az úthoz képest,
% legyezési szöghiba sebessége
% phides = desired yaw rate determined from road radius R Vx/R

A = [0 1                            0                      0;...
     0 -(2*Caf+2*Car)/(m*Vx)        (2*Caf+2*Car)/m        (-2*lf*Caf+2*lr*Car)/(m*Vx);...
     0 0                            0                      1;...
     0 -(2*lf*Caf-2*lr*Car)/(Iz*Vx) (2*lf*Caf-2*lr*Car)/Iz -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*Vx)];

% steering angle as input
B1 = [0;...
     2*Caf/m;...
     0;...
     2*lf*Caf/Iz];
% desired yaw rate as disturbance from road radius Vx/R
B2 = [0;...
     -(2*lf*Caf-2*lr*Car)/(m*Vx)-Vx;...
     0;...
     -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*Vx)];

C = [1 0 0 0];
% Irányíthatóság vizsgalata 
isControllable = (size(A,1) == rank(ctrb(A,B1)));

Q = diag([10 0.1 0.1 10]);
r = 1;
K = lqr(A,B1,Q,r);
%% STEADY STATE ERROR FROM DYNAMIC EQUATIONS
% Szeretnénk hogy a steady state error nullához tartson, ha van a kanyar
% miatti külső zavaró hatásunk, ezt egy előre csatolt (feedforward) taggal
% kompenzáljuk
% delta = -K * x(t) + delta_ff => feedforward miatt a kompenzálás a
% kormányszögben
% Ekkor a closed loop ilyen alakban szerepel
% x'(t) = (A-B1*K)*x(t) + B1*delta_ff + B2*vphides

% 1. Teljes tengelytáv (Wheelbase)
L = lf + lr;

% 2. Tengelyterhelések számítása (statikus terheléseloszlás)
% m_f: első tengelyre eső tömeg, m_r: hátsó tengelyre eső tömeg
m_f = m * (lr / L); 
m_r = m * (lf / L);

% 3. Alulkormányzottsági gradiens (Understeer Gradient - Kv)
% Ez határozza meg a jármű kanyarodási viselkedését (Rajamani Eq. 3.12/Table 3-1)
Kv = (m_f / (2 * Caf)) - (m_r / (2 * Car));

% 4. Oldalgyorsulás (Lateral Acceleration - ay)
% Mivel vphides = Vx / R, ezért R = Vx / vphides
% ay = Vx^2 / R = Vx * (Vx/R) = Vx * vphides
% A szimulációban a kanyar sugara R = 1000m, így:
ay = Vx^2 / R;

% e2 steady state hibaja
e2_ss = -lr/R + (lf*m*Vx*Vx)/((2*Car*L)*R);

% 5. Feedforward kormányszög (delta_ff)
% Ackermann szög (L/R) + dinamikus kompenzáció (Kv * ay)
delta_ff_val = (L / R) + Kv * ay + K(3)*e2_ss;

% Kiírás ellenőrzésképpen (fokban)
fprintf('Szükséges Feedforward szög: %.4f fok\n', rad2deg(delta_ff_val));

%% Állapotvisszacsatolás
% delta = -K * x(t)
% Az openloopra zavaras nelkul szamolom a K erositeseket
% A zavarasnak csak a zart korben lesz szerepe
% Állapotvisszacsatolás

Acl = A-B1*K;...
Bcl = [B1,B2]; % A zavarás mint bemenet jelenik meg a rendszerben
Ccl = [1 0 0 0];
sys_cl = ss(Acl,Bcl,Ccl,0);

%% 4. JAVÍTOTT SZIMULÁCIÓ
t = 0:0.01:15;
vphides_val = Vx / R;

u_vphides = zeros(size(t));
u_delta_ff = zeros(size(t));
u_vphides(t > 2) = vphides_val;
u_delta_ff(t > 2) = delta_ff_val;

% Szimuláció
u_sim = [u_delta_ff', u_vphides'];
[~, ~, x_states] = lsim(sys_cl, u_sim, t);

e1 = x_states(:,1);
e2 = x_states(:,3);
u_feedback = -K * x_states';
u_total = u_feedback' + u_delta_ff';

%% 5. VIZUALIZÁCIÓ (Radián alapú feliratokkal)
figure('Units', 'normalized', 'Position', [0.2 0.1 0.6 0.8], 'Name', 'LQR Zavarválasz (Radián)', 'Color', 'w');

% 1. részábra: Desired Yaw Rate
subplot(4,1,1);
plot(t, u_vphides, 'k', 'LineWidth', 1.5);
grid on; ylabel('[rad/s]');
title('Bemeneti zavarás: \psi_{des}^\prime (Yaw rate)');

% 2. részábra: e1 pozícióhiba
subplot(4,1,2);
plot(t, e1, 'b', 'LineWidth', 1.5);
hold on; yline(0, '--k', 'Alpha', 0.5);
yline(e1(end), '--r', ['Maradó hiba: ', num2str(e1(end), '%.4f'), ' m']);
grid on; ylabel('e_1 [m]');
title('Laterális pozícióhiba');

% 3. részábra: e2 orientációs hiba (Radiánban)
subplot(4,1,3);
plot(t, e2, 'r', 'LineWidth', 1.5);
hold on;
yline(e2(end), '--k', ['Maradó hiba: ', num2str(e2(end), '%.4f'), ' rad']);
grid on; ylabel('e_2 [rad]');
title('Orientációs hiba');

% 4. részábra: Kormányszög (Radiánban)
subplot(4,1,4);
plot(t, u_total, 'g', 'LineWidth', 1.5);
grid on; xlabel('Idő [s]'); ylabel('\delta [rad]');
title('Szabályozó beavatkozása (Kormányszög)');
legend('Visszacsatolt vezérlőjel');

sgtitle(['LQR Zavar elnyomás (Integrátor nélkül) - Minden mértékegység radián: R = ', num2str(R), 'm']);