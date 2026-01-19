%% DYNAMIC MODEL IN TERMS OF ERROR WITH RESPECT TO ROAD
R = 1000;  % radious of the road
Caf = 80000; % Cornering stiffnes az első kerékhez
Car = 80000; % Cornering stiffnes a hátsó kerékhez
lf = 1.1; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 1.58; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 1573; % Az autó tömege
Vx = 30; % Az autó haladási sebessége a saját koordinátarendszerében
Iz = 2873; % Az autó tehetetlenségi nyomatéka
vphides = Vx / R;

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

Q = diag([30 0 0 5]);
r = 1;
K = lqr(A,B1,Q,r);
%% Állapotvisszacsatolás
% delta = -K * x(t)
% Az openloopra zavaras nelkul szamolom a K erositeseket
% A zavarasnak csak a zart korben lesz szerepe

Acl = A-B1*K;...
Bcl = B2; % A zavarás mint bemenet jelenik meg a rendszerben, ez marado hiat fog okozni
Ccl = [1 0 0 0];
sys_cl = ss(Acl,Bcl,Ccl,0);
%% 4. Szimuláció (Minden radiánban)
t = 0:0.01:15;
vphides_val = Vx / R; 

% Bemeneti zavarás (kanyar t=2s-nél)
u_vphides = zeros(size(t));
u_vphides(t > 2) = vphides_val; 

sys_sim = ss(Acl, Bcl, eye(4), 0);
[~, ~, x_states] = lsim(sys_sim, u_vphides, t);

% Adatok kinyerése
e1 = x_states(:,1);
e2 = x_states(:,3); % Radiánban marad
u_rad = -K * x_states'; 

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
plot(t, u_rad, 'g', 'LineWidth', 1.5);
grid on; xlabel('Idő [s]'); ylabel('\delta [rad]');
title('Szabályozó beavatkozása (Kormányszög)');
legend('Visszacsatolt vezérlőjel');

sgtitle(['LQR Zavar elnyomás (Integrátor nélkül) - Minden mértékegység radián: R = ', num2str(R), 'm']);