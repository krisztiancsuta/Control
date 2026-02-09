%% DYNAMIC MODEL IN TERMS OF ERROR WITH RESPECT TO ROAD
% position and orientation error with respect to the road
% Phides' = Vx/R => rate of change of the desired orientation
% Vx*Vx/R = Vx*phides' => desired acceleration


% e1: the distance of the c.g. of the vehicle from the center line of
% the lane
% e1 = y - ydes
% e1' = y' + Vx(phi-phides) speed if Vx is constant (time invariant)
% else (e1' = y' + integral of Vx*e2 respect to t) time varying 
% e1'' = y'' + Vx(phi'-phides') acceleration

% e2: the orientation error of the vehicle with respect to the road.
% e2 = phi - phides


%% Jelkövető szabályzó struktúra módosítással
R = 50;  % radious of the road
Caf = 222685.8 / 2; % Cornering stiffnes az első kerékhez
Car = 136242.8 / 2; % Cornering stiffnes a hátsó kerékhez
lf = 1236e-3; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 2789e-3 - lf; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 2300; % Az autó tömege
Vx = 10; % Az autó haladási sebessége a saját koordinátarendszerében
Iz = 2873; % Az autó tehetetlenségi nyomatéka
vphides = Vx / R;
%% Az állapotvektorunk e1, e1', e2, e2' == 
% laterális pozíciócióhiba a sávközéptől,
% laterális pozíciócióhiba sebessége,
% legyezési szöghiba az úthoz képest,
% legyezési szöghiba sebessége
% A bemenet u = delta => Első kerék kormányszöge
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
% desired yaw rate as input from road radius Vx/R
B2 = [0;...
     -(2*lf*Caf-2*lr*Car)/(m*Vx)-Vx;...
     0;...
     -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*Vx)];

C = [1 0 0 0];
%% Irányíthatóság vizsgalata 
isControllable = (size(A,1) == rank(ctrb(A,B1)));

%% State feedback
% x'(t) = A*x(t) + B*u(t)
% y(t) = CT–x(t)
% u(t) = - K*x(t) + ki*z(t)
% z'(t) = r(t)- y(t)

% x'(t) = A*x(t) + B*u(t)
% z'(t) = -c*x(t) + r(t)
% u = -KT *x(t) + ki*z(t)

% Bővített állapottér
A_ext = [A, zeros(4, 1);... 
        -C, 0];
B_ext = [B1;...
         0 ];

Q = diag([4 0.1 1 0.01 60]);
r = 1;
K = lqr(A_ext,B_ext,Q,r);

%% Állapotvisszacsatolás

Acl = [(A-B1*K(1:4)) -B1*K(5);...
       -C       0];
Bcl = [B2;...
       0];

Ccl = [1 0 0 0 0;  % e1
       0 0 1 0 0]; % e2

sys_cl = ss(Acl,Bcl,Ccl,0);
%% 4. SZIMULÁCIÓ IDŐBEN VÁLTOZÓ PÁLYAPROFILLAL
% Idővektor (0-tól 30 másodpercig, hogy minden szakasz ráférjen)
t = 0:0.01:30;

max_yaw_rate = Vx / R;
% A bemeneti jel kiszámítása (Desired Yaw Rate)
% Ez hat a B2 mátrixon keresztül a rendszerre
vphides_t = generate8likePath(t,max_yaw_rate);
% Szimuláció a kanyarodási sebességgel mint bemenettel
[y_out, t_out, x_states] = lsim(sys_cl, vphides_t, t);

e1 = y_out(:,1);
e2 = y_out(:,2);

% Kormányszög kiszámítása a visszacsatolásból
K_pd = K(1:4); 
K_i  = K(5);   
% u = -K*x (mivel a referencia implicit 0 a sáv közepén)
u_total = -(K_pd * x_states(:, 1:4)' + K_i * x_states(:, 5)')';

%% 5. VIZUALIZÁCIÓ (Pályakövetésre optimalizálva)
figure('Units', 'normalized', 'Position', [0.1 0.1 0.5 0.8], 'Name', 'LQR Pályakövetés', 'Color', 'w');

% 1. részábra: Bemeneti jel (Desired Yaw Rate)
subplot(4,1,1);
plot(t, vphides_t, 'k', 'LineWidth', 1.5);
grid on; ylabel('\phi_{des} [rad/s]');
title('Bemeneti jel: Pálya görbülete (Desired Yaw Rate)');

% 2. részábra: e1 pozícióhiba
subplot(4,1,2);
plot(t, e1, 'b', 'LineWidth', 1.5);            
hold on;
yline(0, '--k', 'Cél (Sávközép)', 'LabelHorizontalAlignment', 'right');
grid on; ylabel('e_1 [m]');
title('Laterális pozícióhiba (Távolság a sávközéptől)');
legend('Aktuális hiba');

% 3. részábra: e2 orientációs hiba
subplot(4,1,3);
plot(t, e2, 'r', 'LineWidth', 1.5); % Átváltva fokba a jobb olvashatóságért
grid on; ylabel('e_2 [rad]');
title('Orientációs hiba (Szögeltérés az úthoz képest)');

% 4. részábra: Kormányszög
subplot(4,1,4);
plot(t, u_total, 'g', 'LineWidth', 1.5);
grid on; xlabel('Idő [s]'); ylabel('\delta [rad]');
title('Beavatkozó jel (Első kerék kormányszög)');
legend('\delta_{total}');