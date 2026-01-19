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
R = 1000;  % radious of the road
Caf = 80000; % Cornering stiffnes az első kerékhez
Car = 80000; % Cornering stiffnes a hátsó kerékhez
lf = 1.1; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 1.58; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 1573; % Az autó tömege
Vx = 100; % Az autó haladási sebessége a saját koordinátarendszerében
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

% A = A + 0.001;

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

ki = 1;

% Bővített állapottér
A_ext = [A, zeros(4, 1);... 
        -C, 0];
B_ext = [B1 ;...
         0  ];

Q = diag([0.001 0.0001 0.0001 0.0001 300]);
R = 200;
K = lqr(A_ext,B_ext,Q,R);

%% Állapotvisszacsatolás

Acl = [(A-B1*K(1:4)) -B1*K(5);...
       -C       0];
Bcl = [zeros(4,1);
       1];
Ccl = [1 0 0 0 0];
sys_cl = ss(Acl,Bcl,Ccl,0);
%% ------------------------------------------------
%% 4. Szimuláció (Időben változó referenciajel)
t = 0:0.01:15;
r = zeros(size(t));
r(t > 2) = 1;      % 1s után 2m kitérés
r(t > 7) = -1;     % 7s után -1m kitérés

[y, t, x_states] = lsim(sys_cl, r, t);

% Kormányszög utólagos kiszámítása: u = -K * x_ext
u_deg = rad2deg(-K * x_states');

%% 5. VIZUALIZÁCIÓ (Egy ablakban)
figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], 'Name', 'LQR Jelkövetés Analízis');

% 5.1. Pólus-zérus térkép
subplot(2,2,1);
pzmap(sys_cl);
grid on; title('Zárt kör pólusai (Stabilitás)');

% 5.2. Jelkövetés (Pozíció)
subplot(2,2,2);
plot(t, r, 'r--', 'LineWidth', 1.2); hold on;
plot(t, y(:,1), 'b', 'LineWidth', 2);
grid on; xlabel('Idő [s]'); ylabel('Pozíció [m]');
legend('Referencia (r)', 'Jármű (e1)');
title('Laterális jelkövetés');

% 5.3. Beavatkozó jel (Kormányszög)
subplot(2,2,3);
plot(t, u_deg, 'g', 'LineWidth', 1.5);
grid on; xlabel('Idő [s]'); ylabel('Kormányszög [deg]');
title('Szabályozási parancs (\delta)');

% 5.4. Orientációs hiba (e2)
subplot(2,2,4);
plot(t, rad2deg(x_states(:,3)), 'm', 'LineWidth', 1.5);
grid on; xlabel('Idő [s]'); ylabel('Szöghiba [deg]');
title('Legyezési szöghiba (e2)');

% Szöveges visszajelzés
if isempty(findall(0,'Type','Figure')) == 0
    disp('Szimuláció sikeresen lefutott.');
end

%% 1. VÉLETLEN NÉGYSZÖGJEL GENERÁLÁSA
t_end = 30; dt = 0.01;
t = 0:dt:t_end;
r = zeros(size(t));

curr_t = 0;
while curr_t < t_end
    duration = 1.5 + rand() * 3;     % Véletlen időtartam: 1.5 - 4.5 mp
    amplitude = (rand() - 0.5) * 4;  % Véletlen amplitúdó: -2 és +2 méter között
    
    idx = (t >= curr_t) & (t < curr_t + duration);
    r(idx) = amplitude;
    curr_t = curr_t + duration;
end

%% 2. SZABÁLYOZÓK HANGOLÁSA ÉS SZIMULÁCIÓ
% Itt két különböző "karaktert" hasonlítunk össze
% Q = [e1, e1_dot, e2, e2_dot, x_int]
configs = {
    'Komfort (Lágy)', [5, 1, 5, 1, 100], 2;
    'Sport (Agresszív)', [50, 2, 50, 2, 5000], 0.05
};

figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], 'Name', 'Random Square Wave Tracking');

for i = 1:size(configs, 1)
    Q = diag(configs{i, 2});
    R = configs{i, 3};
    K = lqr(A_ext, B_ext, Q, R);
    
    % Zárt kör
    Acl = [(A - B1*K(1:4)), -B1*K(5); -C, 0];
    Bcl = [zeros(4,1); 1];
    sys_cl = ss(Acl, Bcl, [1 0 0 0 0], 0);
    
    % Szimuláció
    [y, t_sim, x_states] = lsim(sys_cl, r, t);
    u_deg = rad2deg(-K * x_states');
    
    % --- Megjelenítés ---
    subplot(2,1,1);
    plot(t_sim, y, 'LineWidth', 2); hold on;
    
    subplot(2,1,2);
    plot(t_sim, u_deg, 'LineWidth', 1.5); hold on;
end

% Grafikon szépítése
subplot(2,1,1);
plot(t, r, 'k--', 'LineWidth', 1.2);
grid on; ylabel('Pozíció (e1) [m]');
legend(configs{1,1}, configs{2,1}, 'Referencia');
title('Jelkövetés véletlen amplitúdójú és kitöltésű négyszögjellel');

subplot(2,1,2);
grid on; ylabel('Kormányszög [deg]');
xlabel('Idő [s]');
title('Beavatkozó jel (Szabályozási parancs)');
legend(configs{1,1}, configs{2,1});