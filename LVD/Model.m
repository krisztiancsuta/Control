% The two major elements of the longitudinal vehicle model are the
% vehicle dynamics and the powertrain dynamics.
% 1.The vehicle dynamics are influenced by longitudinal tire forces, aerodynamic
% drag forces, rolling resistance forces and gravitational forces
% 2.The longitudinal powertrain system of the vehicle
% consists of the internal combustion engine, the torque converter, the
% transmission and the wheels
%% Parameters
% Fxf is the longitudinal tire force at the front tires
% Fxr, is the longitudinal tire force at the rear tires
% Faero,, is the equivalent longitudinal aerodynamic drag force
% Rxf is the force due to rolling resistance at the front tires
% Rxr is the force due to rolling resistance at the rear tires
% m is the mass of the vehicle
% g is the acceleration due to gravity
% Phi is the angle of inclination of the road on which the vehicle is travelingm = 1000;

% state vector is deltavx = vx-vx0
% ax = -b/m * Vx + 1/m * Fv
% FFWD tag -> Faero at Vx0 speed

% u = the force that pulls the vehicle forward

% we want to find ax
air_density = 1.225;
Cd = 0.3;
Am = 2.2;
v0 = 25;
b_aero = air_density*Cd*Am*v0;
m  = 1000;
b = 50;

A = -b/m;
B = 1/m;
C = 1; % Output is vx
D = 0;

Aaug = [A 0;
        -C 0];
Baug = [B;...
        0];

Q = diag([100, 50]);
R = 1/1e10; % 1/1e5

K = lqr(Aaug,Baug,Q,R);

%% Állapotvisszacsatolás

Acl = [(A-B*K(1)) -B*K(2);...
       -C       0];
Bcl = [0;...
       1];

Ccl = [1 0];

sys_cl = ss(Acl,Bcl,Ccl,0);
%% 4. SZIMULÁCIÓ IDŐBEN VÁLTOZÓ PÁLYAPROFILLAL
t = 0:0.01:30;
v_ref_val = 10;
v_ref = zeros(size(t));          
u_ff = zeros(size(t));

v_ref(t >= 1) = v_ref_val;       
% FFWD számítása: a célsebességhez tartozó légellenállás
u_ff(t >= 1) = -1/2 * Cd * air_density * Am * v_ref_val^2;

% Referencia mátrix: 1. oszlop az integrátornak (v_ref), 2. oszlop a gyorsulásnak (u_ff)
r = [v_ref']; 

[y_out, t_out, x_states] = lsim(sys_cl, r, t);
vx = y_out(:,1);
%% 5. ÁBRÁZOLÁS - JAVÍTOTT PEDÁL SZÁMÍTÁSSAL
figure('Color', 'w', 'Name', 'LQR + FFWD Sebességkövetés');

% 1. Alábra: Sebesség
subplot(2,1,1);
plot(t, vx, 'LineWidth', 2, 'DisplayName', 'Mért sebesség (v_x)');
hold on;
plot(t, v_ref, '--r', 'LineWidth', 1.5, 'DisplayName', 'Referencia (r)');
grid on; ylabel('Sebesség [m/s]'); title('Jármű sebességkövetése');
legend('Location', 'southeast');

% 2. Alábra: Pedálpozíció (Feedback + FFWD összege!)
% u_feedback = -K1 * delta_v - K2 * integral_hiba
F_feedback = -K(1)*x_states(:,1) - K(2)*x_states(:,2);
F_total = F_feedback + u_ff'; % Itt adjuk hozzá a FFWD ágat!

pedal = F_total / 10000;

subplot(2,1,2);
plot(t, pedal, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);
grid on; ylabel('Pedálpozíció'); xlabel('Idő [s]');
title('Teljes beavatkozó jel (LQR feedback + FFWD)');
ylim([-1.1 1.1]); % Realitás tartomány