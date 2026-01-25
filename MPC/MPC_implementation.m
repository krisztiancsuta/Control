%% DYNAMIC MODEL IN TERMS OF ERROR WITH RESPECT TO ROAD
R = 1000;  % radious of the road
Caf = 80000; % Cornering stiffnes az első kerékhez
Car = 80000; % Cornering stiffnes a hátsó kerékhez
lf = 1.1; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 1.58; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 1573; % Az autó tömege
Vx = 15; % Az autó haladási sebessége a saját koordinátarendszerében
Iz = 2873; % Az autó tehetetlenségi nyomatéka
%% Contious time state space model 
A = [0 1                            0                      0;...
     0 -(2*Caf+2*Car)/(m*Vx)        (2*Caf+2*Car)/m        (-2*lf*Caf+2*lr*Car)/(m*Vx);...
     0 0                            0                      1;...
     0 -(2*lf*Caf-2*lr*Car)/(Iz*Vx) (2*lf*Caf-2*lr*Car)/Iz -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*Vx)];
% steering angle as input
B = [0;...
     2*Caf/m;...
     0;...
     2*lf*Caf/Iz];

C = [1 0 0 0];

c_ss = ss(A,B,C,0);
%% TODO
% Extend model with disturbance 

%% Discretise state space model
Ts = 0.01; % 10 ms sampling time
d_ss = c2d(c_ss,Ts);

Ap = d_ss.A;
Bp = d_ss.B;
Cp = d_ss.C;

Nc = 10;
Np = 40;
rw = 2;

%% Finding optimal solution for delta U
% Using function for calculating following matrices:
% Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e
[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e]=mpcgain(Ap,Bp,Cp,Nc,Np);
% ΔU = inv((Φ'Φ+ R))(Φ'*Rs− Φ'*F*x(ki)) 

%% Receding Horizon control
% xm(k + 1) = Ap*xm(k) + Bp*u(k)

% Simulation Setup
N_sim = 100; % Simulation steps
[m1, n1] = size(Cp);
[n, n_in] = size(Bp);

xm = [0.5; 0; 0.1; 0]; % Initial state: 0.5m lateral error, 0.1 rad yaw error
y = Cp * xm;  
u = 0;  % Initial control signal (u(k-1))

% The augmented state vector: x_e = [delta_xm; y]
% delta_xm = xm(k) - xm(k-1). Since k=0, assume delta_xm = 0
Xf = zeros(n1 + m1, 1); 

r = ones(N_sim,1); % Set point vector

% For plotting
Y = zeros(N_sim, m1);
U = zeros(N_sim, n_in);
delta_x = zeros(N_sim, n1);

R = rw*eye(Nc*n_in);
%% Receding Horizon Control Loop
% We define R outside the loop as it is constant in this case
for kk = 1:N_sim
    % 1. Calculate the optimal control sequence (Delta U)
    % deltaU sequence for the entire control horizon
    % We use the current setpoint r(kk) to scale Phi_R
    deltaU = (Phi_Phi + R) \ (Phi_R * r(kk) - Phi_F * Xf);
    
    % 2. Receding Horizon principle: Apply ONLY the first control move
    delta_u_k = deltaU(1:n_in);
    
    % 3. Update the actual control signal: u(k) = u(k-1) + delta_u(k)
    u = u + delta_u_k;
    
    % 4. Apply to the PLANT (Original State Space)
    xm_old = xm; % Store previous state for delta_xm calculation
    xm = Ap * xm + Bp * u;
    y = Cp * xm;
    
    % 5. Update the augmented state for the next iteration (Xf)
    % This is the core of the receding horizon feedback
    % Xf(k+1) = [xm(k+1)-xm(k); y(k+1)]
    Xf = [(xm - xm_old); y];
    
    % Save data for plotting
    Y(kk, :) = y;
    U(kk, :) = u;
    delta_x(kk,:) = xm - xm_old;
end

%% Plotting Results (Polished Table Layout)
k = 0:(N_sim-1);
figure('Name', 'MPC Járműdinamika Részletes Analízis'); 

% --- FELSŐ RÉSZ: Beavatkozó jel és Pozíció (Széles ábrák) ---

% 1. ábra: Beavatkozó jel (u) - A 4x2-es rács első két helyét foglalja el
subplot(4,2,[1,2])
stairs(k, U, 'k', 'LineWidth', 1.5) 
grid on;
ylabel('\delta_f [rad]')
title('Optimális kormányszög')

% 2. ábra: Kimenet (y) és Referencia - A 4x2-es rács második két helyét foglalja el
subplot(4,2,[3,4])
plot(k, Y, 'b-', 'LineWidth', 2)
hold on;
plot(k, r, 'r--', 'LineWidth', 1)
grid on;
ylabel('y [m]')
title('Oldalirányú pozícióhiba (y)')
legend('Kimenet', 'Referencia', 'Location', 'northeast')

% --- ALSÓ RÉSZ: Delta X értékek (2x2-es táblázatban) ---

% Delta y (Pozícióváltozás)
subplot(4,2,5)
plot(k, delta_x(:,1), 'g', 'LineWidth', 1.2)
grid on;
ylabel('\Delta y')
title('\Delta x_1: Pozíció vált.')

% Delta dy (Oldalirányú sebesség változása)
subplot(4,2,6)
plot(k, delta_x(:,2), 'm', 'LineWidth', 1.2)
grid on;
ylabel('\Delta dy')
title('\Delta x_2: Oldalir. seb. vált.')

% Delta psi (Irányszög hiba változása)
subplot(4,2,7)
plot(k, delta_x(:,3), 'c', 'LineWidth', 1.2)
grid on;
ylabel('\Delta \psi')
title('\Delta x_3: Irányszög vált.')

% Delta dpsi (Legyezési sebesség változása)
subplot(4,2,8)
plot(k, delta_x(:,4), 'r', 'LineWidth', 1.2)
grid on;
ylabel('\Delta d\psi')
title('\Delta x_4: Legyezési seb. vált.')

% Főcím és esztétikai javítások
sgtitle('MPC Járműirányítás - Állapotváltozók és Inkrementumok');

% Ablak méretezése, hogy minden kényelmesen elférjen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.6, 0.85]);