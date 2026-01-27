
%% Constraints types
% rate of change of the control variables (Δu(k))
%   Δu(k)min <= Δu(k)<= Δu(k)max
%   we use it where the rate of change is limited

% Constraints on the Amplitude of the Control Variable
% u(k)min <= u(k)<= u(k)max
% These are the physical hard constraints on the system.
% we need to pay particular attention to the fact that u(k) is an incre-
% mental variable, not the actual physical variable.
% The actual physical control
% variable equals the incremental variable u plus its steady-state value uss.
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

Nc = 4;
Np = 20;
rw = 0.1;

%% Finding optimal solution for delta U
% Using function for calculating following matrices:
% Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e
[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e,F,Phi]=mpcgain(Ap,Bp,Cp,Nc,Np);
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
%% Calculate K vector for comparison with LQR
K_full = (Phi_Phi + R) \ Phi_F;
Kmpc = K_full(1:n_in, :);
Ky = Kmpc(:, end-m1+1:end);
Acl = A_e - B_e * Kmpc;
lambda=eig(Acl);



%% Constraints
% deltau vector called as decision variable
% Δumin <= Δu(ki)<= Δumax
% u(ki) = u(ki-1) + Δu(ki) = u(ki-1) + [1 0 ... 0]ΔU
% u(ki + 1) = u(ki) + Δu(ki + 1) = u(ki) + [0 1... 0]ΔU
% u(ki + 1) = u(ki-1) + [1 1 ... 0]ΔU
% U = C1*u(ki-1) + C2*ΔU
% These matrices can handle MIMO systems
% Steady state steering angle is 0;
u_max = pi/4; % Max steering angle
u_min = -pi/4; % Min steering angle

Umax = ones(Nc, 1) * u_max;
Umin = ones(Nc, 1) * (u_min);

du_max = 0.2; % Max rate of steering angle change
du_min = -0.2; % Min rate of steering angle change

dUmax = ones(Nc, 1) * du_max;
dUmin = ones(Nc, 1) * (du_min);

% Kimeneti korlátok definiálása (ha az út széle pl. +/- 1 méter)
Ymax = ones(Np, 1) * 1; 
Ymin = ones(Np, 1) * -1;

I_n_in = eye(n_in);
C1 = repmat(I_n_in, Nc, 1);
C2 = kron(tril(ones(Nc)), I_n_in);

% The 3 contstraints enaquility can be described as:
M1 = [-C2;C2];
M2 = [-eye(Nc*n_in); eye(Nc*n_in)];
M3 = [-Phi;Phi];
M= [M1;M2;M3];
% M*ΔU <= gamma
%% Receding Horizon Control Loop
% We define R outside the loop as it is constant in this case
for kk = 1:N_sim
 
    N1 = [-Umin + C1*u;... 
           Umax - C1*u];
    N2 = [-dUmin;dUmax];
    N3 = [-Ymin + F*Xf;... 
           Ymax - F*Xf];
    
    gamma = [N1;N2;N3];

    % 1. Calculate the optimal control sequence (Delta U)
    % deltaU sequence for the entire control horizon
    % We use the current setpoint r(kk) to scale Phi_R
    deltaU_unconstrained = (Phi_Phi + R) \ (Phi_R * r(kk) - Phi_F * Xf);
    

    H = (Phi_Phi + R);
    f = -(Phi_R * r(kk) - Phi_F * Xf);

    deltaU = QPhild(H, f, M, gamma);
   

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

%% Plotting Results (Kibővített Analízis Pólusokkal)
k = 0:(N_sim-1);
figure('Name', 'MPC Járműdinamika Részletes Analízis'); 

% 1. ábra: Beavatkozó jel (u)
subplot(5,2,[1,2])
stairs(k, U, 'k', 'LineWidth', 1.5) 
grid on;
ylabel('\delta_f [rad]')
title('Optimális kormányszög')

% 2. ábra: Kimenet (y) és Referencia
subplot(5,2,[3,4])
plot(k, Y, 'b-', 'LineWidth', 2)
hold on;
plot(k, r, 'r--', 'LineWidth', 1)
grid on;
ylabel('y [m]')
title('Oldalirányú pozícióhiba (y)')
legend('Kimenet', 'Referencia', 'Location', 'northeast')

% --- ÚJ RÉSZ: Zárt körű sajátértékek (Pólusok) megjelenítése ---
subplot(5,2,[5,6])
% Egységkör megrajzolása
th = 0:pi/50:2*pi;
plot(cos(th), sin(th), 'k--', 'LineWidth', 0.5); hold on;
% Pólusok ábrázolása
plot(real(lambda), imag(lambda), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
axis equal;
xlabel('Valós rész'); ylabel('Képzetes rész');
title('Zárt körű sajátértékek (Pólusok) az egységkörön');
legend('Egységkör', 'Pólusok', 'Location', 'northeast');

% --- ALSÓ RÉSZ: Delta X értékek ---
% Delta y
subplot(5,2,7)
plot(k, delta_x(:,1), 'g', 'LineWidth', 1.2)
grid on; ylabel('\Delta y'); title('\Delta x_1: Pozíció vált.')

% Delta dy
subplot(5,2,8)
plot(k, delta_x(:,2), 'm', 'LineWidth', 1.2)
grid on; ylabel('\Delta dy'); title('\Delta x_2: Oldalir. seb. vált.')

% Delta psi
subplot(5,2,9)
plot(k, delta_x(:,3), 'c', 'LineWidth', 1.2)
grid on; ylabel('\Delta \psi'); title('\Delta x_3: Irányszög vált.')

% Delta dpsi
subplot(5,2,10)
plot(k, delta_x(:,4), 'r', 'LineWidth', 1.2)
grid on; ylabel('\Delta d\psi'); title('\Delta x_4: Legyezési seb. vált.')

sgtitle('MPC Járműirányítás - Állapotok, Növekmények és Stabilitás');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.05, 0.6, 0.95]);