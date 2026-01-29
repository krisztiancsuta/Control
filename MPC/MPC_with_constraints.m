
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
Ts = 0.033; % 33 ms sampling time
d_ss = c2d(c_ss,Ts);

Ap = d_ss.A;
Bp = d_ss.B;
Cp = d_ss.C;

Nc = 40;% Control Horizon up to 1320ms 
Np = 66;% Prediction Horizon up to 2178ms 
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

xm = [0; 0; 0; 0]; % Initial state
y = Cp * xm;  
u = 0;  % Initial control signal (u(k-1))

% The augmented state vector: x_e = [delta_xm; y]
% delta_xm = xm(k) - xm(k-1). Since k=0, assume delta_xm = 0
Xf = zeros(n1 + m1, 1); 
%% Sawtooth signal as reference
period = 40;       % Hány lépés alatt érjen körbe egy teljes háromszög ciklus
amplitude = 0.5;   % A háromszög csúcsértéke (méterben)

t_vec = (0:N_sim-1)';
r = amplitude * (2/pi) * asin(sin(2*pi*t_vec/period));

%% For plotting
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

du_max = 0.05; % Max rate of steering angle change
du_min = -0.05; % Min rate of steering angle change

dUmax = ones(Nc, 1) * du_max;
dUmin = ones(Nc, 1) * (du_min);

% Kimeneti korlátok definiálása (ha az út széle pl. +/- 1 méter)
Ymax = ones(Np, 1) * 4; 
Ymin = ones(Np, 1) * -4;

I_n_in = eye(n_in);
C1 = repmat(I_n_in, Nc, 1);
C2 = kron(tril(ones(Nc)), I_n_in);

% The 3 contstraints enaquility can be described as:
M1 = [-C2;C2];
M2 = [-eye(Nc*n_in); eye(Nc*n_in)];
M3 = [-Phi;Phi];
M= [M1;M2;M3];
% M*ΔU <= gamma

DU = zeros(N_sim, n_in); % Storing dU values 
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
    DU(kk, :) = delta_u_k;
    
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

%% Plotting Results (Javított számozással)
k = 0:(N_sim-1);
figure('Name', 'MPC Járműdinamika Analízis Korlátokkal'); 

% 1. ábra: Beavatkozó jel (u) - 2 oszlop széles
subplot(6,2,[1,2])
stairs(k, U, 'k', 'LineWidth', 1.5) 
hold on;
plot(k, ones(size(k))*u_max, 'r--', 'LineWidth', 1);
plot(k, ones(size(k))*u_min, 'r--', 'LineWidth', 1);
grid on; ylabel('\delta_f [rad]'); title('Optimális kormányszög (u) és korlátok');
legend('u', 'u_{max/min}');

% 2. ábra: Kimenet (y) és Referencia (Háromszögjel)
subplot(6,2,[3,4])
plot(k, Y, 'b-', 'LineWidth', 2); hold on;
plot(k, r, 'r--', 'LineWidth', 1);
grid on; ylabel('y [m]'); title('Pályakövetés (Háromszög referencia)');
legend('Kimenet', 'Referencia');

% 3. ábra: Delta U diagram (Itt már nem lesz felülírva)
subplot(6,2,[5,6]) % 2 oszlop szélesre vettem, hogy jól látszódjon a dinamika
stairs(k, DU, 'b', 'LineWidth', 1.2); hold on;
plot(k, ones(size(k))*du_max, 'r:', 'LineWidth', 1.5);
plot(k, ones(size(k))*du_min, 'r:', 'LineWidth', 1.5);
grid on; ylabel('\Delta u'); title('\Delta u: Kormányszög sebesség és korlátai');
legend('\Delta u', 'du_{max/min}');

% 4. ábra: Zárt körű sajátértékek (Stabilitás)
subplot(6,2,[7,8])
th = 0:pi/50:2*pi; plot(cos(th), sin(th), 'k--'); hold on;
plot(real(lambda), imag(lambda), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on; axis equal; title('Pólusok az egységkörön');

% 5. ábra: Állapotváltozók növekményei (Delta X)
subplot(6,2,9); plot(k, delta_x(:,1), 'g'); grid on; title('\Delta y');
subplot(6,2,10); plot(k, delta_x(:,2), 'm'); grid on; title('\Delta dy');
subplot(6,2,11); plot(k, delta_x(:,3), 'c'); grid on; title('\Delta \psi');
subplot(6,2,12); plot(k, delta_x(:,4), 'r'); grid on; title('\Delta d\psi');

sgtitle('MPC Járműirányítás - Szabályos háromszögjel követése');