air_density = 1.225;
Cd = 0.3;
Am = 2.2;
v0 = 25;
b_aero = air_density*Cd*Am*v0;
m  = 1000;
b = 50;

A = -b/m;
B = 1/m;
C = 1; %

c_ss = ss(A,B,C,0);
%% Discrete model
Ts = 33e-3; % 33 ms sampling time
d_ss = c2d(c_ss,Ts);

Ap = d_ss.A;
Bp = d_ss.B; % B matrix for Fv angle
Cp = d_ss.C;

Nc = 20;
Np = 60; 
rw = 0.00001;

[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e,F,Phi]=mpcgain(Ap,Bp,Cp,Nc,Np);
%% Receding Horizon control
% xm(k + 1) = Ap*xm(k) + Bp*u(k)
t_vec = 0:Ts:35.7; 

N_sim = length(t_vec);
[m1, n1] = size(Cp);
[n, n_in] = size(B_e);

xm = 0; % Initial state
y = Cp * xm;  
u = 0;  % Initial control signal (u(k-1))

% The augmented state vector: x_e = [delta_xm; y]
% delta_xm = xm(k) - xm(k-1). Since k=0, assume delta_xm = 0
Xf = zeros(n1 + m1, 1);
%% Reference Vx signal
r = zeros(size(t_vec));
u_ff = zeros(size(t_vec));
u_ff(t_vec >= 10) = 1/2 * Cd * air_density * Am *25^2;
u_ff(t_vec >= 20) = 1/2 * Cd * air_density * Am *5^2;
r(t_vec >= 10) = 25; 
r(t_vec >= 20) = 5; 

%% For plotting
Y = zeros(N_sim, m1);
U = zeros(N_sim, n_in);
DU = zeros(N_sim, n_in); % Storing dU values 
delta_x = zeros(N_sim, n1);

R = rw*eye(Nc*n_in)
%% Calculate K vector for comparison with LQR
K_full = (Phi_Phi + R) \ Phi_F;
Kmpc = K_full(1:n_in, :)
Ky = Kmpc(:, end-m1+1:end)
Acl = A_e - B_e * Kmpc;
lambda=eig(Acl);
%% Constraints
u_max = 10000; 
u_min = -10000; 

Umax = ones(Nc, 1) * u_max;
Umin = ones(Nc, 1) * u_min;

du_max = 2000;
du_min = -2000; 

dUmax = ones(Nc, 1) * du_max;
dUmin = ones(Nc, 1) * du_min;

I_n_in = eye(n_in);
C1 = repmat(I_n_in, Nc, 1);
C2 = kron(tril(ones(Nc)), I_n_in);

% The 3 contstraints enaquility can be described as:
M1 = [-C2;C2];
M2 = [-eye(Nc*n_in); eye(Nc*n_in)];
M= [M1;M2];
% M*ΔU <= gamma

DU = zeros(N_sim, n_in); % Storing dU values 

H = (Phi_Phi + R);
%% Receding Horizon Control Loop
% We define R outside the loop as it is constant in this case
for kk = 1:N_sim
 
    N1 = [-Umin + C1*u;... 
           Umax - C1*u];
    N2 = [-dUmin;dUmax];
    gamma = [N1;N2];

    % 1. Calculate the optimal control sequence (Delta U)
    % deltaU sequence for the entire control horizon
    % We use the current setpoint r(kk) to scale Phi_R
    deltaU_unconstrained = (Phi_Phi + R) \ (Phi_R * r(kk) - Phi_F * Xf);
    
    f = -(Phi_R * r(kk) - Phi_F * Xf);

    deltaU = QPhild(H, f, M, gamma);
   

    % 2. Receding Horizon principle: Apply ONLY the first control move
    delta_u_k = deltaU(1:n_in);
    DU(kk, :) = delta_u_k;
    
    % 3. Update the actual control signal: u(k) = u(k-1) + delta_u(k)
    u = u + delta_u_k + u_ff(kk);
    
    u = max(min(u, u_max), u_min);
  

    % 4. Apply to the PLANT (Original State Space)
    xm_old = xm; % Store previous state for delta_xm calculation
    xm = Ap * xm + Bp * u ;
    y = Cp * xm;
    
    % 5. Update the augmented state for the next iteration (Xf)
    % This is the core of the receding horizon feedback
    % Xf(k+1) = [xm(k+1)-xm(k); y(k+1)]
    Xf = [(xm - xm_old); y];
    
    % Save data for plotting
    Y(kk, :) = y;
    U(kk) = u;
    delta_x(kk,:) = xm - xm_old;
end
%% Ábrázolás
figure('Color', 'w', 'Name', 'MPC Járműirányítás');

% 1. Alábra: Sebesség
subplot(2,1,1);
stairs(t_vec, Y, 'LineWidth', 2, 'DisplayName', 'MPC Sebesség (v_x)');
hold on;
stairs(t_vec, r, '--r', 'LineWidth', 1.5, 'DisplayName', 'Referencia');
grid on;
ylabel('Sebesség [m/s]');
title('MPC alapú sebességkövetés');
legend('Location', 'southeast');

% 2. Alábra: Pedálpozíció
% Mappelés a konzulens szerint: F / 10000
pedal = U / 10000;

subplot(2,1,2);
stairs(t_vec, pedal, 'Color', [0.466 0.674 0.188], 'LineWidth', 2);
grid on;
ylabel('Pedálpozíció (F/10000)');
xlabel('Idő [s]');
title('Beavatkozó jel (Gáz/Fék)');
ylim([-1.1 1.1]);