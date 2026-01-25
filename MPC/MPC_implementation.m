%% State space model 
Ap = 0.8;
Bp = 0.1;
Cp = 1;
Nc = 4;
Np = 10;
rw = 0.1;

%% Finding optimal solution for delta U
% Using function for calculating following matrices:
% Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e
[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e]=mpcgain(Ap,Bp,Cp,Nc,Np);

% ΔU = inv((Φ'Φ+ R))(Φ'*Rs− Φ'*F*x(ki)) 
[~ , n_in] = size(Bp);
R = rw*eye(Nc*n_in);

%% Receding Horizon control
% xm(k + 1) = Ap*xm(k) + Bp*u(k)

% Simulation Setup
N_sim = 50; % Simulation steps
[m1, n1] = size(Cp);
[~, n_in] = size(Bp);

xm = 0; % Initial plant state (original system)
y = 0;  % Initial output
u = 0;  % Initial control signal (u(k-1))

% The augmented state vector: x_e = [delta_xm; y]
% delta_xm = xm(k) - xm(k-1). Since k=0, assume delta_xm = 0
Xf = zeros(n1 + m1, 1); 

% For plotting
Y_history = zeros(N_sim, m1);
U_history = zeros(N_sim, n_in);
r = ones(N_sim,1); % Set point vector

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
    Y_history(kk, :) = y;
    U_history(kk, :) = u;
end

%% Plotting Results
figure;
subplot(2,1,1);
plot(1:N_sim, Y_history, 'b', 'LineWidth', 2);
hold on;
plot(1:N_sim, r, 'r--', 'LineWidth', 1.5);
title('Output (y) - Setpoint Tracking');
legend('Actual Output', 'Setpoint');
grid on;

subplot(2,1,2);
stairs(1:N_sim, U_history, 'k', 'LineWidth', 2);
title('Control Signal (u)');
grid on;