
%% DYNAMIC MODEL IN TERMS OF ERROR WITH RESPECT TO ROAD
Caf = 222685.8 / 2; % Cornering stiffnes az első kerékhez
Car = 136242.8 / 2; % Cornering stiffnes a hátsó kerékhez
lf = 1236e-3; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 2789e-3 - lf; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 2300; % Az autó tömege
Vx = 10; % Az autó haladási sebessége
Iz = 2873; % Az autó tehetetlenségi nyomatéka

%% Contious time state space model 
Av = [0 1                            0                      0;...
     0 -(2*Caf+2*Car)/(m*Vx)        (2*Caf+2*Car)/m        (-2*lf*Caf+2*lr*Car)/(m*Vx);...
     0 0                            0                      1;...
     0 -(2*lf*Caf-2*lr*Car)/(Iz*Vx) (2*lf*Caf-2*lr*Car)/Iz -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*Vx)];

Bv1 = [0; 2*Caf/m; 0; 2*lf*Caf/Iz];
Bv2 = [0; -(2*lf*Caf-2*lr*Car)/(m*Vx)-Vx; 0; -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*Vx)];

% Szűrő
num_wf = [0.02633, 0.0238, 2.286, 0.2335, 0.02902];
den_wf = [1, 2.527, 4.584, 2.993, 1.373];
[Af, Bf, Cf, Df] = tf2ss(num_wf, den_wf);
% v vehicle
% f filter
% Kimenetek kinyerése
Ca  = Av(2, :); % Gyorsulashoz a = Ca*xv + Da1*delta   
Da1 = Bv1(2, :); %gyorsulashoz
Da2 = Bv2(2, :);
C = [1 0 0 0]; %eredeti poz hiba kimenete

% Kiterjesztett Állapotmátrixok
Ac = [Av        zeros(4, 4); ...
      Bf * Ca   Af];
Bc1 = [Bv1;...
       Bf * Da1];
Bc2 = [Bv2;...
       zeros(4,1)];
Bc = [Bc1 Bc2];

Cc = [C        zeros(1, 4); ... 
      Df * Ca  Cf];             
Dc = [0 0;...
      Df * Da1 0];

% Folytonosból diszkrét
c_ss = ss(Ac, Bc, Cc, Dc);
Ts = 33e-3; 
d_ss = c2d(c_ss, Ts);

%% MPC Setup
Nc = 4;
Np = 30;
% Q mátrix: 1. elem: e_y súlya, 2. elem: a_f súlya
Q = diag([10, 1]); 
Rw = 100;

Ap = d_ss.A;
Bp1 = d_ss.B(:,1); % Kormányzás bemenet
Bp2 = d_ss.B(:,2); % Zavarás (pályagörbület)
Cp = d_ss.C;

% Feltételezem, hogy az mpcgain a kiterjesztett [dx; y] formára hozza a mátrixokat
[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e, F, Phi] = mpcgain(Ap, Bp1, Cp, Nc, Np, Q);

%% Receding Horizon control init
t_vec = 0:Ts:35; 
N_sim = length(t_vec);
n_in = 1; % Bemenetek szama
m1 = 2;   % Kimenetek szama (e_y és a_f)
n = 8;    % Allapotok szama

xm = zeros(n, 1); 
y = Cp * xm;
u = 0;
Xf = [zeros(n, 1); y]; % Kiterjesztett állapot: 10x1 méret

Y = zeros(N_sim, m1);
U = zeros(N_sim, n_in);
delta_x = zeros(N_sim, n);
DU = zeros(N_sim, n_in);

R_mat = kron(eye(Nc), Rw);
% Zavarás vektor (kívánt szögsebesség a pályából)
r = generate8likePath(t_vec, 10/20); 

%% Constraints
u_max = 0.5; u_min = -0.5;
du_max = 0.08; du_min = -0.08;

Umax = ones(Nc, 1) * u_max; Umin = ones(Nc, 1) * u_min;
dUmax = ones(Nc, 1) * du_max; dUmin = ones(Nc, 1) * du_min;

% Kimeneti korlátok: [Max pozíció hiba; Max szűrt gyorsulás]
y_max = [1.0; 100000]; 
y_min = [-1.0; -100000]; 

Ymax = repmat(y_max, Np, 1);
Ymin = repmat(y_min, Np, 1);

I_n_in = eye(n_in);
C1 = repmat(I_n_in, Nc, 1);
C2 = kron(tril(ones(Nc)), I_n_in);

M1 = [C2; -C2];
M2 = [eye(Nc*n_in); -eye(Nc*n_in)];
M3 = [-Phi; Phi];
% JAVÍTVA: M3 hozzáadása a korlátmátrixhoz!
M = [M1; M2; M3]; 

lambda0 = zeros(size(M, 1), 1); % Lagrange szorzók inicializálása
H = (Phi_Phi + R_mat);
last_ucon = 0;


A_raw = zeros(N_sim, 1); % Tároló a nem szűrt gyorsulásnak
%% Receding Horizon Control Loop
for kk = 1:N_sim
    
    % Korlátvektorok frissítése az aktuális állapottal
    N1 = [Umax - C1*last_ucon;...
         -Umin + C1*last_ucon];
    N2 = [dUmax; -dUmin];
    N3 = [-Ymin + F*Xf;... 
           Ymax - F*Xf];
    
    % JAVÍTVA: N3 hozzáadása a gamma vektorhoz!
    gamma = [N1; N2; N3];
    
    % A hiba minimalizálása. A referenciánk e_y = 0 és a_f = 0, így Phi_R kiesik (0).
    deltaU_unconstrained = - H \ (Phi_F * Xf);
    
    f = Phi_F * Xf;
    % QPhild függvény meghívása a korlátos megoldáshoz
    deltaU = QPhild(H, f, M, gamma, lambda0);
   
    % Csak az első vezérlő jelet alkalmazzuk
    delta_u_k = deltaU(1:n_in);
    DU(kk, :) = delta_u_k;
    
    % Új u kiszámítása
    u = last_ucon + delta_u_k;
  
    % Alkalmazás a növényen (PLANT) - Beleértve a zavarást (r)
    xm_old = xm;
    xm = Ap * xm + Bp1 * u + Bp2 * r(kk);
    y = Cp * xm;
    
    A_raw(kk) = Ca * xm(1:4) + Da1 * u + Da2 * r(kk);
    % Kiterjesztett állapot frissítése a következő iterációhoz: Xf = [Δx; y]
    Xf = [(xm - xm_old); y];
    
    % Adatok mentése
    Y(kk, :) = y';
    U(kk) = u;
    delta_x(kk,:) = (xm - xm_old)';
    last_ucon = u;
end

%% Plotting Results 
k_time = t_vec; % Használjuk az igazi idővektort

figure('Name', 'MPC Járműdinamika Analízis', 'Position', [100 100 1200 800]); 

% 1. ábra: Beavatkozó jel (u)
subplot(3, 2, 1)
stairs(k_time, U, 'k', 'LineWidth', 1.5); hold on;
plot(k_time, ones(size(k_time))*u_max, 'r--', 'LineWidth', 1);
plot(k_time, ones(size(k_time))*u_min, 'r--', 'LineWidth', 1);
grid on; ylabel('\delta_f [rad]'); title('Kormányszög (u) és korlátok');

% 2. ábra: Delta u
subplot(3, 2, 2)
stairs(k_time, DU, 'b', 'LineWidth', 1.2); hold on;
plot(k_time, ones(size(k_time))*du_max, 'r:', 'LineWidth', 1.5);
plot(k_time, ones(size(k_time))*du_min, 'r:', 'LineWidth', 1.5);
grid on; ylabel('\Delta u [rad/s]'); title('Kormányzás változási sebessége');

% 3. ábra: Pozíció hiba (e_y)
subplot(3, 2, 3)
plot(k_time, Y(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(k_time, ones(size(k_time))*y_max(1), 'r--', 'LineWidth', 1);
plot(k_time, ones(size(k_time))*y_min(1), 'r--', 'LineWidth', 1);
grid on; ylabel('e_y [m]'); title('Oldalirányú pozícióhiba (1. kimenet)');

% 4. ábra: Szűrt vs Nyers gyorsulás összehasonlítása
subplot(3, 2, 4)
plot(k_time, A_raw, 'k:', 'LineWidth', 1); hold on; % Nyers gyorsulás pontozva
plot(k_time, Y(:, 2), 'm-', 'LineWidth', 2);      % Szűrt gyorsulás folytonos vonallal
grid on; 
ylabel('Gyorsulás [m/s^2]'); 
title('ISO 2631-1 Szűrt vs Nyers Gyorsulás');
legend('Nyers (a_{raw})', 'Szűrt (a_f)');

% 5. ábra: Néhány fontosabb járműállapot (Delta X)
subplot(3, 2, 5); 
plot(k_time, delta_x(:,1), 'g', 'LineWidth', 1.5); hold on;
plot(k_time, delta_x(:,3), 'c', 'LineWidth', 1.5);
grid on; legend('\Delta e_y', '\Delta e_\psi'); title('Pozíció és irányszög hiba változása');

subplot(3, 2, 6); 
plot(k_time, delta_x(:,2), 'm', 'LineWidth', 1.5); hold on;
plot(k_time, delta_x(:,4), 'r', 'LineWidth', 1.5);
grid on; legend('\Delta dy', '\Delta d\psi'); title('Sebességek változása');

sgtitle('MPC Járműirányítás ISO 2631-1 komfortszűrővel');