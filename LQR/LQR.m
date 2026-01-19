%% Allapotok definialasa
% 5.18 Példa egyenleteit felhaszlva
Za = -1.397;
Zq = 1.0;
Ma = -5.47;
Mq = -3.27;
Zd = -0.124;
Md = -13.2;
u = 400;

A = [Za Zq 0 0;...
     Ma Mq 0 0;...
     0  1  0 0;...
     -u 0  u 0];

B = [Zd;...
     Md;...
     0;...
     0];
C = [1 0 0 0];

sys = ss(A,B,C,0)
%% Iranyithatosag vizsgalata
Co = ctrb(sys.A,sys.B);
rang = rank(Co);
n = size(A,1);
isControllable = (n == rang)
% Ha iranyithato akkor biztos van optimalis LQ controller hozza
%% LQR súlyozó mátrixai 
R = 25;
Q = diag([100 0 0 0.0001]);


%% Ricatti egyenlet megoldasa
% At = transpose(A);
% Bt = transpose(B);
% rinv = inv(r);
%CARE => At*P + P*A-P*b*rinv*Bt*P + Q = 0;

% P megoldas
% L Closed Loop sajatertekei
% K Visszacsatolás erősitései
[P,L,K]=care(A,B,Q,R);

%% Zárt rendszer állapotegyenletei 
syscl = ss(A-B*K,B,C,0)

%% A rendszer analízise.

Acl = A - B*K;
syscl_x = ss(Acl, B, eye(4), zeros(4,1));

t  = 0:0.01:10;
dt = t(2) - t(1);

state_names = {'x_1','x_2','x_3','x_4'};

% --- 1) Impulzusválasz (x(0)=0) ---
u_imp = zeros(size(t));
u_imp(1) = 1/dt;                 % egység-területű impulzus közelítés
x0_zero = zeros(4,1);

[x_imp, t_out] = lsim(syscl_x, u_imp, t, x0_zero);

% --- 2) Kezdeti állapot válasz (u=0, x4(0)=200) ---
u_zero = zeros(size(t));
x0_ic  = [0; 0; 0; 200];

[x_ic, t_ic] = lsim(syscl_x, u_zero, t, x0_ic);

% --- Plot: 1 ablak, 5 panel ---
figure;
tl = tiledlayout(3,2, "TileSpacing","compact", "Padding","compact");
title(tl, 'Zárt kör állapotválaszok: impulzus (x(0)=0) és kezdeti eltérés (x_4(0)=200)');

% Panelek 1-4: impulzus -> külön állapotok
for i = 1:4
    nexttile;
    plot(t_out, x_imp(:,i), 'LineWidth', 1.5);
    grid on;
    xlabel('t [s]');
    ylabel(state_names{i});
    title([state_names{i} ' impulzusválasz (x(0)=0)']);
end

% Panel 5: kezdeti állapot -> összes állapot együtt
nexttile;
plot(t_ic, x_ic, 'LineWidth', 1.5);
grid on;
xlabel('t [s]');
ylabel('Állapotok');
title('Kezdeti állapot válasz (x_4(0)=200 m, u(t)=0)');
legend('x_1','x_2','x_3','x_4','Location','best');
