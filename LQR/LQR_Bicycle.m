%% Allapotok definialasa
% Dynamics of lateral vehicle motion
% Nagyobb sebességek esetén már nem lehet azt állítani hogy mindegyik
% kerék sebességének iránya megegyezik a kerék állásának irányával. Ezért
% dinamkius modellt kell alkalmazni. 
% 2 DOF a modell: A jármű laterális pozíciója y
% (a jármű forgásközéppontjához viszonyitva),
% a jármű globális koordinátarendszerben X tengelyhez mért haladási szöge
% phi.
% Az ay gyorsulás a acp és y'' gyorsulás összegéből adódik
% F = m*a alkalmzásával és a nyomatékokra felírható egyenletekből
% Az Fkerék laterális erők arányos a kerekek !kicsi! slip-anglejevel

% A példa értékeke egy toyota corolla E210 hez tartoznak 
Caf = 11600; % Cornering stiffnes az első kerékhez
Car = 14200; % Cornering stiffnes a hátsó kerékhez
lf = 1.1; % Az autó tömegközéppontjától mért elülső tengelytáv
lr = 1.5; % Az autó tömegközéppontjától mért hátsó tengelytáv
m = 1250; % Az autó tömege
Vx = 30; % Az autó haladási sebessége a saját koordinátarendszerében
Iz = 3000; % Az autó tehetetlenségi nyomatéka
% Vy = y', ay = y''
%% Az állapotvektorunk y, y', phi, phi' == 
% laterális pozíció,
% laterális sebesség,
% legyezési szög,
% legyezési szögsebesség
% A bemenet u = delta => Első kerék kormányszöge

A = [0 1                            0  0;...
     0 -(2*Caf+2*Car)/(m*Vx)        0 -Vx-(2*lf*Caf-2*lr*Car)/(m*Vx);...
     0 0                            0  1;...
     0 -(2*lf*Caf-2*lr*Car)/(Iz*Vx) 0 -(2*lf*lf*Caf+2*lr*lr*Car)/(Iz*Vx)];

A = A + 0.001;

B = [0;...
     2*Caf/m;...
     0;...
     2*lf*Caf/Iz];

C = eye(4);

sys = ss(A,B,C,0);
%% Iranyithatosag vizsgalata
Co = ctrb(sys.A,sys.B);
rang = rank(Co);
n = size(A,1);
isControllable = (n == rang)
% Ha iranyithato akkor biztos van optimalis LQ controller hozza
%% LQR súlyozó mátrixai 
R = 1 ;
Q = diag([1 0 1 0]);


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

t  = 0:0.001:4;
dt = t(2) - t(1);

state_names = {'y','vy','phi','vphi'};

% --- 1) Impulzusválasz (x(0)=0) ---
u_imp = zeros(size(t));
u_imp(1) = 1/dt;                 % egység-területű impulzus közelítés
x0_zero = zeros(4,1);

[x_imp, t_out] = lsim(syscl_x,u_imp,t);

% --- 2) Kezdeti állapot válasz (u=0, x4(0)=200) ---
u_zero = zeros(size(t));
x0_ic  = [0; 0; 45; 0];

[x_ic, t_ic] = lsim(syscl_x, u_zero, t, x0_ic);

% --- Plot: 1 ablak, 5 panel ---
figure;
tl = tiledlayout(3,2, "TileSpacing","compact", "Padding","compact");

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
zeros(size(t));
legend('x_1','x_2','x_3','x_4','Location','best');
