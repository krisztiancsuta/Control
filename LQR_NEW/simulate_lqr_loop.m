function [res] = simulate_lqr_loop(sys, K_gain, signals)
% SIMULATE_LQR_LOOP - Manuális zárt körű szimuláció (for-loop)
%
% Bemenetek:
%   sys.A, sys.B1, sys.B2, sys.C - Folytonos rendszermátrixok
%   K_gain.K_pd, K_gain.K_i     - Szabályozó erősítések
%   signals.t, signals.dist, signals.u_ff, signals.r - Idősorok

    % --- 1. Előkészítés ---
    dt = signals.t(2) - signals.t(1);
    n = length(signals.t);
    nx = size(sys.A, 1);
    
    % Diszkretizálás a pontos léptetéshez (ZOH)
    % A bemeneti mátrixot összevonjuk [B1, B2] a diszkretizáláshoz
    tmp_sys = c2d(ss(sys.A, [sys.B1, sys.B2], sys.C, 0), 0.033);
    Ad = tmp_sys.A;
    Bd1 = tmp_sys.B(:, 1); % Bemeneti csatorna
    Bd2 = tmp_sys.B(:, 2); % Zavarás csatorna

    % --- 2. Memória foglalás ---
    x = zeros(nx, n);       % Állapotok: [e1; e1_dot; e2; e2_dot]
    z = zeros(1, n);        % Integrátor állapota (z' = r - y)
    u_tot = zeros(1, n);    % Teljes kormányszög
    y = zeros(1, n);        % Kimenet (e1)

    % --- 3. Szimulációs ciklus (A zárt kör) ---
    for k = 1:n-1
        % A kimenet mérése (Feedback alapja)
        y(k) = sys.C * x(:, k);
        
        % Integrátor léptetése (Numerikus integrálás)
        % z(k+1) = z(k) + hiba * dt
        z(k+1) = z(k) + (signals.r(k) - y(k)) * dt;
        
        % SZABÁLYOZÁSI TÖRVÉNY (Closed-Loop + Feedforward)
        u_fb = -K_gain.K_pd * x(:, k) + K_gain.K_i * z(k);
        u_tot(k) = u_fb + signals.u_ff(k); % FF kompenzáció hozzáadása
        
        % RENDSZER DINAMIKA (Zárt kör visszacsatolása)
        % Az állapotfrissítésnél az u_tot-ot használjuk
        x(:, k+1) = Ad * x(:, k) + Bd1 * u_tot(k) + Bd2 * signals.dist(k);
    end
    
    % Utolsó pont számítása
    y(n) = sys.C * x(:, n);
    u_tot(n) = -K_gain.K_pd * x(:, n) + K_gain.K_i * z(n) + signals.u_ff(n);

    % --- 4. Eredmények ---
    res.t = signals.t;
    res.x = x';
    res.y = y';
    res.u = u_tot';
end