function drawpath(t, x_states_ff, x_states_no, Vx, vphides_t)
    dt = t(2) - t(1);
    N = length(t);
    
    % Állapotok kinyerése (e1 az első oszlop)
    e1_ff = x_states_ff(:, 1);
    e1_no = x_states_no(:, 1);
    
    % Tömbök előkészítése
    psi_des = zeros(N, 1);
    X_des = zeros(N, 1);
    Y_des = zeros(N, 1);
    
    % 1. Referencia út integrálása (Frenet-keret alapja)
    for i = 1:N-1
        psi_des(i+1) = psi_des(i) + vphides_t(i) * dt;
        X_des(i+1) = X_des(i) + (Vx * cos(psi_des(i))) * dt;
        Y_des(i+1) = Y_des(i) + (Vx * sin(psi_des(i))) * dt;
    end
    
    % 2. Globális pozíciók számítása (Rajamani transzformáció)
    % Feedback + Feedforward eset
    X_car_ff = X_des - e1_ff .* sin(psi_des);
    Y_car_ff = Y_des + e1_ff .* cos(psi_des);
    
    % Csak Feedback eset
    X_car_no = X_des - e1_no .* sin(psi_des);
    Y_car_no = Y_des + e1_no .* cos(psi_des);
    
    % 3. Megjelenítés
    figure('Color', 'w', 'Name', 'Globális Pályakövetés Összehasonlítása');
    hold on;
    
    % Referencia út (szaggatott fekete)
    plot(X_des, Y_des, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Referencia út');
    
    % Csak Feedback (piros)
    plot(X_car_no, Y_car_no, 'r', 'LineWidth', 1, 'DisplayName', 'Csak Feedback');
    
    % Feedback + Feedforward (kék)
    plot(X_car_ff, Y_car_ff, 'b', 'LineWidth', 2, 'DisplayName', 'LQR + Feedforward');
    
    axis equal; grid on;
    xlabel('X [m]'); ylabel('Y [m]');
    title('Jármű globális trajektóriája a szinuszos pályán');
    legend('Location', 'best');
    
    % Opcionális: Kezdőpont jelölése
    plot(X_des(1), Y_des(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
end