function drawpath(t, x_states, Vx, vphides_t)
    dt = t(2) - t(1); % Időlépés (0.01)
    N = length(t);
    
    % Állapotok kinyerése
    e1 = x_states(:, 1);
    
    % Tömbök előkészítése a globális koordinátáknak
    psi_des = zeros(N, 1);
    X_des = zeros(N, 1);
    Y_des = zeros(N, 1);
    
    % Manuális integrálás FOR ciklussal
    for i = 1:N-1
        % 1. Út szögének integrálása
        psi_des(i+1) = psi_des(i) + vphides_t(i) * dt;
        
        % 2. Út középvonalának integrálása (X_des, Y_des)
        X_des(i+1) = X_des(i) + (Vx * cos(psi_des(i))) * dt;
        Y_des(i+1) = Y_des(i) + (Vx * sin(psi_des(i))) * dt;
    end
    
    % 3. Az autó tényleges pozíciója (Rajamani transzformáció)
    % Ez nem integrálás, hanem algebrai eltolás az út normálisa mentén
    X_car = X_des - e1 .* sin(psi_des);
    Y_car = Y_des + e1 .* cos(psi_des);

    % Rajzolás
    figure('Color', 'w', 'Name', 'Globalis palya és az auto pályája');
    plot(X_des, Y_des, 'k--', 'LineWidth', 0.5); hold on;
    
    % Csak pontokkal rajzoljuk az autó útját
    plot(X_car, Y_car, 'b-', 'MarkerSize', 3);

    axis equal; grid on;
    xlabel('X [m]'); ylabel('Y [m]');
    title('Globalis palya és az auto pályája');

end