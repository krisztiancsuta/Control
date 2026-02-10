function [phides, u_delta_ff] = generate8likePath(t, max_yaw_rate, m, lf, lr, Caf, Car, Vx,K3)
    % t: idővektor
    % max_yaw_rate: maximális legyezési sebesség (Vx/R)
    % m, lf, lr, Caf, Car, Vx: Jármű paraméterek a feedforward számításhoz
    
    phides = zeros(size(t)); 
    
    % Idővektor felosztása szakaszokra (egyenes, kanyar, egyenes, kanyar)
    num_segments = 4;  
    segment_length = floor(length(t) / num_segments);  
    
    % 1. Egyenes szakasz (yaw rate = 0)
    phides(1:segment_length) = 0;  
    
    % 2. Első kanyar (pozitív yaw rate)
    phides(segment_length+1 : 2*segment_length) = max_yaw_rate;  
    
    % 3. Második egyenes szakasz (yaw rate = 0)
    phides(2*segment_length+1 : 3*segment_length) = 0;  
    
    % 4. Második kanyar (negatív yaw rate - ellenkező irány)
    phides(3*segment_length+1 : end) = -max_yaw_rate;  

%% Feedforward kormányszög számítása (Végtelen R kezelése)
    L = lf + lr;
    
    % Görbület (1/R) kiszámítása: Ha phides=0, a görbület is 0 (egyenes)
    % Nincs osztás R-rel, így elkerüljük az Inf/NaN hibákat
    invR = phides ./ Vx; 
    
    % Alulkormányzottsági gradiens (Rajamani 3.14)
    Kv = (m * lr) / (2 * Caf * L) - (m * lf) / (2 * Car * L);
    
    % Oldalgyorsulás (ay = Vx^2 / R)
    ay = Vx^2 .* invR;
    
    % Állandósult állapotú orientációs hiba (Rajamani 3.15 / 41. oldal környéke)
    % e2_ss = hiba, ami a kanyar íve miatt elkerülhetetlen feedforward nélkül
    e2_ss = invR .* (-lr + (lf * m * Vx^2) / (2 * Car * L));
    
    % Teljes feedforward kormányszög (delta_ff)
    % delta_ff = L/R + Kv*ay - K3*e2_ss  (Figyelem: K3*e2_ss levonásra kerül a szabályzó kompenzálásához)
    u_delta_ff = (L .* invR) + (Kv .* ay) + (K3 .* e2_ss);
end