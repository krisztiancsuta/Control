% Paraméterek (4A szűrő aszimmetrikus esete)
L = 0.75e-3; % 1.5mH / 2 [cite: 752]
C = 4.4e-9;  % 2.2nF * 2 [cite: 752]

% Átviteli függvény létrehozása: H(s) = 1 / (L*C*s^2 + 0*s + 1)
num = [1];
den = [L*C 0 1];
H = tf(num, den);

% Karakterisztika kirajzolása (Bode-diagram)
bode(H);
grid on;
title('Ideális LC szűrő átviteli függvénye');