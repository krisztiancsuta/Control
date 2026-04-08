% ISO 2631-1 Wf (Motion Sickness) szűrő közelítéseinek összehasonlítása
clear; clc;close all;
%% 1. Diszkrét adatok kinyerése a táblázatból (image_9b450b.jpg)
% Frekvenciák [Hz]
f_std = [0.02, 0.025, 0.0315, 0.04, 0.05, 0.063, 0.08, 0.1, 0.125, 0.16, ...
         0.2, 0.25, 0.315, 0.4, 0.5, 0.63, 0.8, 1, 1.25, 1.6, 2, 2.5, 3.15, 4];

% Wf faktorok (faktor / 1000 = magnitúdó)
Wf_factor = [24.2, 37.7, 59.7, 97.1, 157, 267, 461, 695, 895, 1006, ...
             992, 854, 619, 384, 224, 116, 53.0, 23.5, 9.98, 3.77, 1.55, 0.64, 0.25, 0.097];
Wf_std_mag = Wf_factor / 1000;

%% 2. Közelítő átviteli függvények definiálása (Zuo & Nayfeh )
s = tf('s');

Wf2 = tf([0.8892, 0], [1, 0.8263, 1.163]);
Wf3 = tf([0.05726, 0, 3.876, 0], [1, 4.263, 4.777, 4.396]);
Wf4 = tf([0.02633, 0.0238, 2.286, 0.2335, 0.02902], [1, 2.527, 4.584, 2.993, 1.373]);
Wf5 = tf([0.1457, 0.2331, 13.75, 1.705, 0.3596], [1, 7.757, 19.06, 28.37, 18.52, 7.230]);

%% 3. Ábrázolás
f_range = logspace(-2, 1, 1000); 
w_range = 2 * pi * f_range;

figure('Color', 'w', 'Position', [100, 100, 900, 600]);

% Közelítő görbék kiszámítása
[m2, ~] = bode(Wf2, w_range);
[m3, ~] = bode(Wf3, w_range);
[m4, ~] = bode(Wf4, w_range);
[m5, ~] = bode(Wf5, w_range);

% Folytonos görbék rajzolása
semilogx(f_range, squeeze(m2), ':', 'LineWidth', 1.2, 'Color', [0.7 0.7 0.7]); hold on;
semilogx(f_range, squeeze(m3), '--', 'LineWidth', 1.2, 'Color', 'b');
semilogx(f_range, squeeze(m4), '-', 'LineWidth', 2, 'Color', 'r');
semilogx(f_range, squeeze(m5), '-.', 'LineWidth', 1.2, 'Color', 'k');

% Szabványos pontok (diszkrét értékek a táblázatból)
plot(f_std, Wf_std_mag, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);

% Formázás
grid on;
xlabel('Frekvencia (Hz)');
ylabel('Magnitúdó (Súlyozó tényező)');
title('W_f (Motion Sickness) szűrő: Elméleti közelítések vs. ISO 2631-1 pontok');
legend('2. rendű közelítés', '3. rendű közelítés', '4. rendű közelítés', ...
       '5. rendű közelítés', 'ISO 2631-1 szabványos pontok', 'Location', 'northeast');
axis([0.01 10 0 1.2]);