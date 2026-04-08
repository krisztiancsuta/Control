function [C, D] = myOutputJacobian(x, u, psi_dot_des)
% MYOUTPUTJACOBIAN kiszámítja a kimeneti függvény analitikus deriváltjait.
% Bemenetek:
%   x = [e1; de1; e2; de2; vx]
%   u = [steering_angle; force]
%
% Kimenetek:
%   C = dy/dx (2x5 mátrix)
%   D = dy/du (2x2 mátrix)

    % C mátrix: A kimenetek deriváltjai az állapotokra nézve
    % y_1 = x_1 -> deriváltja x_1 szerint 1, a többi állapot szerint 0
    % y_2 = x_5 -> deriváltja x_5 szerint 1, a többi állapot szerint 0
    C = [1, 0, 0, 0, 0;
         0, 0, 0, 0, 1];

    % D mátrix: A kimenetek deriváltjai a beavatkozó jelekre nézve
    % A kimenetünk közvetlenül nem függ a bemenettől (nincs feedthrough)
    D = [0, 0;
         0, 0];
end