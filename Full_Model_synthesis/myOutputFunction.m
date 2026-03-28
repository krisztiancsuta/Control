function y = myOutputFunction(x, u, psi_dot_des)
% MYOUTPUTFUNCTION kiszámítja a járműmodell kimeneteit.
% Bemenetek:
%   x = [e1; de1; e2; de2; vx] (5x1 állapotvektor)
%   u = [steering_angle; force] (2x1 bemeneti vektor)
%
% Kimenet:
%   y = [e1; vx] (2x1 kimeneti vektor)

    % Memória foglalása a kimenetnek (2 kimenetünk van)
    y = zeros(2,1);
    
    % 1. kimenet: Keresztirányú pozícióhiba (lateral position error)
    y(1) = x(1);
    
    % 2. kimenet: Hosszirányú sebesség (longitudinal velocity)
    y(2) = x(5);
    
end