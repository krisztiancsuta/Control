function [eta, km, lambda, du_con] = QPhild(H, f, A_cons, b, lambda0, use_nn_prediction, Xf, r, last_ucon, kk)
    % Inputs:
    % - H: Hessian matrix
    % - f: Gradient vector
    % - A_cons: Constraint matrix
    % - b: Constraint vector
    % - lambda0: Initial Lagrange multipliers
    % - use_nn_prediction: Flag to use neural network for initial guess
    % - Xf: State feedback vector
    % - r: Reference signal
    % - last_ucon: Last control input
    %
    % Outputs:
    % - eta: Optimized solution
    % - km: Number of iterations
    % - lambda: Lagrange multipliers

    [n1, ~] = size(A_cons);
    
    eta = -H \ f;  
    
    du_con = 0;
    lambda = zeros(size(lambda0));

    kk = 0;
    for i = 1:n1
        if (A_cons(i,:) * eta > b(i))
            kk = kk + 1;
            break;
        end
    end
    
    km = 1;
    if kk ~= 0
        P = A_cons * (H \ A_cons');
        d = (A_cons * (H \ f) + b);
        [n, ~] = size(d);
        lambda = lambda0;
        
        for km = 1:4000
            lambda_p = lambda;
            for i = 1:n
                w = P(i,:) * lambda - P(i,i) * lambda(i,1);
                w = w + d(i,1);
                la = -w / P(i,i);
                lambda(i,1) = max(0, la);
            end
            
            al = (lambda - lambda_p)' * (lambda - lambda_p);
            if al < 1e-12
                break;
            end
        end
        
        eta = -H \ f - H \ A_cons' * lambda;

        du_con = - H \ A_cons' * lambda;
        du_con = du_con(1);
    end
     
end
