function [max_lambda, min_lambda] = evaluateLambda(P1, Ka, Ma)

% -------------------------------------------------------------------------
% Function: evaluate_lambda
%
% Purpose:
%   Evaluates the maximum and minimum values of the lambda parameter over a 
%   range of rho, c_min, and P2 values, given a fixed P1 (power), number of 
%   active users (Ka), and active messages (Ma). Lambda is the Chernoff parameter in the TUMA bound.
%   This is used to aid choose the delta parameter that comes up in the
%   bound.
%
% Inputs:
%   - P1 : Fixed transmit power (codebook variance)
%   - Ka : Number of active users
%   - Ma : Number of active messages
%
% Outputs:
%   - max_lambda : Maximum value of lambda over the parameter grid
%   - min_lambda : Minimum value of lambda over the parameter grid
%
% -------------------------------------------------------------------------



    % Define the range of rho (discrete values from 0.01 to 1)
    rho = linspace(0.01, 1, 100); % Adjust number of points as needed

    % Define the range of c_min_values (integer values from 1 to 2Ka^2/Ma)
    c_min_values = 1:floor(2 * Ka^2 / Ma);

    % Define the range of P2 (discrete values from 0 to P1)
    P2_values = linspace(0, P1, 50); % Adjust number of points as needed

    % Initialize lambda_values for all combinations of P2
    lambda_values = zeros(length(c_min_values), length(rho), length(P2_values));

    % Compute lambda values for each combination of c_min_values, rho, and P2
    for k = 1:length(P2_values)
        P2 = P2_values(k);
        for i = 1:length(c_min_values)
            for j = 1:length(rho)
                lambda_values(i, j, k) = ((P1 * c_min_values(i) - 2) + sqrt((P1 * c_min_values(i) - 2)^2 + 4 * (P1 * c_min_values(i) * (1 + rho(j))))) / ...
                                        (4 * (P1 * c_min_values(i)) * (1 + rho(j)));
            end
        end
    end

    % Flatten lambda_values to find global max and min across all P2 values
    lambda = lambda_values(:);
    max_lambda = max(lambda);
    min_lambda = min(lambda);

    
end
