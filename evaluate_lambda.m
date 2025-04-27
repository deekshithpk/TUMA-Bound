function [max_lambda, min_lambda] = evaluate_lambda(P1, Ka, Ma)
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
