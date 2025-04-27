function [n_vec, Ka_new] = generateMultiplicityVector(Ma, Ka)
    %Generate Zipf pmf
    N = Ma;
    s = 1;
    zipf_pmf = (1:N).^(-s);
    zipf_pmf = zipf_pmf / sum(zipf_pmf); %Normalize to create a PMF
    
    %Scale by Ka and round to nearest integer
    scaled_values = round(zipf_pmf * Ka);
    
    %Identify zeros and adjust multiplicities
    n_vec = scaled_values;
    zero_indices = find(n_vec == 0);
    excess_indices = find(n_vec > 1);
    
    for loop = 1:length(zero_indices)
        if ~isempty(excess_indices)
            n_vec(zero_indices(loop)) = 1; %Assign 1 to the zero index
            n_vec(excess_indices(1)) = n_vec(excess_indices(1)) - 1; % Reduce excess by 1
            if n_vec(excess_indices(1)) == 1
                excess_indices(1) = []; %Remove if it becomes 1
            end
        else
            n_vec(zero_indices(loop)) = 1; %Just assign 1 if no excess exists
        end
    end
    
    % Step 4: Calculate Ka_new
    Ka_new = sum(n_vec);
end
