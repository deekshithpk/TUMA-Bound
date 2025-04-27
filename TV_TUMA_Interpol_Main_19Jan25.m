tic
clear variables
clc

% Step 1: Initialize variables
k = 128;
n = 38400;
M = 2^k;
Ka = 30;
Ma = 25; % Fixed Ma value
EbN0_min =1.3; % Minimum EbN0 in dB
EbN0_max =1.5; % Maximum EbN0 in dB
num_EbN0_vals =5 ; % Number of EbN0 values
EbN0_vec = linspace(EbN0_min, EbN0_max, num_EbN0_vals); % Predefined EbN0 vector
eps_vec = zeros(1, num_EbN0_vals); % To store epsilon values
eps_target = 0.05; % Target epsilon

target_closeness = 0.01; % Define closeness threshold for eps_target

R = log(M) / n;
num_discr_rho = 100;
Zero_prox = 0;
P_frac =0;
num_discr_P = 100;

% Generate the multiplicity vector and Ka_new
[n_vec, Ka_new] = generateMultiplicityVector(Ma, Ka);


% Compute TV distance between scaled multiplicity and Uniform(Ka)
P = n_vec / Ka_new;                % Scaled multiplicity vector
Q = ones(1, Ka_new) / Ka_new;          % Uniform distribution
overlap_distance = 0.5 * sum(abs(P - Q(1:length(P))));
remaining_distance = 0.5 * sum(Q(length(P) + 1:end));
TV_distance = overlap_distance + remaining_distance


% Dynamically load the appropriate .mat file based on the value of Ma
filename = sprintf('Subset_Data_Ma_%d_Ka_%d_25Jan2025.mat', Ma, Ka);
if exist(filename, 'file')
    data = load(filename);
    fieldname = fieldnames(data);
    subset_data = data.(fieldname{1}); % Extract the data structure
    
    % Extract the required fields into vectors
    n_S_vec = arrayfun(@(x) x.n_S, subset_data);
    Length_S_vec = arrayfun(@(x) x.Length_S, subset_data);
   
    
    % Determine the maximum lengths for padding
max_n_N_length = max(arrayfun(@(x) size(x.n_N_ell1_norm_vec, 1), subset_data));
max_Length_N_length = max(arrayfun(@(x) size(x.SubsetLengths, 1), subset_data));

% Initialize matrices with zeros
n_N_mat = zeros(length(subset_data), max_n_N_length);
Length_N_mat = zeros(length(subset_data), max_Length_N_length);

% Populate matrices row by row
for i = 1:length(subset_data)
    % Extract vectors for the current row
    n_N_current = subset_data(i).n_N_ell1_norm_vec(:, 1)';
    Length_N_current = subset_data(i).SubsetLengths(:, 1)';
    
    % Assign to the respective row in the matrices
    n_N_mat(i, 1:length(n_N_current)) = n_N_current;
    Length_N_mat(i, 1:length(Length_N_current)) = Length_N_current;
end


% Determine the maximum length for padding
max_ComplementSet_length = max(arrayfun(@(x) length(x.ComplementSet), subset_data));

% Initialize the ComplementSet matrix with zeros
ComplementSet_mat = zeros(length(subset_data), max_ComplementSet_length);

% Populate the ComplementSet matrix row by row
for i = 1:length(subset_data)
    % Extract the current ComplementSet
    ComplementSet_current = subset_data(i).ComplementSet;
    
    % Assign to the respective row in the matrix
    ComplementSet_mat(i, 1:length(ComplementSet_current)) = ComplementSet_current;
end




    Multiplicity_vec = arrayfun(@(x) x.Multiplicity, subset_data);
else
    error('File %s does not exist.', filename);
end



% Step 2: Loop over predefined EbN0 values
for j = 1:num_EbN0_vals
    EbN0_dB = EbN0_vec(j);

    % Compute P_input for this EbN0
    P_input = 2 * k * 10^(EbN0_dB / 10) / n;

    % Compute epsilon using TUMA_Bound_OptP_9Jan25
    [~, eps_vec(j), ~] = TUMA_Bound_OptP_19Jan25_Vera(k, n, Ma, Ka, R, ...
                                                   num_discr_rho, P_input, ...
                                                   P_frac, num_discr_P, ...
                                                   Zero_prox, n_vec, n_S_vec, Length_S_vec, n_N_mat, Length_N_mat, Multiplicity_vec, ComplementSet_mat);
end

% Step 3: Analyze epsilon values
valid_indices = find(eps_vec > 0 & eps_vec < 1);
if isempty(valid_indices)
    if min(eps_vec) > eps_target
        disp('All epsilon values are greater than eps_target. Consider increasing the EbN0 range.');
    elseif max(eps_vec) < eps_target
        disp('All epsilon values are less than eps_target. Consider decreasing the EbN0 range.');
    else
        disp('No valid epsilon values in the current range. Check the input parameters.');
    end
else
    refined_EbN0_vec = EbN0_vec(valid_indices);
    refined_eps_vec = eps_vec(valid_indices);
    % Display the refined EbN0 and epsilon vectors
disp('Refined EbN0 values:');
disp(refined_EbN0_vec);

disp('Refined epsilon values:');
disp(refined_eps_vec);

Ma
Ka
Ka_new
    % Find closest values to eps_target
    close_indices = find(abs(refined_eps_vec - eps_target) < target_closeness);
    if length(close_indices) <= 2
        disp('Interpolation is done from a few values. Consider refining the EbN0 range.');
    end

    % Interpolation to estimate EbN0 for eps_target
    EbN0_est = interp1(refined_eps_vec, refined_EbN0_vec, eps_target, 'linear', 'extrap');
    fprintf('Estimated EbN0 for eps_target = %.5f: %.2f dB\n', eps_target, EbN0_est);
end

% Step 4: Plot epsilon versus EbN0
figure;
semilogy(EbN0_vec, eps_vec, '-o');
xlabel('EbN0 (dB)');
ylabel('Epsilon');
title(sprintf('Epsilon vs EbN0 for Ma = %d, TV Distance = %.3f', Ma, TV_distance));
grid on;
TV_distance
% disp('Results:');
% disp('EbN0 (dB) values:'); disp(EbN0_vec);
% disp('Epsilon values:'); disp(eps_vec);

toc
