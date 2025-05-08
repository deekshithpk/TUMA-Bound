function [eps_sum_nN_j_opt_rho,min_optimal_rho_value] = computeTUMABoundCore(t, ell, n_S_ell_1norm, eye, n_N_eye_1norm_vec, jay_vec, rho, P1, n, Ma, Ka, R)
    
% -------------------------------------------------------------------------
% Function: computeTUMABoundCore
%
% Purpose:
%   Computes the core expression of the TUMA bound for a given value of rho 
%   and fixed transmit power P1. This function performs the main calculation 
%   contributing to the bound value and is called internally by 
%   computeTUMABoundAtP.
%
% Inputs:
%   - t                 : Ranges upto Ka (a possible distance between transmitted and decoded multiplicity) 
%   - ell               : Support set size of impostor set
%   - n_S_ell_1norm     : Sum multiplicity of miss-detected set
%   - eye               : Support size of inflated set
%   - n_N_eye_1norm_vec : Sum of transmitted multiplicity of inflated set
%   - jay_vec           : Sum multiplicty of impostor set
%   - rho               : Gallager rho-trick parameter
%   - P1                : Transmit power (codebook variance)
%   - n                 : Block length
%   - Ma, Ka            : Number of active messages and users
%   - R                 : Transmission rate
%
% Outputs:
%   - eps_sum_nN_j_opt_rho      : Evaluated bound contribution when summed
%   over N, j for optimized rho
%   - min_optimal_rho_value     : Minimum rho that minimizes the bound 
%
% Notes:
%   This function contains the core logic for bound evaluation. 
% -------------------------------------------------------------------------




    % Create grids for n_N_eye_1norm_vec, jay_vec, and rho
    [n_N_eye_1norm_grid, jay_grid, rho_grid] = ndgrid(n_N_eye_1norm_vec, jay_vec, rho);

    % Initialize a combined validity mask for all dimensions
    valid_mask = true(size(n_N_eye_1norm_grid));

    % Condition 1: n_N_eye_1norm >= t + eye - n_S_ell_1norm
    condition1 = n_N_eye_1norm_grid >= (t + eye - n_S_ell_1norm);
    valid_mask = valid_mask & condition1;

    % Apply the validity mask to each grid
    n_N_eye_1norm_grid(~valid_mask) = NaN;
    jay_grid(~valid_mask) = NaN;
    rho_grid(~valid_mask) = NaN;

    % Compute c_min_values
    c_min_values = computeCminParameter(t, ell, eye, jay_grid, Ma, n_S_ell_1norm);


    % Compute lambda_values with broadcasting
    lambda_values = ((P1 * c_min_values - 2) + sqrt((P1 * c_min_values - 2).^2 + 4 .* (P1 * c_min_values .* (1 + rho_grid)))) ./ ...
                    (4 * (P1 * c_min_values) .* (1 + rho_grid));
    lambda_values(isnan(lambda_values)) = 0;

    % Define b_function and calculate b_values
    b_function = @(lambda_vals, c_min_vals, P1) ...
        lambda_vals - lambda_vals ./ (1 + 2 * P1 * c_min_vals .* lambda_vals);
    b_values = b_function(lambda_values, c_min_values, P1);
    b_values(isnan(b_values)) = 0;

    % Define a_function and calculate a_values
    a_function = @(lambda_vals, c_min_vals, P1) ...
        0.5 * log(1 + 2 * P1 * c_min_vals .* lambda_vals);
    a_values = a_function(lambda_values, c_min_values, P1);
    a_values(isnan(a_values)) = 0;

    % Define E_function and calculate E_values
    E_function = @(a_vals, b_vals, rho_vals) rho_vals .* a_vals + 0.5 * log(1 - 2 .* b_vals .* rho_vals);
    E_values = E_function(a_values, b_values, rho_grid);
    E_values(isnan(E_values)) = Inf;

    % Compute intermediate terms for R2
    Ma_minus_ell_minus_i = Ma - ell - eye;
    log_factorial_Ma_minus_ell_minus_i_minus_1 = gammaln(Ma_minus_ell_minus_i);
    
    if eye > 1
    log_factorial_i_minus_1 = gammaln(eye - 1);
    else
     log_factorial_i_minus_1=0;   
    end

    log_factorial_Ma_minus_ell_minus_i = gammaln(Ma_minus_ell_minus_i + 1);

    R2_term_a1 = (gammaln(Ma+1) - gammaln(ell + 1) - log_factorial_Ma_minus_ell_minus_i) / n;
    R2_term_a2 = (gammaln(Ma-ell+1) - gammaln(eye + 1) - gammaln(Ma-ell-eye+1)) / n;
    R2_term_a = 0;%R2_term_a1 + R2_term_a2;

    % Compute t_minus_nS_plus_i and handle invalid cases
t_minus_nS_plus_i = t - n_S_ell_1norm + eye;
log_factorial_t_minus_nS_plus_i_minus_1 = NaN; % Initialize with NaN
if t_minus_nS_plus_i > 0
    log_factorial_t_minus_nS_plus_i_minus_1 = gammaln(t_minus_nS_plus_i);
end
    log_factorial_t_minus_nS = gammaln(t - n_S_ell_1norm + 1);
    
    if eye<=1
    R2_term_c = 0;
    else
    R2_term_c = (log_factorial_t_minus_nS_plus_i_minus_1 - log_factorial_i_minus_1 - log_factorial_t_minus_nS) / n;
    end

    % Compute nS_plus_nN_minus_t and handle invalid cases
nS_plus_nN_minus_t = n_S_ell_1norm + n_N_eye_1norm_grid - t;
log_factorial_nS_plus_nN_minus_t_minus_1 = NaN(size(nS_plus_nN_minus_t)); % Initialize with NaN
valid_nS_plus_nN_minus_t = nS_plus_nN_minus_t > 0; % Valid cases where nS_plus_nN_minus_t > 0
log_factorial_nS_plus_nN_minus_t_minus_1(valid_nS_plus_nN_minus_t) = ...
gammaln(nS_plus_nN_minus_t(valid_nS_plus_nN_minus_t));

% Compute nS_plus_nN_minus_t_minus_eye_plus_1 and handle invalid cases
nS_plus_nN_minus_t_minus_eye_plus_1 = nS_plus_nN_minus_t - eye + 1;
log_factorial_nS_plus_nN_minus_t_minus_i = NaN(size(nS_plus_nN_minus_t)); % Initialize with NaN
valid_nS_plus_nN_minus_t_minus_eye_plus_1 = nS_plus_nN_minus_t_minus_eye_plus_1 > 0; % Valid cases where > 0
log_factorial_nS_plus_nN_minus_t_minus_i(valid_nS_plus_nN_minus_t_minus_eye_plus_1) = ...
gammaln(nS_plus_nN_minus_t_minus_eye_plus_1(valid_nS_plus_nN_minus_t_minus_eye_plus_1));

if eye<=1
R2_term_b = 0;
else
R2_term_b = (log_factorial_nS_plus_nN_minus_t_minus_1 - log_factorial_i_minus_1 - log_factorial_nS_plus_nN_minus_t_minus_i) / n;
end

    R2_final = R2_term_a + min(R2_term_b, R2_term_c);

    


    % Compute intermediate terms for R1
% Initialize with NaN
log_factorial_j_minus_1 = NaN(size(jay_grid)); 
log_factorial_j_minus_ell = NaN(size(jay_grid));
log_factorial_ell_minus_1 = NaN(size(jay_grid)); 

% Safe computations with conditions
valid_j_minus_1 = jay_grid > 0; % jay_grid must be > 0
log_factorial_j_minus_1(valid_j_minus_1) = gammaln(jay_grid(valid_j_minus_1));

valid_j_minus_ell = (jay_grid - ell + 1) > 0; % jay_grid - ell + 1 must be > 0
log_factorial_j_minus_ell(valid_j_minus_ell) = gammaln(jay_grid(valid_j_minus_ell) - ell + 1);

% Handle ell separately as it is scalar
if ell > 0
    log_factorial_ell_minus_1(:) = gammaln(ell); % Assign the same value across valid indices
else
    log_factorial_ell_minus_1(:) = NaN; % Invalid case
end


    R1_term_b = (log_factorial_j_minus_1 - log_factorial_ell_minus_1 - log_factorial_j_minus_ell) / n;

    % Compute t_minus_j and handle invalid cases
t_minus_j = t - jay_grid;
log_factorial_t_minus_j_minus_1 = NaN(size(t_minus_j)); % Initialize with NaN
valid_t_minus_j = t_minus_j > 0; % Valid cases where t_minus_j > 0
log_factorial_t_minus_j_minus_1(valid_t_minus_j) = gammaln(t_minus_j(valid_t_minus_j));

% Compute t_minus_j - Ma_minus_ell_minus_i + 1 and handle invalid cases
t_minus_j_minus_Ma_minus_ell_minus_i_plus_1 = t_minus_j - Ma_minus_ell_minus_i + 1;
log_factorial_t_minus_j_minus_Ma_minus_ell_minus_i = NaN(size(t_minus_j_minus_Ma_minus_ell_minus_i_plus_1)); % Initialize with NaN
valid_t_minus_j_minus_Ma_minus_ell_minus_i_plus_1 = t_minus_j_minus_Ma_minus_ell_minus_i_plus_1 > 0; % Valid cases
log_factorial_t_minus_j_minus_Ma_minus_ell_minus_i(valid_t_minus_j_minus_Ma_minus_ell_minus_i_plus_1) = ...
    gammaln(t_minus_j_minus_Ma_minus_ell_minus_i_plus_1(valid_t_minus_j_minus_Ma_minus_ell_minus_i_plus_1));

% Compute R1_term_c with updated values
R1_term_c = (log_factorial_t_minus_j_minus_1 - log_factorial_Ma_minus_ell_minus_i_minus_1 - log_factorial_t_minus_j_minus_Ma_minus_ell_minus_i) / n;

    R1_term_a = ell * R - gammaln(ell + 1) / n;
    R1_final = R1_term_a + R1_term_b + R1_term_c;

    %R_max=(gammaln(Ka)-gammaln(Ka-Ma+1)-gammaln(Ma))/n;

    % Combine R1 and R2 for Rate Values
    Rate_values = R1_final + R2_final;


    Rate_values(isnan(Rate_values)) = -Inf;


    % Calculate f_values
    f_values = exp(n * (rho_grid .* Rate_values - E_values));
    f_values(isnan(f_values)) = 0;

    % Take the minimum of f_values over the rho dimension
min_f_values_rho = min(f_values, [], 3);

% Compute the minimum values and optimal indices along the rho dimension
[min_f_values_rho, optimal_rho_indices] = min(f_values, [], 3);

% Extract the rho values corresponding to the optimal indices
optimal_rho_values = rho(optimal_rho_indices);

% Save the minimum of the optimal rho values
min_optimal_rho_value = min(optimal_rho_values(:));

% Sum over the j dimension (jay_grid)
sum_f_values_j = sum(min_f_values_rho, 2);

% Sum over the n_N_eye_1norm dimension
final_sum = sum(sum_f_values_j, 1);

% Assign final results
eps_sum_nN_j_opt_rho = final_sum;


    
end
