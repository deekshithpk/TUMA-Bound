function [P_opt, eps_opt, EbN0_dB_opt] = computeOptimalBoundOverP(k, n, Ma, Ka, R, num_discr_rho, P_input, P_frac, num_discr_P, Zero_prox,n_vec, n_S_vec, Length_S_vec, n_N_mat, Length_N_mat, Multiplicity_vec, ComplementSet_mat)%,n1,id_proxy)

% -------------------------------------------------------------------------
% Function: computeOptimalBoundOverP
%
% Purpose:
%   Computes the bound for a discretized set of power values P in [0, P_bound],
%   and returns the optimal power P_opt that minimizes the bound, along with
%   the corresponding minimal bound value (eps_opt) and the associated EbN0.
%
% Inputs:
%   - k, n       : Message length and block length
%   - Ma, Ka     : Number of active messages and users
%   - R          : Transmission rate
%   - num_discr_rho : Discretization count for rho (internal optimization)
%   - P_input    : Upper bound on power sweep range (P ∈ [0, P_input])
%   - P_frac     : Fraction of P_input for discretization granularity
%   - num_discr_P: Number of points in P sweep
%   - Zero_prox  : Threshold for zero proximity (numerical stability)
%   - n_vec, n_S_vec, Length_S_vec, n_N_mat, Length_N_mat
%   - Multiplicity_vec, ComplementSet_mat : Problem-specific subset structures
%
% Outputs:
%   - P_opt      : Optimal transmit power minimizing the bound
%   - eps_opt    : Minimum bound value achieved
%   - EbN0_dB_opt: Eb/N0 (in dB) corresponding to P_opt
%
% Notes:
%   This function is used to optimize over the power values for codebook for
%   the TUMA bound on total variation distance and is a helper function
%   to main_compute_tv_bound_zipf_type.m
% -------------------------------------------------------------------------



[n_vec, Ka_new] = generateMultiplicityVector(Ma, Ka);
Ka=Ka_new;

% Define P1_vec                                             
    P1_vec = linspace(P_frac * P_input, P_input - Zero_prox, num_discr_P);

    % Initialize eps_vec
    eps_vec = zeros(1, length(P1_vec));
    p0_vec  = zeros(1, length(P1_vec));

    c_val=1;
    rho_temp=0;
    alpha_chosen=0.5; % Arbitrarily fixed; to be optimized, ideally
    
    % Define rho range
    rho = linspace(1e-2, 1, num_discr_rho);

    % Step 2: Define vectors
    t_vec = 1:Ka;           % Vector [1, ..., Ka]
   rho_min_opt_P_vec=zeros(1, length(P1_vec));

   %P1_vec=0.9538*P_input;

    % Loop over P1_vec and compute epsilon values
    parfor i = 1:length(P1_vec)
        P1 = P1_vec(i); % Assign current P1
        [pt,rho_min_opt]= computeTUMABoundAtP(n, Ma, Ka, R, P1, n_vec, rho, t_vec,n_S_vec, Length_S_vec, n_N_mat, Length_N_mat, Multiplicity_vec, ComplementSet_mat); % Compute epsilon
        rho_min_opt_P_vec(i)=rho_min_opt;
         p0 = Ma * gammainc(n * P_input / P1, n, 'upper');
         
         p0_vec(i)=p0;
    % Note that the multiplication factor is Ma (instead of Ka in UMA) and the
    % absence of the probability of message collision

% Compute an upper bound on the optimal lambda, to choose delta value
    % appropriately
    lambda_ub = ((P1 * c_val - 2) + sqrt((P1 * c_val - 2).^2 + 4 .* (P1 * c_val .* (1+rho_temp)))) ./ ...
                    (4 * (P1 * c_val) .* (1 + rho_temp));

    delta_ub=0.5-lambda_ub;
    delta_prime_ub=1/(1-2*delta_ub)-1;
    delta_chosen=alpha_chosen*delta_prime_ub;

    [lmax,lmin]=evaluateLambda(P1, Ka, Ma);
    delta_chosen=0.5-lmin;
    p1=exp(-n*delta_chosen^2/8);

    eps_vec(i) =pt+p0+p1;
    end

    % Find the minimum epsilon and corresponding P1
    [eps_opt, min_idx] = min(eps_vec); % Find minimum epsilon and its index
    P_opt = P1_vec(min_idx); % Corresponding optimal P1

    % Compute EbN0_dB_opt
    EbN0_dB_opt = 10 * log10(n * P_opt / (2 * k));

    p0_Popt=p0_vec(min_idx);
    










end
