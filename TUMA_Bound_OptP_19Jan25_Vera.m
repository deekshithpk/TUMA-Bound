function [P_opt, eps_opt, EbN0_dB_opt] = TUMA_Bound_OptP_19Jan25_Vera(k, n, Ma, Ka, R, num_discr_rho, P_input, P_frac, num_discr_P, Zero_prox,n_vec, n_S_vec, Length_S_vec, n_N_mat, Length_N_mat, Multiplicity_vec, ComplementSet_mat)%,n1,id_proxy)
tic

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
        [pt,rho_min_opt]= TUMA_Bound_GivenP1_19Jan25(n, Ma, Ka, R, P1, n_vec, rho, t_vec,n_S_vec, Length_S_vec, n_N_mat, Length_N_mat, Multiplicity_vec, ComplementSet_mat); % Compute epsilon
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

    [lmax,lmin]=evaluate_lambda(P1, Ka, Ma);
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
    % Return the results

  % min(rho_min_opt_P_vec)
% Return the results
% Compute EbN0_dB_opt




% EbN0_dB_input = 10 * log10(n * P_input / (2 * k));
% toc
% epsilon_final_opt=eps_opt;
% data.Ka=Ka;
% data.EbN0_dB_input=EbN0_dB_input;
% data.eps_opt=epsilon_final_opt;
% data.date_index=n1;
% data.P_input=P_input;
% data.k=k;
% data.n=n;
% data.Ma=Ma;
% sim_t=toc;
% data.sim_t_min=sim_t/60;
% 
% 
% n2=n1+id_proxy-1;




% % Define base folder path
% base_folder = '/cephyr/users/deepat/Vera/DPK';
% % Construct the folder name based on Ma
% folder_name = sprintf('TV_Ma%d_13Jan25', Ma); % Appends Ma value to folder name
% % Combine base folder and folder name to form the full path
% folder_path = fullfile(base_folder, folder_name);
% % Create the folder if it does not already exist
% if ~exist(folder_path, 'dir')
% mkdir(folder_path);
% fprintf('Created folder: %s\n', folder_path);
% else
% % fprintf('Folder already exists: %s\n', folder_path);
% end
% %folder_path = '/cephyr/users/deepat/Vera/DPK/TV_Ma40_17Dec24';
% filename = fullfile(folder_path, ['TV_' num2str(n2) '.mat']);
% save(filename,'data');
end
