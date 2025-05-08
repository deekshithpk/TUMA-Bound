function [eps, rho_min_opt] = computeTUMABoundAtP(n, Ma, Ka, R, P1, n_vec, rho, t_vec, n_S_vec, Length_S_vec, n_N_mat, Length_N_mat, Multiplicity_vec, ComplementSet_mat)
 
% -------------------------------------------------------------------------
% Function: computeTUMABoundAtP
%
% Purpose:
%   Computes the TUMA bound for a fixed transmit power level P1,
%   given system parameters.
%   This function is used as a helper within computeOptimalBoundOverP to
%   evaluate the bound across a range of P values.
%
% Inputs:
%   - n          : Block length
%   - Ma, Ka     : Number of active messages and users
%   - R          : Transmission rate
%   - P1         : Fixed transmit power level (manifests as codebook variance in the bound)
%   - n_vec, rho : Multiplicity vector, optimization parameter,
%                  respectively
%   - t_vec, n_S_vec, Length_S_vec, n_N_mat, Length_N_mat
%   - Multiplicity_vec, ComplementSet_mat : Subset info
%
% Outputs:
%   - eps        : Computed bound value (e.g., TV distance or related metric)
%   - rho_min_opt: Optimal rho value used in the bound computation
%
% Notes:
%   Intended for internal use within power optimization
% -------------------------------------------------------------------------




    rho_min_vec=[];

    
    % Initialize eps
    eps = 0;

    % Step 3: Create loops for t_vec and ell_vec
    for t = t_vec
       %disp(['Processing t = ', num2str(t)]); % Display current t for debugging

        % Initialize p_t
        p_t = 0;

        if t == Ka
        ell_vec = Ma;           % Set ell_vec to Ma (a scalar) if t equals Ka
        else
        ell_vec = 0:Ma;     % Define the vector [0, ..., Ma] otherwise
        end


        for ell = ell_vec

            if ell > t
                    continue;
            end            
           

            % Find the indices of entries in Length_S_vec that match ell
matching_indices = find(Length_S_vec == ell);

% Use these indices to find the corresponding entries in n_S_vec
matching_n_S_values = n_S_vec(matching_indices);


            p_t_ell = 0;

            % Step 10-11: Loop through rows of n_S_ell_mat
            for row_index = matching_indices
                n_S_ell_1norm = n_S_vec(row_index);

% Initialize the vectors for the current row
temp_n_N_vec = n_N_mat(row_index, :);
temp_n_len_vec = Length_N_mat(row_index, :);

% Remove the first element, if the length is greater than 1
if length(temp_n_N_vec) > 1
    temp_n_N_vec = temp_n_N_vec(2:end);
end
if length(temp_n_len_vec) > 1
    temp_n_len_vec = temp_n_len_vec(2:end);
end

% Keep only non-zero elements from the truncated vectors
curr_n_sum_vec = temp_n_N_vec(temp_n_N_vec > 0);
curr_n_len_vec = temp_n_len_vec(temp_n_len_vec > 0);

% Append a zero at the beginning of the vectors
curr_n_sum_vec = [0, curr_n_sum_vec];
curr_n_len_vec = [0, curr_n_len_vec];



             



                if n_S_ell_1norm > t || ell > n_S_ell_1norm
                continue;
                end


                %Ma_min_S_ell_vec = ComplementSet_mat(row_index, :); % Locations of 0s
                
                if ell > 0 && t == n_S_ell_1norm
                    eye_start = 0;         % Set eye_start to 0
                    eye_end = 0;           % Set eye_end to 0
                elseif ell > 0
                     eye_start = max(0,Ma-t);         % Start from 0 when ell > 0
                     eye_end = Ma - ell;    % End at Ma - ell
                else
                     eye_start = max(1,Ma-t);         % Start from 1 when ell == 0
                     eye_end = Ma - 1;      % End at Ma - 1
                end


                % Initialize eps_eye
                p_t_ell_nS = 0;

                for eye = eye_start:eye_end

                  matching_indices1 = find(curr_n_len_vec == eye);

n_N_eye_1norm_vec = curr_n_sum_vec(matching_indices1);




                   % Define jay_vec with condition for ell, eye
if ell > 0 && eye == Ma - ell
    jay_vec = t;  % Set jay_vec to t if ell > 0 and eye equals Ma - ell
elseif ell == 0
    jay_vec = 0;  % If ell equals 0, set jay_vec to 0
else
    jay_vec = ell:(t - Ma + ell + eye);  % Otherwise, define jay_vec with the range
end



                    % Call the  function computeTUMABoundCore
                    [p_t_ell_nS_eye_sum_nN_j_opt_rho,min_optimal_rho_value] = computeTUMABoundCore(t, ell, n_S_ell_1norm, eye, n_N_eye_1norm_vec, jay_vec, rho, P1, n, Ma, Ka, R);
rho_min_vec=[rho_min_vec' min_optimal_rho_value ]';
                    % Update eps_eye
                    p_t_ell_nS = p_t_ell_nS + p_t_ell_nS_eye_sum_nN_j_opt_rho;
                end

                % Update eps_nS
                p_t_ell = p_t_ell + Multiplicity_vec(row_index)*p_t_ell_nS;
            end

            % Update eps_ell
            p_t = p_t + p_t_ell;
        end

        % Update eps
        eps = eps + t / Ka * p_t;
    end
    rho_min_opt=min(rho_min_vec);
end
