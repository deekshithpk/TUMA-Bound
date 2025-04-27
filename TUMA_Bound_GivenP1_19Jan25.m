function [eps, rho_min_opt] = TUMA_Bound_GivenP1_19Jan25(n, Ma, Ka, R, P1, n_vec, rho, t_vec, n_S_vec, Length_S_vec, n_N_mat, Length_N_mat, Multiplicity_vec, ComplementSet_mat)
                                         
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
           %disp(['  Processing ell = ', num2str(ell)]); % Display current ell for debugging

            % Create S_ell_mat
           % if ell == 0
              %  S_ell_mat = zeros(1, Ma); % Special case for ell = 0
           % else
               % subsets = nchoosek(1:Ma, ell); % Generate all subsets of size ell
               % S_ell_mat = zeros(size(subsets, 1), Ma); % Initialize matrix

                %for row = 1:size(subsets, 1)
                %    S_ell_mat(row, subsets(row, :)) = 1; % Set locations to 1
               % end
            %end

            % Multiply n_vec with each row of S_ell_mat
            %n_S_ell_mat = S_ell_mat .* n_vec; % Element-wise multiplication

            % Initialize eps_nS

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


%if isempty(Ma_min_S_ell_vec)
   % disp('Ma_min_S_ell_vec is empty. No 0s in S_ell_vec.');
   % N_eye_mat = zeros(1, 1); % Special case for empty Ma_min_S_ell_vec
   % disp('N_eye_mat set to 1X1 zero matrix.');
%else
                    % Create N_eye_mat
                   % if eye == 0
                       % N_eye_mat = zeros(1, length(Ma_min_S_ell_vec)); % Special case for eye = 0
                    %else
                      %  eye_subsets = nchoosek(1:length(Ma_min_S_ell_vec), eye);
                       % N_eye_mat = zeros(size(eye_subsets, 1), length(Ma_min_S_ell_vec));

                       % for row = 1:size(eye_subsets, 1)
                           % N_eye_mat(row, eye_subsets(row, :)) = 1; % Set locations to 1
                       % end
                   % end
 %end

                   
                   % Obtain n_Ma_min_S_ell_vec
%if isempty(Ma_min_S_ell_vec)
   % disp('Ma_min_S_ell_vec is empty. No 0s in S_ell_vec.');
   % n_Ma_min_S_ell_vec = 0; % Set to an empty vector for this case
  %  disp('n_Ma_min_S_ell_vec set to 0.');
%else
%    n_Ma_min_S_ell_vec = n_vec(Ma_min_S_ell_vec); % Extract relevant elements from n_vec
%end
                    

                    % Multiply to create n_N_eye_mat
                    %n_N_eye_mat = N_eye_mat .* n_Ma_min_S_ell_vec;


                    % Step 13: Loop through rows of n_N_eye_mat
                   % n_N_eye_1norm_vec = zeros(size(n_N_eye_mat, 1), 1); % Initialize vector
                   % n_Nhat_Ma_min_ell_min_eye_1norm_vec = Ka-n_N_eye_1norm_vec-n_S_ell_1norm; 

%                     for n_row = 1:size(n_N_eye_mat, 1)
%                         % Fetch current row of n_N_eye_mat
%                         n_N_eye_vec = n_N_eye_mat(n_row, :);
%                         n_N_eye_1norm = sum(n_N_eye_vec); % Sum of the row
% 
%                         % Compute n_Nhat_Ma_min_ell_min_eye_1norm directly
%                         n_Nhat_Ma_min_ell_min_eye_1norm = Ka - n_N_eye_1norm - n_S_ell_1norm;
% 
%                         % Store the results in the respective vectors
%                         n_N_eye_1norm_vec(n_row) = n_N_eye_1norm;
%                         n_Nhat_Ma_min_ell_min_eye_1norm_vec(n_row) = n_Nhat_Ma_min_ell_min_eye_1norm;
%                     end

                   % Define jay_vec with condition for ell, eye
if ell > 0 && eye == Ma - ell
    jay_vec = t;  % Set jay_vec to t if ell > 0 and eye equals Ma - ell
elseif ell == 0
    jay_vec = 0;  % If ell equals 0, set jay_vec to 0
else
    jay_vec = ell:(t - Ma + ell + eye);  % Otherwise, define jay_vec with the range
end



                    % Call the new function TUMA_Bound_9Jan25
                    [p_t_ell_nS_eye_sum_nN_j_opt_rho,min_optimal_rho_value] = TUMA_Bound_19Jan25(t, ell, n_S_ell_1norm, eye, n_N_eye_1norm_vec, jay_vec, rho, P1, n, Ma, Ka, R);
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
