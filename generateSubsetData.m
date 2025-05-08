
% -------------------------------------------------------------------------
% Script Name: generateSubsetData.m
%
% Purpose:
%   Generates and saves subsets of a multiplicity vector. This script serves as 
%   a helper to main_compute_tv_bound_zipf_type.m by precomputing and storing 
%   subset structures in a .mat file used during bound evaluation. 
%   Pre-generating these subsets improves runtime efficiency during simulations.


% Usage:
%   - Set Ma (number of active messages) and Ka (total active users).
%   - Choose between manual or automatic generation of n_vec.
%   - Run the script to generate and save the structured data.
%
% Output:
%   - A .mat file containing a struct array `Subset_Data_Ma_<Ma>_Ka_<Ka>_<Date>.mat`
%     with subset sums, complement sets, multiplicities, and related metadata.
%
% Notes:
%   - Ensure generateMultiplicityVector.m is available in the path.
%   - This script is typically run once per (Ma, Ka) pair before bound evaluation.
%   - Some steps (e.g., subset matrix generation) may be computationally intensive.
%
% Last modified: 2025-05-06
% -------------------------------------------------------------------------






tic; % Start timer

clc
clear variables

% Step 1: Initialize parameters and generate n_vec
use_manual_initialization = false; % Set to false for automatic generation

if use_manual_initialization
    % Manual Initialization
    Ma = 3; % Number of elements
    Ka = 5; % Total sum
    n_vec = [3, 1, 1]; % Example n_vec
else
    % Automatic generation
    Ma = 2; % Example value
    Ka = 5; % Example value
    [n_vec, Ka_new] = generateMultiplicityVector(Ma, Ka);
end

% Ensure n_vec is sorted in decreasing order
n_vec = sort(n_vec, 'descend');

% Step 4: Split n_vec into two parts: n_vec_a and n_vec_ones
n_vec_a = n_vec(n_vec > 1); % Elements greater than 1
n_vec_ones = ones(1, sum(n_vec == 1)); % All ones

% Generate pairs for n_vec_a
pairs_n_vec_a = generatePairs(n_vec_a);

% Generate pairs for n_vec_ones
pairs_n_vec_ones = generatePairsForOnes(n_vec_ones);

% Step 5: Combine results from n_vec_a and n_vec_ones
pairs_combined = combinePairs(pairs_n_vec_a, pairs_n_vec_ones, n_vec);



% Step 7: Update n_N_ell1_norm_vec for each row based on ComplementSet
for row_num = 1:length(pairs_combined)
    % Update n_N_ell1_norm_vec with subset sums of the ComplementSet
    [subset_sums, subset_lengths, multiplicities, length_multiplicities]=computeSubsetSumsOptimized(pairs_combined(row_num).ComplementSet);
    pairs_combined(row_num).n_N_ell1_norm_vec = [subset_sums multiplicities];
    pairs_combined(row_num).SubsetLengths = [subset_lengths length_multiplicities] ;
end

% Step 6: Rearrange rows of pairs variable
pairs_combined_sorted = sortPairs(pairs_combined);

% Measure and display elapsed time in appropriate units
elapsed_time = toc; % Stop timer

if elapsed_time >= 60
    elapsed_time_minutes = elapsed_time / 60; % Convert to minutes
    fprintf('Elapsed Time: %.2f minutes\n', elapsed_time_minutes); % Display in minutes
else
    elapsed_time_seconds = elapsed_time; % Keep in seconds
    fprintf('Elapsed Time: %.2f seconds\n', elapsed_time_seconds); % Display in seconds
end

% Generate today's date
todays_date = datestr(now, 'ddmmmyyyy'); % Format: 17Jan2025

% Create variable name with Ma, Ka, and date
output_variable_name = sprintf('Subset_Data_Ma_%d_Ka_%d_%s', Ma, Ka, todays_date);

% Save the variable with the generated name
eval(sprintf('%s = pairs_combined_sorted;', output_variable_name)); % Rename variable
save(sprintf('%s.mat', output_variable_name), output_variable_name); % Save to .mat file


%% ---------------- Helper Functions ---------------- %%
% Function to generate pairs for any general n_vec
function pairs = generatePairs(n_vec)
    subset_matrix = generateSubsetMatrix(length(n_vec));
    subset_complement_matrix = 1 - subset_matrix;

    pairs = struct('n_S', {}, 'Length_S', {}, 'Subset_S', {}, 'ComplementSet', {}, ...
                   'n_N_ell1_norm_vec', {}, 'SubsetLengths', {}, 'Multiplicity', {});

    for r = 1:size(subset_matrix, 1)
        subset_S = subset_matrix(r, :);
        n_S = sum(subset_S .* n_vec);
        complement_indices = find(~subset_S);
        complement_set = n_vec(complement_indices);

        subset_matrix_for_S_complement = generateSubsetMatrix(length(complement_indices));
        subset_values_complement = subset_matrix_for_S_complement .* complement_set;

        n_N_ell1_norm_vec = sum(subset_values_complement, 2);
        subset_lengths_complement = sum(subset_matrix_for_S_complement, 2);

        pairs(r).n_S = n_S;
        pairs(r).Length_S = sum(subset_S);
        pairs(r).Subset_S = n_vec(logical(subset_S));
        pairs(r).ComplementSet = complement_set;
        pairs(r).n_N_ell1_norm_vec = n_N_ell1_norm_vec';
        pairs(r).SubsetLengths = subset_lengths_complement';
        pairs(r).Multiplicity = 1; % Uniform multiplicity initially
    end
end

% Function to update multiplicities in pairs_full
function pairs_updated = updateMultiplicities(pairs)
    % Dictionary to track unique pairs and their multiplicities
    pair_dict = containers.Map('KeyType', 'char', 'ValueType', 'double');

    % Iterate through pairs to calculate multiplicities
    for r = 1:length(pairs)
        current_n_S = pairs(r).n_S;
        current_norm_vec = pairs(r).n_N_ell1_norm_vec;

        % Create a unique key for the pair
        pair_key = sprintf('n_S=%d, n_N_ell1_norm_vec=[%s]', current_n_S, num2str(current_norm_vec));

        % Update multiplicity in the dictionary
        if isKey(pair_dict, pair_key)
            pair_dict(pair_key) = pair_dict(pair_key) + 1;
        else
            pair_dict(pair_key) = 1;
        end
    end

    % Iterate again to populate updated pairs with new multiplicities
    pairs_updated = pairs;
    keys = pair_dict.keys;

    for r = 1:length(pairs)
        current_n_S = pairs(r).n_S;
        current_norm_vec = pairs(r).n_N_ell1_norm_vec;

        % Create a unique key for the pair
        pair_key = sprintf('n_S=%d, n_N_ell1_norm_vec=[%s]', current_n_S, num2str(current_norm_vec));

        % Assign updated multiplicity
        pairs_updated(r).Multiplicity = pair_dict(pair_key);
    end

    % Filter pairs to keep only unique rows
    pairs_updated = uniqueStructArray(pairs_updated, {'n_S', 'n_N_ell1_norm_vec'});
end

% Function to filter unique rows in a struct array
function unique_array = uniqueStructArray(struct_array, fields)
    % Generate unique keys for the fields
    keys = arrayfun(@(x) sprintf('%s', struct2str(x, fields)), struct_array, 'UniformOutput', false);
    [~, unique_indices] = unique(keys);
    unique_array = struct_array(unique_indices);
end

% Function to combine pairs from n_vec_a and n_vec_ones
function pairs_combined = combinePairs(pairs_a, pairs_ones,n_vec)
    pairs_combined = [];
    counter = 1;
    for r = 1:length(pairs_ones)
        num_ones = r - 1;
        for row = 1:length(pairs_a)
            combined_subset_s = [pairs_a(row).Subset_S, ones(1, num_ones)];
            combined_n_s = pairs_a(row).n_S + num_ones;

combined_complement_set = multisetDiff(n_vec, combined_subset_s);


            % Step 1: Generate all subsets as a zero-one matrix
num_elements = length(combined_complement_set); % Number of elements in the set
subset_matrix = [];% generateSubsetMatrix(num_elements); % Generate subset matrix

% Step 2: Multiply zero-one matrix with `combined_complement_set`
%subset_values_matrix = subset_matrix .* combined_complement_set;

% Step 3: Compute the sum of elements of each subset (row-wise sum)
combined_n_N_ell1_norm_vec = 0;% sum(subset_values_matrix, 2);
SubsetLengths1 = sum(subset_matrix, 2);

            combined_multiplicity = pairs_a(row).Multiplicity * pairs_ones(r).Multiplicity;

            pairs_combined(counter).n_S = combined_n_s;
            pairs_combined(counter).Length_S = pairs_a(row).Length_S + num_ones;
            pairs_combined(counter).Subset_S = combined_subset_s;
            pairs_combined(counter).ComplementSet = combined_complement_set;
            pairs_combined(counter).n_N_ell1_norm_vec = combined_n_N_ell1_norm_vec';
            pairs_combined(counter).SubsetLengths = SubsetLengths1';
            pairs_combined(counter).Multiplicity = combined_multiplicity;
            counter = counter + 1;
        end
    end
end

% Function to generate the subset matrix for {1, ..., Ma}
function subset_matrix = generateSubsetMatrix(u)
    num_subsets = 2^u;
    subset_matrix = false(num_subsets, u);
    for i = 0:num_subsets-1
        binary_representation = dec2bin(i, u) - '0';
        subset_matrix(i + 1, :) = binary_representation;
    end
end

% Utility function to create string from struct fields
function str = struct2str(struct, fields)
    values = arrayfun(@(f) struct.(f{1}), fields, 'UniformOutput', false);
    str = jsonencode(values);
end
% Function to generate pairs for n_vec consisting only of ones
function pairs = generatePairsForOnes(n_vec_ones)
    len = length(n_vec_ones);
    pairs = struct('n_S', {}, 'Length_S', {}, 'Subset_S', {}, 'ComplementSet', {}, ...
                   'n_N_ell1_norm_vec', {}, 'SubsetLengths', {}, 'Multiplicity', {});

    % Generate rows for each possible subset length
    for row = 0:len
        pairs(row + 1).n_S = row; % Subset sum (row index)
        pairs(row + 1).Length_S = row; % Length of subset S
        pairs(row + 1).Subset_S = ones(1, row); % Subset S values
        pairs(row + 1).ComplementSet = ones(1, len - row); % Complement set values
        pairs(row + 1).n_N_ell1_norm_vec = 0:(len - row); % Vector of sums
        pairs(row + 1).SubsetLengths = 0:(len - row); % Subset lengths
        pairs(row + 1).Multiplicity = nchoosek(len, row); % Binomial coefficient
    end
end

function complement_set = multisetDiff(n_vec, subset_s)
    % Custom implementation of set difference for multisets
    complement_set = n_vec; % Start with the full set
    for elem = subset_s
        % Remove one occurrence of `elem` from `complement_set`
        idx = find(complement_set == elem, 1, 'first');
        if ~isempty(idx)
            complement_set(idx) = []; % Remove the first occurrence
        end
    end
end

%% ---------------- Helper Function to Sort Pairs ---------------- %%
function sorted_pairs = sortPairs(pairs)
    % Extract n_S and Length_S for sorting
    n_S_values = arrayfun(@(x) x.n_S, pairs);
    length_S_values = arrayfun(@(x) x.Length_S, pairs);

    % Create sorting indices
    [~, sort_indices] = sortrows([n_S_values(:), length_S_values(:)]);

    % Reorder pairs based on sorting indices
    sorted_pairs_temp = pairs(sort_indices);

    % Combine rows with the same (n_S, n_N_ell1_norm_vec)
    pair_dict = containers.Map('KeyType', 'char', 'ValueType', 'double');
    sorted_pairs = struct('n_S', {}, 'Length_S', {}, 'Subset_S', {}, ...
                          'ComplementSet', {}, 'n_N_ell1_norm_vec', {}, ...
                          'SubsetLengths', {}, 'Multiplicity', {});

    for i = 1:length(sorted_pairs_temp)
        % Create a unique key based on (n_S, n_N_ell1_norm_vec)
        current_n_S = sorted_pairs_temp(i).n_S;
        current_norm_vec = sorted_pairs_temp(i).n_N_ell1_norm_vec;

        % Create key for the row
        key = sprintf('n_S=%d, n_N_ell1_norm_vec=[%s]', current_n_S, num2str(current_norm_vec));

        if isKey(pair_dict, key)
            % Update multiplicity if key already exists
            pair_dict(key) = pair_dict(key) + sorted_pairs_temp(i).Multiplicity;
        else
            % Add a new entry
            pair_dict(key) = sorted_pairs_temp(i).Multiplicity;

            % Copy other fields for the new row
            sorted_pairs(end + 1).n_S = current_n_S;
            sorted_pairs(end).Length_S = sorted_pairs_temp(i).Length_S;
            sorted_pairs(end).Subset_S = sorted_pairs_temp(i).Subset_S;
            sorted_pairs(end).ComplementSet = sorted_pairs_temp(i).ComplementSet;
            sorted_pairs(end).n_N_ell1_norm_vec = current_norm_vec;
            sorted_pairs(end).SubsetLengths = sorted_pairs_temp(i).SubsetLengths;
            sorted_pairs(end).Multiplicity = 0; % Initialize to update later
        end
    end

    % Update the multiplicities for unique rows
    for i = 1:length(sorted_pairs)
        % Create the key for each unique row
        key = sprintf('n_S=%d, n_N_ell1_norm_vec=[%s]', ...
                      sorted_pairs(i).n_S, ...
                      num2str(sorted_pairs(i).n_N_ell1_norm_vec));
        sorted_pairs(i).Multiplicity = pair_dict(key);
    end
end


