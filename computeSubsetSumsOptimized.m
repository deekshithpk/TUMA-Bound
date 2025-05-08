function [subset_sums, subset_lengths, multiplicities, length_multiplicities] = computeSubsetSumsOptimized(n_vec)
    % -------------------------------------------------------------------------
% Function: computeSubsetSumsOptimized
%
% Purpose:
%   Efficiently computes the sums of all possible subsets of the input 
%   vector n_vec, along with corresponding subset lengths and their 
%   multiplicities. Optimized for cases where n_vec contains many 1's.
%
% Inputs:
%   - n_vec: A vector of positive integers (e.g., [2, 3, 1, 1, 1])
%
% Outputs:
%   - subset_sums        : Unique values of subset sums (sorted)
%   - subset_lengths     : Corresponding lengths of subsets
%   - multiplicities     : Number of times each subset sum appears
%   - length_multiplicities : Number of times each subset length appears
% -------------------------------------------------------------------------


    % Separate the vector into non-ones and count of ones
    non_ones = n_vec(n_vec > 1); % Non-one elements
    ones_count = sum(n_vec == 1); % Count of ones

    % Step 1: Compute subset sums for non-ones
    num_subsets = 2^length(non_ones); % Total number of subsets for non-ones
    subset_sums_non_ones = zeros(num_subsets, 1);
    subset_lengths_non_ones = zeros(num_subsets, 1);

    % Generate all subset sums and lengths for non-ones
    for i = 1:num_subsets
        % Logical indices corresponding to the binary representation of i-1
        indices = logical(bitget(i - 1, 1:length(non_ones)));

        % Sum the elements corresponding to the indices
        subset_sums_non_ones(i) = sum(non_ones(indices));
        
        % Length of the subset (number of selected elements)
        subset_lengths_non_ones(i) = sum(indices);
    end

    % Step 2: Compute sums for subsets of ones
    % The sum of any subset of `1`s is just its cardinality (0 to ones_count)
    subset_sums_ones = 0:ones_count;
    subset_lengths_ones = 0:ones_count; % Length is the cardinality for ones

    % Step 3: Combine the two subset sums and lengths
    % Use broadcasting-compatible addition and concatenation
    [grid_non_ones, grid_ones] = ndgrid(subset_sums_non_ones, subset_sums_ones);
    combined_sums = grid_non_ones + grid_ones;

    [grid_lengths_non_ones, grid_lengths_ones] = ndgrid(subset_lengths_non_ones, subset_lengths_ones);
    combined_lengths = grid_lengths_non_ones + grid_lengths_ones;



    % Flatten the results
    combined_sums = combined_sums(:);
    combined_lengths = combined_lengths(:);

    % Step 4: Compute unique subset sums, lengths, and multiplicities
    [subset_sums, ~, idx] = unique(combined_sums, 'sorted');
     [subset_lengths, ~, idx1] = unique(combined_lengths, 'sorted');
    multiplicities = accumarray(idx, 1); % Count occurrences of each unique sum
    length_multiplicities = accumarray(idx1, 1); % Aggregate lengths
end

