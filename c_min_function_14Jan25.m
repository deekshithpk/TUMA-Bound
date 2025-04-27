function c_min_values = c_min_function_14Jan25(t, ell, eye, jay_grid, Ma, n_S_ell_1norm)
    
    % First term: (ell ~= 0) .* (n_S_ell_1norm^2) ./ ell
    if ell ~= 0
        first_term = (n_S_ell_1norm^2) / ell;
    else
        first_term = 0;
    end
    
    % Second term: (eye ~= 0) .* ceil((t - (n_S_ell_1norm)^2) ./ eye)
    if eye ~= 0
        second_term = ceil((t - (n_S_ell_1norm))^2 / eye);
    else
        second_term = 0;
    end
    
    % Third term: ((Ma - ell - eye) ~= 0) .* ceil((t - jay_grid).^2 ./ (Ma - ell - eye))
    if (Ma - ell - eye) ~= 0
        third_term = ceil((t - jay_grid).^2 ./ (Ma - ell - eye));
    else
        third_term = zeros(size(jay_grid)); % Set to zeros with the same size as jay_grid
    end
    
    % Fourth term: (ell ~= 0) .* ceil(jay_grid.^2 ./ ell)
    if ell ~= 0
        fourth_term = ceil(jay_grid.^2 ./ ell);
    else
        fourth_term = zeros(size(jay_grid)); % Set to zeros with the same size as jay_grid
    end
    
    % Combine all terms
    c_min_values = first_term + second_term + third_term + fourth_term;
end
