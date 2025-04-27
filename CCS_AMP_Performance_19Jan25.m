% MATLAB code to process data and generate the desired results

% Data (EbN0, TV Distance) for each Ka and Ma
% Organized as a structure for clarity
data = struct( ...
    'Ka100_Ma100', [3.4, 0.14610000000000004; 3.8, 0.06749999999999996; 4.0, 0.049599999999999964], ...
    'Ka92_Ma80', [3.0, 0.13883110884036498; 3.4, 0.08468733787986174; 3.6, 0.06712718272821458; 3.8, 0.04533025779253421; 4.0, 0.03708091698667428], ...
    'Ka98_Ma60', [4.0, 0.08561035692106184; 4.4, 0.06482594027150647; 4.6, 0.05840677497307235; 4.8, 0.05357418412979602], ...
    'Ka102_Ma50', [1.0, 0.46912607167534737; 3.0, 0.09970755261240036; 3.4, 0.08160975234594561; 4.0, 0.055287397209886284; 4.2, 0.052226799559560895; 4.4, 0.042396808790774836], ...
    'Ka102_Ma30', [3.0, 0.07677950064181673; 3.6, 0.0600092985550107; 3.8, 0.04121754016013672; 4.0, 0.037845974565248755], ...
    'Ka100_Ma10', [-10.0, 0.10815755702397001; -9.4, 0.0873412782015134; -9.0, 0.08097871255663236; -8.0, 0.060490582970946206; -7.4, 0.059048202657295906; -7.0, 0.043263611112618686], ...
    'Ka100_Ma2', [-25.0, 0.07691721137987899; -23.0, 0.05495658388945679; -22.0, 0.045738529014662556; -21.0, 0.04054294453691564; -20.0, 0.03464223139936609] ...
);

% Target TV distance
target_TV_distance = 0.05;

% Initialize vectors to store results
Ma_values = [100, 80, 60, 50, 30, 10, 2];
Ka_values = [100, 92, 98, 102, 102, 100, 100];
EbN0_values = [];

% Process data to find EbN0 corresponding to TV distance = 0.05
field_names = {'Ka100_Ma100', 'Ka92_Ma80', 'Ka98_Ma60', 'Ka102_Ma50', 'Ka102_Ma30', 'Ka100_Ma10', 'Ka100_Ma2'};
for i = 1:length(field_names)
    table = data.(field_names{i});
    EbN0 = interp1(table(:, 2), table(:, 1), target_TV_distance, 'linear', 'extrap');
    EbN0_values = [EbN0_values; EbN0];
end

% TV distances between UMA profile and TUMA profile for each Ma value
UMA_TUMA_TV_distances = [0, 0.1304, 0.3878, 0.5098, 0.7, 0.9, 0.98];

% Generate the plot
figure;
plot(UMA_TUMA_TV_distances, EbN0_values, '-o', 'LineWidth', 1.5);
grid on;
hold on;

% Add annotations
text(UMA_TUMA_TV_distances, EbN0_values, string(Ma_values), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Labels and title
xlabel('TV distance (UMA profile, TUMA profile)');
ylabel('TV distance (transmitted profile, decoded profile)');
title('CCS-AMP Performance vs Single \rho-trick Bound');

% Add a box with general information
annotation('textbox', [0.3, 0., 0.2, 0.15], 'String', {
    'General Info:',
    'k = 128',
    'n = 38400',
    'Target TV distance = 0.05'
}, 'EdgeColor', 'black', 'BackgroundColor', 'white');

% Save the results
Ka_prime = 100;
Ka_vec = Ka_values;
save('results.mat', 'EbN0_values', 'UMA_TUMA_TV_distances', 'Ka_prime', 'Ka_vec');

% Additional data for "Single rho trick bound"
EbN0_db_single_rho = [0.0452, 1.51, 1.44, 1.4, 1.3, -7.11, -27.32];

% Hold on to the previous figure and plot the new data
hold on;
plot(UMA_TUMA_TV_distances, EbN0_db_single_rho, '--s', 'LineWidth', 1.5, 'MarkerSize', 8);

% Add legend to differentiate between the plots
legend('CCS-AMP Performance', 'Single \rho-trick bound', 'Location', 'best');

% Ensure grid is still visible
grid on;

% Save the updated plot
saveas(gcf, 'CCS_AMP_Performance_Updated.png');

