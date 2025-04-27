function generateSquareGrid()
    % Generate a 4x5 grid of square boxes with random numbers (1-20) in each box.

    % Randomize a 4x5 grid with unique numbers from 1 to 20
    gridNumbers = reshape(randperm(20), 5, 4)';

    % Create figure
    figure('Color', 'w', 'Position', [100, 100, 500, 400]);
    hold on;
    axis off;

    % Define box dimensions
    rows = 4; % Number of rows
    cols = 5; % Number of columns
    boxSize = 1; % Size of each box

    % Draw grid and display numbers
    for row = 1:rows
        for col = 1:cols
            % Compute box position
            xPos = col - 1; % X-coordinate
            yPos = -(row - 1); % Y-coordinate (negative for downward placement)

            % Draw the box with a light color
            rectangle('Position', [xPos, yPos, boxSize, boxSize], ...
                'FaceColor', [0.9, 0.95, 1], 'EdgeColor', 'k');

            % Display the number at the center of the box
            text(xPos + boxSize / 2, yPos - boxSize / 2, ...
                num2str(gridNumbers(row, col)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
        end
    end

    % Set the axis limits to match the grid
    xlim([0, cols]);
    ylim([-rows, 0]);

    % Title
    title('4x5 Grid with Random Numbers', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
end
