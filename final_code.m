function main()
    % Initialize the environment
    [gridSize, startPos, mapGrid,radiationMap,radiationGrid, occupancyGrid,logOddsOGM,ax] = initialize_environment();

    % Call the exploration function
    explore_grid_multi(gridSize, startPos, mapGrid,radiationMap,radiationGrid, occupancyGrid,logOddsOGM,ax);
end



%/////////////////////////////////////////////////////
%moving




function [sensorData, visibleCells, costMap, parentMap, radiationMap, occupancyGrid, logOddsOGM] = ...
    sensor_scan(robotPos, mapGrid, gridSize, costMap, parentMap, radiationMap, radiationGrid, occupancyGrid, logOddsOGM, directions)

    % Sensor radius
    radius = 3;
    [rows, cols] = meshgrid(-radius:radius, -radius:radius);
    validCells = rows.^2 + cols.^2 <= radius^2;
    rows = rows(validCells);
    cols = cols(validCells);

    % Initialize outputs
    visibleCells = [];
    sensorData = [];
    % Iterate over boundary cells only
    for i = 1:length(rows)
        % Calculate target position
        targetRow = robotPos(1) + rows(i);
        targetCol = robotPos(2) + cols(i);
        

        % Check if within bounds
        if targetRow > 0 && targetRow <= gridSize(1) && ...
           targetCol > 0 && targetCol <= gridSize(2)
            

            % Use line-of-sight to update OGM and detect visibility
            [visible, occupancyGrid, logOddsOGM] = has_line_of_sight(robotPos, [targetRow, targetCol], mapGrid, occupancyGrid, logOddsOGM);
            % If visible, process the cell
            if visible
                % Add to visible cells
                visibleCells = [visibleCells; targetRow, targetCol]; %#ok<AGROW>
                % Add to sensor data (environmental info)

                % Update radiation OGM (direct mapping)
                radiationMap(targetRow, targetCol) = radiationGrid(targetRow, targetCol);

                % Update Dijkstra's cost map for shortest path
                currentCost = costMap(robotPos(1), robotPos(2)) + abs(rows(i)) + abs(cols(i));
                if occupancyGrid(targetRow,targetCol) < 0.9 
                    if currentCost < costMap(targetRow, targetCol)
                        costMap(targetRow, targetCol) = currentCost;
                        parentMap(targetRow, targetCol) = sub2ind(gridSize, robotPos(1), robotPos(2));
                    end
                else
                    costMap(targetRow, targetCol) = inf;
                    parentMap(targetRow, targetCol) = sub2ind(gridSize, robotPos(1), robotPos(2));
                end
            
                
            end
        end
    end
end



function [visible, occupancyGrid, logOddsOGM] = has_line_of_sight(robotPos, targetPos, mapGrid, occupancyGrid, logOddsOGM)
    % Bresenham's line algorithm to check for line-of-sight
    x1 = robotPos(2); % Column index of robot position
    y1 = robotPos(1); % Row index of robot position
    x2 = targetPos(2); % Column index of target position
    y2 = targetPos(1); % Row index of target position
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    sx = sign(x2 - x1);
    sy = sign(y2 - y1);

    err = dx - dy;
    sensorFailureRate =0.1;

    % Define log-odds updates
    logOddsOccupied = log(0.9 / 0.1); % Log-odds for high confidence in occupancy
    logOddsFree = log(0.1 / 0.9);     % Log-odds for high confidence in free space

    % Track whether any obstacles were encountered

    while true
        % Check if we've reached the target
        if x1 == x2 && y1 == y2
            % If the target is a wall but there are no obstacles in the line of sight
            perceivedState = mapGrid(y1, x1);
            if rand() < sensorFailureRate
                perceivedState = 1 - perceivedState; % Flip state (free becomes occupied, occupied becomes free)
            end

            % Perform Bayesian update for the target cell
            if perceivedState == 1 % Wall detected
                visible = true; % Still mark it as visible
                logOddsOGM(y1, x1) = logOddsOGM(y1, x1) + logOddsOccupied; % Mark as an obstacle
                occupancyGrid(y1, x1) = 1 / (1 + exp(-logOddsOGM(y1, x1))); % Update probability
            else
                visible = true; % Mark as visible for free space
                logOddsOGM(y1, x1) = logOddsOGM(y1, x1) + logOddsFree; % Decrease confidence in occupancy
                occupancyGrid(y1, x1) = 1 / (1 + exp(-logOddsOGM(y1, x1)));

            end
            return;
        end

        % Check if the current point is an obstacle
        if mapGrid(y1, x1) == 1

            % Update OGM for the detected obstacle using log-odds
            logOddsOGM(y1, x1) = logOddsOGM(y1, x1) + logOddsOccupied; % Increase confidence for obstacle
            occupancyGrid(y1, x1) = 1 / (1 + exp(-logOddsOGM(y1, x1))); % Convert back to probability

            visible = false; % Line of sight is blocked
            return;
        end

        % Update cells along the way as free space
         % Convert back to probability
        % Step to the next point along the line
        e2 = 2 * err;
        if e2 > -dy
            err = err - dy;
            x1 = x1 + sx;
        end
        if e2 < dx
            err = err + dx;
            y1 = y1 + sy;
        end
    end
end








function path = construct_bfs_path(parent, target, start)
    % Reconstruct the path from BFS parent map
    path = target;


    while ~isequal(path(1, :), start)
        currentStr = mat2str(path(1, :));
        parentStr = mat2str(parent(currentStr));
        parentCell = sscanf(parentStr, '[%d %d]'); % Convert string back to array
        path = [reshape(parentCell, 1, 2); path];
    end

    path = path(2, :); % Return the next step
end










function path = reconstruct_path_dijkstra(costMap, start, goal, parentMap)
    % Reconstruct path from goal to start using parentMap
    current = goal;
    path = goal; % Initialize path with the goal
    while ~isequal(current, start)
        parentIdx = parentMap(current(1), current(2));
        [row, col] = ind2sub(size(parentMap), parentIdx);
        current = [row, col];
        path = [current; path]; % Prepend to path
    end
end


function pathsToStart = compute_return_paths(gridSize, startPosList, robotPositions, occupancyGrid)
    % Compute Dijkstra paths for all robots to return to their start positions
    numRobots = size(startPosList, 1);
    pathsToStart = cell(numRobots, 1);

    for i = 1:numRobots
        start = robotPositions(i, :);
        goal = startPosList(i, :);

        % Perform Dijkstra's algorithm
        [costMap, parentMap] = dijkstra(gridSize, start, goal, occupancyGrid);

        % Reconstruct the path from the parentMap
        pathsToStart{i} = reconstruct_path_dijkstra(costMap, start, goal, parentMap);
    end
end

function [costMap, parentMap] = dijkstra(gridSize, start, goal, occupancyGrid)
    % Initialize cost map and parent map
    costMap = inf(gridSize);
    parentMap = zeros(gridSize);
    visited = false(gridSize); % Track visited nodes

    % Start position
    costMap(start(1), start(2)) = 0;

    % Open set (list of nodes to explore)
    openSet = [start, 0]; % Each row: [row, col, cost]

    % Movement directions (4-connectivity)
    directions = [0, 1; 1, 0; 0, -1; -1, 0];

    while ~isempty(openSet)
        % Find the node with the smallest cost
        [~, idx] = min(openSet(:, 3));
        current = openSet(idx, 1:2);
        currentCost = openSet(idx, 3);
        openSet(idx, :) = []; % Remove current node from open set

        % Mark as visited
        visited(current(1), current(2)) = true;

        % Stop if we reach the goal
        if isequal(current, goal)
            break;
        end

        % Explore neighbors
        for d = 1:size(directions, 1)
            neighbor = current + directions(d, :);
            if neighbor(1) > 0 && neighbor(1) <= gridSize(1) && ...
               neighbor(2) > 0 && neighbor(2) <= gridSize(2) && ...
               ~visited(neighbor(1), neighbor(2)) && ...
               occupancyGrid(neighbor(1), neighbor(2)) < 0.9 % Not a wall

                newCost = currentCost + 1; % Uniform cost
                if newCost < costMap(neighbor(1), neighbor(2))
                    costMap(neighbor(1), neighbor(2)) = newCost;
                    parentMap(neighbor(1), neighbor(2)) = sub2ind(gridSize, current(1), current(2));
                    openSet = [openSet; neighbor, newCost]; % Add neighbor to open set
                end
            end
        end
    end
end







%% ///////////////////////////////////
%UI


function [gridSize, robotPositions, mapGrid, radiationMap, radiationGrid, occupancyGrid, logOddsOGM, ax] = initialize_environment()
    % Grid and environment parameters
    gridSize = [30, 30];

    % Create grid and add obstacles
    mapGrid = zeros(gridSize);
    mapGrid(5, 5:25) = 1;
    mapGrid(15, 10:20) = 1;
    mapGrid(25, 1:20) = 1;
    mapGrid(10:25, 10) = 1; %for simulation with hole     mapGrid(10:25, 10) = 1;
    mapGrid(15:25, 20) = 1;
    mapGrid(1:10, 5) = 1;

    % Add random obstacles
    rng(42);
    numRandomObstacles = 50;
    randomPositions = randi([1, gridSize(1)], numRandomObstacles, 2);
    for i = 1:numRandomObstacles
        mapGrid(randomPositions(i, 1), randomPositions(i, 2)) = 1;
    end
    
    mapGrid(16, 13) = 0;
    mapGrid(2, 3) = 0;



    % Initialize log-odds grids
    logOddsOGM = log(0.5 / 0.5) * ones(gridSize); % Neutral log-odds for occupancy
    occupancyGrid = 1 ./ (1 + exp(-logOddsOGM)); % Initial probabilities for occupancy
    radiationMap = nan(gridSize);

    % Generate radiation levels with radioactive objects
    numRadioactiveObjects = 5;
    maxRadiationIntensity = 9;
    [radiationGrid, radioactiveSources] = generate_radiation_map_with_sources(gridSize, numRadioactiveObjects, maxRadiationIntensity);

    % Mark radioactive objects as obstacles
    for i = 1:size(radioactiveSources, 1)
        source = radioactiveSources(i, :);
        mapGrid(source(1), source(2)) = 1; % Mark as obstacles
    end

    % Visualization setup
    figure;

    % Subplot 1: Environment Grid
    ax1 = subplot(1, 3, 1); % Environment grid
    hold on;
    axis equal;
    xlim([0.5, gridSize(2) + 0.5]);
    ylim([0.5, gridSize(1) + 0.5]);
    grid on;
    set(gca, 'XTick', 1:gridSize(2), 'YTick', 1:gridSize(1), ...
             'XTickLabel', '', 'YTickLabel', '', ...
             'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.3);
    title('Environment Grid with Radioactive Objects');

    % Display obstacles and radiation levels
    for row = 1:gridSize(1)
        for col = 1:gridSize(2)
            if mapGrid(row, col) == 1
                rectangle('Position', [col-0.5, row-0.5, 1, 1], 'FaceColor', 'k', 'EdgeColor', 'none');
            else
                text(col, row, num2str(radiationGrid(row, col)), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', 8, 'Color', 'k');
            end
        end
    end

    % Mark starting positions for robots
    

    for i = 1:size(radioactiveSources, 1)
        source = radioactiveSources(i, :);
        plot(source(2), source(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow', ...
             'MarkerEdgeColor', 'black'); % Use a yellow circle as the marker
    end


    % Subplot 2: Radiation OGM
    ax2 = subplot(1, 3, 3); % Radiation probability map
    hold on;
    axis equal;
    xlim([0.5, gridSize(2) + 0.5]);
    ylim([0.5, gridSize(1) + 0.5]);
    grid on;
    set(gca, 'XTick', 1:gridSize(2), 'YTick', 1:gridSize(1), ...
             'XTickLabel', '', 'YTickLabel', '', ...
             'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.3);
    title('Radiation Probability Map');
    colormap(ax2, custom_radiation_colormap());
    c2 = colorbar('Location', 'eastoutside');
    c2.Label.String = 'Radiation Level';
    c2.Ticks = linspace(0, 1, 10);
    c2.TickLabels = 0:9;

    % Subplot 3: Occupancy Grid OGM
    ax3 = subplot(1, 3, 2); % Occupancy probability map
    hold on;
    axis equal;
    xlim([0.5, gridSize(2) + 0.5]);
    ylim([0.5, gridSize(1) + 0.5]);
    grid on;
    set(gca, 'XTick', 1:gridSize(2), 'YTick', 1:gridSize(1), ...
             'XTickLabel', '', 'YTickLabel', '', ...
             'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.3);
    title('Occupancy Probability Map');
    colormap(ax3, flipud(gray));
    c3 = colorbar('Location', 'eastoutside');
    c3.Ticks = linspace(0, 1, 5);
    c3.TickLabels = {'0 (Free)', '0.25', '0.5', '0.75', '1 (Occupied)'};
    % Initialize all cells with 0.5 (unknown state, gray color)
    for row = 1:gridSize(1)
        for col = 1:gridSize(2)
            % Draw a rectangle for each cell, with 0.5 gray color
            rectangle('Position', [col - 0.5, row - 0.5, 1, 1], ...
                      'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');
        end
    end
    ax = [ax1, ax2, ax3];

    % Initial scan

    robotPositions = [1,1; 1,2 ;1,3]; % All robots start at the same position
    plot(ax1,robotPositions(1, 2), robotPositions(1, 1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(ax1,robotPositions(2, 2), robotPositions(2, 1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(ax1,robotPositions(3, 2), robotPositions(3, 1), 'go', 'MarkerSize', 10, 'LineWidth', 2);





end




function [radiationGrid, radioactiveSources] = generate_radiation_map_with_sources(gridSize, numSources, maxIntensity)
    radiationGrid = zeros(gridSize); % Start with a grid of zeros
    rng(42); % For reproducibility

    % Randomly place radioactive sources
    radioactiveSources = randi([1, gridSize(1)], numSources, 2);


    for i = 1:numSources
        source = radioactiveSources(i, :);
        for row = 1:gridSize(1)
            for col = 1:gridSize(2)
                % Calculate Manhattan distance
                distance = abs(row - source(1)) + abs(col - source(2));
    
                if distance == 0
                    % Set maximum intensity for the source cell
                    radiationLevel = inf;
                else
                    % Scale the intensity falloff for other cells
                    effectiveDistance = floor(distance / 2); % Divide by 2 and floor to nearest integer
                    radiationLevel = max(0, maxIntensity - effectiveDistance);
                end
    
                % Assign the radiation level to the cell
                radiationGrid(row, col) = max(radiationGrid(row, col), radiationLevel);
            end
        end
    end



end








function color = get_radiation_color(radiationLevel)
    % Define a custom colormap for radiation levels
    customColormap = [
        1.0, 0.98, 0.8;  % Lemon Chiffon (#FFFACD)
        0.85, 0.65, 0.13;  % Goldenrod (#DAA520)
        0.68, 0.92, 0.63;  % Celery (#ADE792)
        0.83, 0.69, 0.22;  % Metallic Gold (#D4AF37)
        0.99, 0.57, 0.16;  % Sea Buckthorn (#FD9125)
        1.0, 0.51, 0.0;  % Pizazz (#FF8200)
        0.9, 0.17, 0.31;  % Amaranth (#E52B50)
        0.7, 0.13, 0.13;  % Fire Brick (#B22222)
        0.5, 0.09, 0.09;  % Falu Red (#801818)
    ];
    

    
    % Ensure radiationLevel is within the valid range [1, 9]
    normalizedLevel = max(1, min(round(radiationLevel), 9));
    
    % Map the normalized radiation level to the corresponding color
    color = customColormap(normalizedLevel, :);
end


function cmap = custom_radiation_colormap()
    % Define a colormap with intermediate levels
    cmap = [
        1.0, 0.98, 0.8;  % Lemon Chiffon (#FFFACD)
        0.85, 0.85, 0.13;  % Goldenrod (#DAA520)
        0.68, 0.92, 0.63;  % Celery (#ADE792)
        0.90, 0.69, 0.22;  % Metallic Gold (#D4AF37)
        0.99, 0.57, 0.16;  % Sea Buckthorn (#FD9125)
        1.0, 0.51, 0.0;  % Pizazz (#FF8200)
        0.9, 0.17, 0.31;  % Amaranth (#E52B50)
        0.7, 0.13, 0.13;  % Fire Brick (#B22222)
        0.5, 0.09, 0.09;  % Falu Red (#801818)
    ];
end
function grayColor = custom_obstacle_colormap(grayValue)
    % Define 20 distinct grayscale levels for better visual separation
    grayscaleLevels = linspace(0.1, 1, 20); % 20 levels from 0.1 (white) to 1.0 (black)
    % Find the closest level to the input value
    [~, idx] = min(abs(grayscaleLevels - grayValue));
    grayColor = repmat(grayscaleLevels(idx), 1, 3); % Repeat the value across R, G, and B
end


%% Multi robot installation 

function voronoiRegions = compute_voronoi_regions(robotPositions, gridSize, occupancyGrid)
    % Compute dynamic Voronoi regions for robots
    [X, Y] = meshgrid(1:gridSize(2), 1:gridSize(1)); % Create grid coordinates
    numRobots = size(robotPositions, 1); % Number of robots
    distances = zeros(size(X, 1), size(X, 2), numRobots);

    for i = 1:numRobots
        % Calculate distance to robot i; add large penalty for obstacles
        distances(:, :, i) = sqrt((X - robotPositions(i, 2)).^2 + (Y - robotPositions(i, 1)).^2) + ...
                             1e6 * (occupancyGrid > 0.7); % High cost for obstacles
    end
    
    % Assign each cell to the closest robot
    [~, voronoiRegions] = min(distances, [], 3);
    
    % Ensure all cells are assigned to a robot (fallback for unassigned cells)
    unassignedCells = isnan(voronoiRegions);
    if any(unassignedCells(:))
        for row = 1:gridSize(1)
            for col = 1:gridSize(2)
                if unassignedCells(row, col)
                    % Assign the cell to the closest robot manually
                    [~, closestRobot] = min(sqrt((robotPositions(:, 1) - row).^2 + ...
                                                 (robotPositions(:, 2) - col).^2));
                    voronoiRegions(row, col) = closestRobot;
                end
            end
        end
    end


    % Assign each cell to the closest robot
    [~, voronoiRegions] = min(distances, [], 3);
end



function frontiers = detect_frontiers(occupancyGrid, regionMask)
    % Identify frontiers in the specified region
    [rows, cols] = find(abs(occupancyGrid - 0.5) < 1e-2 & regionMask); % Unexplored cells in region
    directions = [0, 1; 1, 0; 0, -1; -1, 0]; % Cardinal directions
    frontiers = [];

    % Check for cells with unexplored neighbors
    for i = 1:length(rows)
        hasUnexploredNeighbor = false; % Track if this cell has any unexplored neighbors
        for d = 1:size(directions, 1)
            neighbor = [rows(i), cols(i)] + directions(d, :);
            if neighbor(1) > 0 && neighbor(1) <= size(occupancyGrid, 1) && ...
               neighbor(2) > 0 && neighbor(2) <= size(occupancyGrid, 2) && ...
               abs(occupancyGrid(neighbor(1), neighbor(2)) - 0.5) < 1e-2
               % Add the cell as a frontier if it has at least one unexplored neighbor
               frontiers = [frontiers; rows(i), cols(i)]; %#ok<AGROW>
               hasUnexploredNeighbor = true;
               break; % Stop checking neighbors for this cell
            end
        end

        % If no unexplored neighbors, consider it an isolated unexplored cell
        if ~hasUnexploredNeighbor
            frontiers = [frontiers; rows(i), cols(i)]; %#ok<AGROW>
        end
    end
end






function explore_grid_multi(gridSize, startPosList, mapGrid, radiationMap, radiationGrid, occupancyGrid, logOddsOGM, ax)
    numRobots = size(startPosList, 1); % Number of robots
    robotPositions = startPosList; % Initial positions
    prevRobotPositions = cell(numRobots, 1); % Track previous positions
    currentRobotPosition = robotPositions;
    costMap = inf(gridSize); % Cost map for Dijkstra-like behavior
    parentMap = zeros(gridSize); % Initialize parentMap with zeros
    costMap(robotPositions(2,1), robotPositions(2,2)) = 0; % Cost to reach the start is zero



    % Initialize Voronoi regions

    % Movement directions (row, column offsets)
    directions = [0, 1; 1, 0; 0, -1; -1, 0];

    % Main exploration loop
    while true
        allExplored = true;
        allUpdatedCells = []; % Collect updated cells for all robots
        voronoiRegions = compute_voronoi_regions(robotPositions, gridSize, occupancyGrid);

        for i = 1:numRobots
            % Update Voronoi regions dynamically
            regionMask = (voronoiRegions == i);

            % Sensor updates for walls and radiation
            [sensorData, visibleCells, costMap, parentMap, radiationMap, occupancyGrid, logOddsOGM] = ...
                sensor_scan(robotPositions(i, :), mapGrid, gridSize, costMap, parentMap, radiationMap, radiationGrid, occupancyGrid, logOddsOGM,directions);
            % Collect updated cells
            updatedCells = [visibleCells; robotPositions(i, :)]; % Include current position
            if ~isempty(prevRobotPositions{i})
                updatedCells = [updatedCells; prevRobotPositions{i}]; % Include previous position
            end
            allUpdatedCells = [allUpdatedCells; updatedCells]; %#ok<AGROW>

            % Detect frontiers within the robot's Voronoi region
            frontiers = detect_frontiers(occupancyGrid, regionMask);
            if isempty(frontiers)
                continue;
            end

            % Select the closest frontier within the Voronoi region

            % Plan path to the selected frontier using BFS
            % disp(frontiers)
            % disp(i)
            nextMove = bfs_to_frontier(robotPositions(i, :), frontiers, occupancyGrid);

            if isempty(nextMove)
                disp(['Robot ' num2str(i) ' cannot reach any frontiers.']);
                continue;
            end

            % Update robot position (move one step along the path)
            prevRobotPositions{i} = robotPositions(i, :);
            currentRobotPosition(i, :) = robotPositions(i, :);
            robotPositions(i, :) = nextMove;

            % Mark that at least one robot is still exploring
            allExplored = false;
        end

        % Draw all updated cells for all robots at once
        if ~isempty(allUpdatedCells)
            draw_grid_update_multi(currentRobotPosition, radiationMap, allUpdatedCells, occupancyGrid, ax);
        end
        draw_robot_multi(robotPositions, prevRobotPositions, radiationMap, occupancyGrid, ax);

        % Stop the loop if all robots have completed their exploration
        if allExplored
            disp('All robots have completed their exploration.');
            break;
        end
    end

    % Final sensor scan for completeness
    for i = 1:numRobots
        [sensorData, visibleCells, costMap, parentMap, radiationMap, occupancyGrid, logOddsOGM] = ...
            sensor_scan(robotPositions(i, :), mapGrid, gridSize, costMap, parentMap, radiationMap, radiationGrid, occupancyGrid, logOddsOGM,directions);
        allUpdatedCells = [allUpdatedCells; visibleCells]; %#ok<AGROW>
    end

    % Draw the final updates
    if ~isempty(allUpdatedCells)
        draw_grid_update_multi([], radiationMap, allUpdatedCells, occupancyGrid, ax);
    end
    draw_robot_multi(robotPositions, prevRobotPositions, radiationMap, occupancyGrid, ax);


    % Return robots to their starting positions
    disp('Returning robots to their start positions...');
    % Precompute the paths for all robots
    pathsToStart = compute_return_paths(gridSize, startPosList, robotPositions, occupancyGrid);

    
    % Find the maximum path length among all robots

    maxSteps = max(cellfun(@(path) size(path, 1), pathsToStart));
    
    % Initialize the robots' positions along their paths
    for step = 1:max(cellfun(@length, pathsToStart))
        for i = 1:numRobots
            if step <= length(pathsToStart{i})
                prevRobotPositions{i} = robotPositions(i, :);
                robotPositions(i, :) = pathsToStart{i}(step, :);
            end
        end
    
        % Draw all robots at their updated positions
        draw_robot_multi(robotPositions, prevRobotPositions, radiationMap, occupancyGrid, ax);

    
    end


    validate_exploration_results(gridSize,mapGrid,occupancyGrid, logOddsOGM)

    disp('All robots returned to their start positions.');
end


function path = bfs_to_frontier(robotPos, frontiers, occupancyGrid)
    % Perform BFS to find the nearest reachable frontier
    gridSize = size(occupancyGrid);
    queue = {robotPos}; % Start from the robot's position
    parentMap = containers.Map('KeyType', 'char', 'ValueType', 'any'); % Track paths
    parentMap(mat2str(robotPos)) = NaN; % Mark the start position

    % Convert frontiers to a set for quick lookup
    frontierSet = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    for i = 1:size(frontiers, 1)
        frontierSet(mat2str(frontiers(i, :))) = true;
    end

    % Perform BFS
    while ~isempty(queue)
        current = queue{1};
        queue(1) = [];

        % Check if the current cell is a frontier
        if isKey(frontierSet, mat2str(current))
            % Reconstruct path to this frontier
            % if parentMap.Count == 1
            %     show_occupancy_grid_values(voronoi);
            % end
            path = construct_bfs_path(parentMap, current, robotPos);
            return;
        end

        % Explore neighbors
        for d = [0, 1; 1, 0; 0, -1; -1, 0]' % Cardinal directions
            neighbor = current + d';
            if neighbor(1) > 0 && neighbor(1) <= gridSize(1) && ...
               neighbor(2) > 0 && neighbor(2) <= gridSize(2) && ...
               occupancyGrid(neighbor(1), neighbor(2)) < 0.9 && ... % Not a wall
               ~isKey(parentMap, mat2str(neighbor))
                queue{end + 1} = neighbor; % Add to queue
                parentMap(mat2str(neighbor)) = current; % Mark as visited
            end
        end
    end

    % If no frontier is reachable
    path = [];
    disp('No reachable frontier found.');
end




function draw_grid_update_multi(robotPositions, radiationMap, updatedCells, occupancyGrid, ax)

    % Draw only the updated cells
    for i = 1:size(updatedCells, 1)
        row = updatedCells(i, 1);
        col = updatedCells(i, 2);

        occupied = false;
        for e = 1:size(robotPositions, 1) % Iterate over rows of robotPositions
            if row == robotPositions(e, 1) && col == robotPositions(e, 2) % Check both row and column
                occupied = true;
                break; % Exit the loop early if a match is found
            end
        end
        
        if occupied
            continue; % Skip the current iteration
        end

        % Update based on occupancyGrid
        grayValue = 1 - occupancyGrid(row, col); % Map to grayscale
        grey_color = custom_obstacle_colormap(grayValue);
        rectangle(ax(3), 'Position', [col-0.5, row-0.5, 1, 1], 'FaceColor', grey_color, 'EdgeColor', 'none');
        if occupancyGrid(row, col) > 0.8
            % High occupancy probability (obstacle)
            radiationLevel = radiationMap(row, col);
            rectangle(ax(2), 'Position', [col-0.5, row-0.5, 1, 1], 'FaceColor', 'k', 'EdgeColor', 'none');
            if radiationLevel == inf
                plot(ax(2), col, row, 'x', 'MarkerSize', 10, 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'white');

            end
            
            % Map occupancy probability to a grayscale value for free or unexplored cells

        % Update based on radiationMap
        elseif ~isnan(radiationMap(row, col))
            % Radiation level detected
            radiationLevel = radiationMap(row, col);
            color = get_radiation_color(radiationLevel); % Map radiation level to a color
            rectangle(ax(2), 'Position', [col-0.5, row-0.5, 1, 1], 'FaceColor', color, 'EdgeColor', 'none');
            
        end
    end

    % Update the visualization
    drawnow;
end


function draw_robot_multi(robotPositions, prevRobotPositions, radiationMap, occupancyGrid, ax)



    % Restore the previous robot positions
    for i = 1:size(prevRobotPositions, 1)
        if ~isempty(prevRobotPositions{i})
            row = prevRobotPositions{i}(1);
            col = prevRobotPositions{i}(2);

            % Restore previous cell based on occupancyGrid
            if ~isnan(radiationMap(row, col))
                radiationLevel = radiationMap(row, col);
                color = get_radiation_color(radiationLevel);
                rectangle(ax(2), 'Position', [col-0.5, row-0.5, 1, 1], ...
                          'FaceColor', color, 'EdgeColor', 'none');
            end

            grayValue = 1 - occupancyGrid(row, col); % Grayscale for occupancy
            rectangle(ax(3), 'Position', [col-0.5, row-0.5, 1, 1], ...
                      'FaceColor', [grayValue, grayValue, grayValue], 'EdgeColor', 'none');
        end
    end

    % Draw the current robot positions
    for i = 1:size(robotPositions, 1)
        row = robotPositions(i, 1);
        col = robotPositions(i, 2);

        % Draw robot positions on both layers
        rectangle(ax(2), 'Position', [col-0.5, row-0.5, 1, 1], ...
                  'FaceColor', [0, 1, 1], 'EdgeColor', 'k');
        rectangle(ax(3), 'Position', [col-0.5, row-0.5, 1, 1], ...
                  'FaceColor', [0, 1, 1], 'EdgeColor', 'k');
        text(ax(2),col, row, string(i), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', 8, 'Color', 'k');
        text(ax(3),col, row, string(i), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', 8, 'Color', 'k');
    end


    drawnow;
end

function validate_exploration_results(gridSize, trueMap, occupancyGrid, logOddsOGM)
    % Metrics for validation
    numCells = numel(trueMap);
    correctPredictions = sum((trueMap(:) == 1 & occupancyGrid(:) > 0.7) | (trueMap(:) == 0 & occupancyGrid(:) < 0.3));
    falsePositives = sum(trueMap(:) == 0 & occupancyGrid(:) > 0.7);
    falseNegatives = sum(trueMap(:) == 1 & occupancyGrid(:) < 0.3);
    
    % Compute accuracy, FPR, and FNR
    accuracy = correctPredictions / numCells * 100;
    fpr = falsePositives / nnz(trueMap(:) == 0) * 100;
    fnr = falseNegatives / nnz(trueMap(:) == 1) * 100;
    
    % Display results
    disp('Exploration Validation Results:');
    fprintf('Accuracy: %.2f%%\n', accuracy);
    fprintf('False Positive Rate (FPR): %.2f%%\n', fpr);
    fprintf('False Negative Rate (FNR): %.2f%%\n', fnr);

end


function show_occupancy_grid_values(occupancyGrid)
    % Function to display the occupancyGrid with real values (no color)
    %
    % Parameters:
    %   occupancyGrid - The grid to be displayed (values between 0 and 1)
    %   titleText - Optional title for the figure

    % Create a figure
    figure;
    
    % Get the size of the grid
    [rows, cols] = size(occupancyGrid);
    
    % Create a table-like visualization using a grid of text
    hold on;
    axis equal;
    xlim([0.5, cols + 0.5]);
    ylim([0.5, rows + 0.5]);
    set(gca, 'YDir', 'reverse'); % Reverse y-axis to display like a matrix

    % Draw the grid lines
    for i = 0.5:1:rows
        plot([0.5, cols + 0.5], [i, i], 'k'); % Horizontal lines
    end
    for j = 0.5:1:cols
        plot([j, j], [0.5, rows + 0.5], 'k'); % Vertical lines
    end

    % Write the values inside the grid cells
    for i = rows:-1:1 % Reverse row order
        for j = 1:cols
            text(j, rows - i + 1, sprintf('%.2f', occupancyGrid(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        end
    end

    % Add optional title


    % Set axis labels
    xlabel('X (Grid Columns)');
    ylabel('Y (Grid Rows)');
    set(gca, 'XTick', 1:cols, 'YTick', 1:rows); % Align ticks with grid cells
    grid on; % Enable gridlines for better cell separation
    set(gca, 'GridColor', 'k'); % Set grid color to black
end
