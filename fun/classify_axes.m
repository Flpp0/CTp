function axes = classify_axes(coeff)
% classify_axes Classifies PCA principal components into anatomical axes.
%
% Syntax:
%   axes = classify_axes(coeff)
%
% Description:
%   This function takes a 3x3 PCA coefficient matrix (obtained from PCA performed on brain
%   voxel coordinates) and assigns each principal component to an anatomical axis (Left-Right,
%   Front-Back, Top-Bottom) based on maximum absolute alignment with the corresponding coordinate.
%
% Input:
%   coeff - A 3x3 matrix where each column is a principal component.
%
% Output:
%   axes  - A structure with fields:
%             LeftRight  - Principal component assigned to the Left-Right axis.
%             FrontBack  - Principal component assigned to the Front-Back axis.
%             TopBottom  - Principal component assigned to the Top-Bottom axis.
%
% Example:
%   coeff = [0.8 0.1 0.6; 0.2 0.9 0.1; 0.5 0.3 0.8];
%   axes = classify_axes(coeff);

    % Heuristic way of classifying each principal component to an
    % anatomical axes. 

    % Absolute values
    absCoeff = abs(coeff);
    
    % Anatomical Axes and corresponding directions
    anatomicalAxes = {'LeftRight', 'FrontBack', 'TopBottom'};
    anatomicalDirections = {'X', 'Y', 'Z'};
    
    axes = struct();
    assignedPCs = false(1, 3); % Keep track if the component already selected or not
    
    for i = 1:3
        axisName = anatomicalAxes{i};
        direction = anatomicalDirections{i};
        
        % Determine the corresponding row index based on the coordinate direction
        switch direction
            case 'X'
                axisIdx = 1;
            case 'Y'
                axisIdx = 2;
            case 'Z'
                axisIdx = 3;
        end
        
        % Find the principal component with the highest absolute coefficient in the chosen direction
        [~, maxPC] = max(absCoeff(axisIdx, :));
        
        % If this component has already been assigned, choose the next highest one
        if assignedPCs(maxPC)
            [~, sortedPCs] = sort(absCoeff(axisIdx, :), 'descend');
            for j = 1:length(sortedPCs)
                if ~assignedPCs(sortedPCs(j))
                    maxPC = sortedPCs(j);
                    break;
                end
            end
        end
        
        % Assign the selected principal component to the current anatomical axis
        axes.(axisName) = coeff(:, maxPC);
        assignedPCs(maxPC) = true;
    end
end
