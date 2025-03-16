function plotRegistrationParameters(translations_mm, rotations, time_values, highMovementsTable, baseDir, patientCode, upsample_factor)
% plotRegistrationParameters Plots translation and rotation parameters over time.
%
% Syntax:
%   plotRegistrationParameters(translations_mm, rotations, time_values, highMovementsTable, baseDir, patientCode, upsample_factor)
%
% Description:
%   This function creates plots for the translation (in mm) and Euler angle rotation (in degrees)
%   parameters extracted over time. It highlights time points with high movements. 
%
% Inputs:
%   translations_mm    - Array of translations in millimeters.
%   rotations          - Array of Euler angle rotations in degrees.
%   time_values        - Column vector of time indices.
%   highMovementsTable - Table containing indices and movement magnitudes for high movements.
%   baseDir            - Base directory where patient folders are located.
%   patientCode        - Patient code (folder name).
%   upsample_factor    - Upsampling factor used in processing.
%
% Example:
%   plotRegistrationParameters(translations_mm, rotations, time_values, highMovementsTable, baseDir, patientCode, upsample_factor);

    fontSizeValue = 14;

    % Plot Translations
    fig_translations = figure('Name', 'Translations over Time (mm)', 'Position', [100, 100, 800, 600]);
    subplot(3,1,1);
    plot(time_values, translations_mm(:,1), '-o', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', fontSizeValue); 
    ylabel('Tx (mm)', 'FontSize', fontSizeValue);
    title('Translation X over Time', 'FontSize', fontSizeValue);
    grid on;
    set(gca, 'FontSize', fontSizeValue);

    subplot(3,1,2);
    plot(time_values, translations_mm(:,2), '-o', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', fontSizeValue);
    ylabel('Ty (mm)', 'FontSize', fontSizeValue);
    title('Translation Y over Time', 'FontSize', fontSizeValue);
    grid on;
    set(gca, 'FontSize', fontSizeValue);

    subplot(3,1,3);
    plot(time_values, translations_mm(:,3), '-o', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', fontSizeValue);
    ylabel('Tz (mm)', 'FontSize', fontSizeValue);
    title('Translation Z over Time', 'FontSize', fontSizeValue);
    grid on;
    set(gca, 'FontSize', fontSizeValue);

    % Highlight high movement time points for translations
    hold on;
    for i = 1:3
        subplot(3,1,i);
        hold on;
        plot(time_values(highMovementsTable.Index), translations_mm(highMovementsTable.Index, i), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
    hold off;

    % Plot Rotations
    fig_rotations = figure('Name', 'Rotations over Time (degrees)', 'Position', [100, 100, 800, 600]);
    subplot(3,1,1);
    plot(time_values, rotations(:,1), '-o', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', fontSizeValue);
    ylabel('Rotation X (deg)', 'FontSize', fontSizeValue);
    title('Rotation around X-axis', 'FontSize', fontSizeValue);
    grid on;
    set(gca, 'FontSize', fontSizeValue);

    subplot(3,1,2);
    plot(time_values, rotations(:,2), '-o', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', fontSizeValue);
    ylabel('Rotation Y (deg)', 'FontSize', fontSizeValue);
    title('Rotation around Y-axis', 'FontSize', fontSizeValue);
    grid on;
    set(gca, 'FontSize', fontSizeValue);

    subplot(3,1,3);
    plot(time_values, rotations(:,3), '-o', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontSize', fontSizeValue);
    ylabel('Rotation Z (deg)', 'FontSize', fontSizeValue);
    title('Rotation around Z-axis', 'FontSize', fontSizeValue);
    grid on;
    set(gca, 'FontSize', fontSizeValue);

    % Highlight high movement time points for rotations
    hold on;
    for i = 1:3
        subplot(3,1,i);
        hold on;
        plot(time_values(highMovementsTable.Index), rotations(highMovementsTable.Index, i), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
    hold off;

    % Save the Figures in PDF and high-resolution PNG formats
    paramsSaveDir = fullfile(baseDir, patientCode, 'Registration', '3D', 'Mutual Information', ['Upsample_', num2str(upsample_factor)], 'Registration Parameters');
    if ~exist(paramsSaveDir, 'dir')
        mkdir(paramsSaveDir);
    end
    
    % PDF
    exportgraphics(fig_rotations, fullfile(paramsSaveDir, 'RotationsOverTime.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');
    exportgraphics(fig_translations, fullfile(paramsSaveDir, 'TranslationsOverTime.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');
    % PNG
    exportgraphics(fig_rotations, fullfile(paramsSaveDir, 'RotationsOverTime.png'), 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);
    exportgraphics(fig_translations, fullfile(paramsSaveDir, 'TranslationsOverTime.png'), 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);
    
    % Save Registration Parameters to a MAT file
    save(fullfile(paramsSaveDir, 'registrationParameters.mat'), 'translations_mm', 'rotations', 'time_values', 'highMovementsTable');
end
