%% COMMON PCA ANALYSIS (DATA AND MODEL)
%
% This script fits a single Principal Component Analysis (PCA) model across
% four experimental recordings or four conditions from an Izhikevich model.
% By calculating the PCA on the concatenated data, all conditions share the
% same component space, allowing for direct trajectory comparisons.
%
% Generates:
%   1) A combined 3D figure showing all four trajectories in the same space.
%   2) Individual 2D and 3D views for each condition.
%   3) A MAT file with the numerical PCA results (coefficients, scores, etc.).

clear; 
clc; 
close all;

%% 1. PATHS AND DIRECTORIES
project_folder = pwd;
addpath(fullfile(project_folder, 'Functions'), '-end');

%% 2. SELECT DATA OR MODEL
% ================================================================
% LEAVE ONLY ONE OF THE FOLLOWING TWO LINES ACTIVE:
% mode_type = 'data';
mode_type = 'model';
% ================================================================
% Experimental Data File
data_path = "C:\Users\javie\Desktop\Data\Tesis_Data\Sep0622_serotonin.mat";

% Model File
model_path =  fullfile(project_folder, 'Grid_Results', 'Models', 'Modelo_IZH_seed_009_perm_002.mat');
% ================================================================

%% 3. COMMON CONFIGURATION
dt = 1e-3;                 % Time step [s]
n_conditions = 4;
use_global_zscore = false; % Standardize before PCA
max_plot_points = 4000;    % Downsample for plotting to avoid massive files

% Figure export settings controlled by auxiliary function
save_png = true;
save_svg = true;
export_formats = {'pdf'};
if save_png, export_formats{end+1} = 'png'; end
if save_svg, export_formats{end+1} = 'svg'; end

close_figures = false;

% Smoothing (Gaussian Kernel standard deviation in seconds)
smooth_data = 5;
smooth_model = 2.5;

% Figure sizes [width, height] in cm
size_combined = [30 15];
size_individual = [15 10];

% Custom color gradients for the 4 conditions (Light to Dark)
c_starts = [0.80 0.90 1.00; 1.00 0.90 0.95; 0.95 0.85 1.00; 0.80 1.00 0.80];
c_ends   = [0.00 0.20 0.70; 0.85 0.10 0.45; 0.45 0.05 0.65; 0.00 0.50 0.15];
cmaps = cell(4, 1);
for i = 1:4
    cmaps{i} = [linspace(c_starts(i,1), c_ends(i,1), 256)', ...
                linspace(c_starts(i,2), c_ends(i,2), 256)', ...
                linspace(c_starts(i,3), c_ends(i,3), 256)'];
end

%% 4. SOURCE CONFIGURATION
switch lower(mode_type)
    case 'data'
        input_path = data_path;
        min_time = 605;
        max_time = 1005;
        smooth_sigma_s = smooth_data;
        cond_names = {'No serotonin'; 'Serotonin injection'; ...
                      'Maintained serotonin'; 'Washout, no serotonin'};
        source_name = 'Experimental Data';
        results_folder = fullfile(project_folder, 'Common_PCA_Results', 'Data');
    case 'model'
        input_path = model_path;
        min_time = 30;
        max_time = 270;
        smooth_sigma_s = smooth_model;
        [~, model_file_name] = fileparts(model_path);
        model_file_name = regexprep(model_file_name, '[^a-zA-Z0-9_-]', '_');
        source_name = strrep(model_file_name, '_', ' ');
        results_folder = fullfile(project_folder, 'Common_PCA_Results', 'Model', model_file_name);
    otherwise
        error("mode_type must be 'data' or 'model'.");
end

if ~isfile(input_path), error('File not found: %s', input_path); end
if ~isfolder(results_folder), mkdir(results_folder); end

%% 5. LOAD AND VALIDATE FILE
data = load(input_path);
for g = 1:n_conditions
    if ~isfield(data, sprintf('spks_%d', g))
        error('File is missing variable spks_%d.', g);
    end
end

% Extract Model condition names
combinations = nan(n_conditions, 3);
if strcmpi(mode_type, 'model')
    if isfield(data, 'modelo') && isfield(data.modelo, 'combinaciones')
        combinations = double(data.modelo.combinaciones);
        cond_names = cell(n_conditions, 1);
        for g = 1:n_conditions
            cond_names{g} = sprintf('C%d: a=%.4g, b=%.4g, d=%.4g', ...
                g, combinations(g,1), combinations(g,2), combinations(g,3));
        end
    else
        cond_names = cellstr(compose('Condition %d', (1:n_conditions)'));
    end
end

% Determine total neurons
expected_neurons = 0;
if isfield(data, 'x') && ~isempty(data.x)
    expected_neurons = numel(data.x);
else
    for g = 1:n_conditions
        spks = data.(sprintf('spks_%d', g));
        if ~isempty(spks), expected_neurons = max(expected_neurons, max(round(spks(:,1)))); end
    end
end
if expected_neurons < 1, error('Could not determine the number of neurons.'); end

%% 6. BUILD SMOOTHED ACTIVITY MATRICES (Optimized Native Smoothing)
X_list = cell(n_conditions, 1);
time_list = cell(n_conditions, 1);
n_samples = zeros(n_conditions, 1);

time_edges = min_time:dt:max_time;
n_bins = length(time_edges) - 1;
t_vector = time_edges(1:end-1)' + dt/2;

% Pre-compute Gaussian Kernel
sigma_bins = smooth_sigma_s / dt;
radius = ceil(4 * sigma_bins);
kernel = exp(-0.5 * ((-radius:radius)' / sigma_bins).^2);
kernel = kernel / sum(kernel); % Normalize

for g = 1:n_conditions
    spks = data.(sprintf('spks_%d', g));
    
    % Ensure correct orientation
    if ~isempty(spks)
        spks = double(spks);
        if size(spks, 2) ~= 2 && size(spks, 1) == 2, spks = spks.'; end
    end
    
    % Fast binning
    counts = zeros(n_bins, expected_neurons);
    if ~isempty(spks)
        valid_spks = spks(spks(:,2) >= min_time & spks(:,2) <= max_time, :);
        for n = 1:expected_neurons
            neuron_times = valid_spks(round(valid_spks(:,1)) == n, 2);
            counts(:, n) = histcounts(neuron_times, time_edges);
        end
    end
    
    % Gaussian Convolution
    activity = conv2(counts, kernel, 'same') / dt; % Convert to Hz
    
    X_list{g} = activity;
    n_samples(g) = n_bins;
    time_list{g} = t_vector;
end

%% 7. FIT COMMON PCA
X_combined = vertcat(X_list{:});
global_mean = [];
global_std = [];

if use_global_zscore
    global_mean = mean(X_combined, 1, 'omitnan');
    global_std = std(X_combined, 0, 1, 'omitnan');
    global_std(global_std == 0 | isnan(global_std)) = 1;
    X_combined = (X_combined - global_mean) ./ global_std;
end

[coeff, score_combined, latent, ~, explained, mu] = pca(X_combined);
if size(score_combined, 2) < 3
    error('PCA generated less than 3 principal components.');
end

fprintf('\nPCA (%s):\n', source_name);
fprintf('PC1 = %.2f %%\n', explained(1));
fprintf('PC2 = %.2f %%\n', explained(2));
fprintf('PC3 = %.2f %%\n', explained(3));
fprintf('PC1 + PC2 + PC3 = %.2f %%\n\n', sum(explained(1:3)));

% Split scores back into conditions
score_conds = cell(n_conditions, 1);
start_idx = 1;
for g = 1:n_conditions
    end_idx = start_idx + n_samples(g) - 1;
    score_conds{g} = score_combined(start_idx:end_idx, :);
    start_idx = end_idx + 1;
end

% Determine Common Axes Limits
S_all = score_combined(:, 1:3);
ax_limits = zeros(3, 2);
for pc = 1:3
    margin = 0.05 * (max(S_all(:,pc)) - min(S_all(:,pc)));
    if margin == 0, margin = 1; end
    ax_limits(pc, :) = [min(S_all(:,pc)) - margin, max(S_all(:,pc)) + margin];
end

%% 8. COMBINED 3D FIGURE
fig_combined = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 size_combined], 'Renderer', 'painters', 'Visible', 'off');
ax = axes(fig_combined, 'Position', [0.08 0.12 0.64 0.80]);
hold(ax, 'on'); box(ax, 'on');

for g = 1:n_conditions
    % Downsample for performance
    S = score_conds{g}(:, 1:3);
    if size(S, 1) > max_plot_points
        idx = unique(round(linspace(1, size(S, 1), max_plot_points)));
    else
        idx = 1:size(S, 1);
    end
    S_plot = S(idx, :);
    cmap = cmaps{g};
    
    % Draw colored trajectory using surface hack
    n_pts = size(S_plot, 1);
    color_idx = round(linspace(1, size(cmap, 1), n_pts));
    rgb = cmap(color_idx, :);
    C = zeros(n_pts, 2, 3);
    C(:, 1, :) = reshape(rgb, n_pts, 1, 3);
    C(:, 2, :) = reshape(rgb, n_pts, 1, 3);
    
    surface(ax, [S_plot(:,1) S_plot(:,1)], [S_plot(:,2) S_plot(:,2)], [S_plot(:,3) S_plot(:,3)], C, ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.0);
        
    % Mark Start Point
    plot3(ax, S_plot(1,1), S_plot(1,2), S_plot(1,3), 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', cmap(1,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
end

xlim(ax, ax_limits(1, :)); ylim(ax, ax_limits(2, :)); zlim(ax, ax_limits(3, :));
view(ax, 3); axis(ax, 'vis3d');
xlabel(ax, sprintf('PC1 (%.1f %%)', explained(1)));
ylabel(ax, sprintf('PC2 (%.1f %%)', explained(2)));
zlabel(ax, sprintf('PC3 (%.1f %%)', explained(3)));
set(ax, 'FontSize', 9, 'LineWidth', 0.8, 'TickDir', 'out');

% Add mini colorbars for each condition
bar_positions = [0.78 0.72 0.035 0.16; 0.78 0.51 0.035 0.16; ...
                 0.78 0.30 0.035 0.16; 0.78 0.09 0.035 0.16];
for g = 1:n_conditions
    ax_bar = axes(fig_combined, 'Position', bar_positions(g, :));
    img = reshape(linspace(0, 1, 256), [], 1);
    t_range = [time_list{g}(1) time_list{g}(end)];
    imagesc(ax_bar, [0 1], t_range, img);
    colormap(ax_bar, cmaps{g});
    set(ax_bar, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right', ...
        'YTick', t_range, 'FontSize', 8, 'Box', 'on');
    ylabel(ax_bar, cond_names{g}, 'Rotation', 0, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', 'Interpreter', 'none');
end
annotation(fig_combined, 'textbox', [0.75 0.91 0.23 0.035], 'String', 'Color = Time', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Save Combined Figure using export_figure
base_name_combined = fullfile(results_folder, 'Common_PCA_3D_Combined');
export_figure(fig_combined, base_name_combined, size_combined(1), size_combined(2), export_formats, 600);
close(fig_combined);

%% 9. INDIVIDUAL FIGURES (2D and 3D)
for g = 1:n_conditions
    S = score_conds{g}(:, 1:3);
    t = time_list{g};
    
    if size(S, 1) > max_plot_points
        idx = unique(round(linspace(1, size(S, 1), max_plot_points)));
    else
        idx = 1:size(S, 1);
    end
    S_plot = S(idx, :); t_plot = t(idx);
    cmap = cmaps{g};
    
    fig_ind = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 size_individual], 'Renderer', 'painters', 'Visible', 'off');
    ax_ind = axes(fig_ind, 'Position', [0.11 0.12 0.70 0.81]);
    hold(ax_ind, 'on'); box(ax_ind, 'on');
    
    % Surface Hack for colored line
    n_pts = size(S_plot, 1);
    color_idx = round(linspace(1, size(cmap, 1), n_pts));
    rgb = cmap(color_idx, :);
    C = zeros(n_pts, 2, 3);
    C(:, 1, :) = reshape(rgb, n_pts, 1, 3); C(:, 2, :) = reshape(rgb, n_pts, 1, 3);
    
    surface(ax_ind, [S_plot(:,1) S_plot(:,1)], [S_plot(:,2) S_plot(:,2)], [S_plot(:,3) S_plot(:,3)], C, ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.0);
        
    % Start point
    plot3(ax_ind, S_plot(1,1), S_plot(1,2), S_plot(1,3), 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', cmap(1,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        
    xlim(ax_ind, ax_limits(1, :)); ylim(ax_ind, ax_limits(2, :)); zlim(ax_ind, ax_limits(3, :));
    xlabel(ax_ind, sprintf('PC1 (%.1f %%)', explained(1)));
    ylabel(ax_ind, sprintf('PC2 (%.1f %%)', explained(2)));
    zlabel(ax_ind, sprintf('PC3 (%.1f %%)', explained(3)));
    set(ax_ind, 'FontSize', 9, 'LineWidth', 0.8, 'TickDir', 'out');
    
    colormap(ax_ind, cmap); clim(ax_ind, [t_plot(1) t_plot(end)]);
    cb = colorbar(ax_ind); cb.Label.String = 'Time (s)';
    
    % --- Save 2D View ---
    view(ax_ind, 2); set(ax_ind, 'DataAspectRatioMode', 'auto', 'PlotBoxAspectRatioMode', 'auto');
    drawnow;
    base_name_2d = fullfile(results_folder, sprintf('PCA_2D_Condition_%d', g));
    export_figure(fig_ind, base_name_2d, size_individual(1), size_individual(2), export_formats, 600);
    
    % --- Save 3D View ---
    view(ax_ind, 3); axis(ax_ind, 'vis3d');
    drawnow;
    base_name_3d = fullfile(results_folder, sprintf('PCA_3D_Condition_%d', g));
    export_figure(fig_ind, base_name_3d, size_individual(1), size_individual(2), export_formats, 600);
    
    close(fig_ind);
end

%% 10. SAVE NUMERICAL RESULTS
mat_file = fullfile(results_folder, 'Common_PCA_Results.mat');
save(mat_file, 'coeff', 'score_conds', 'latent', 'explained', 'mu', 'time_list', ...
    'cond_names', 'combinations', 'mode_type', 'input_path', ...
    'min_time', 'max_time', 'dt', 'smooth_sigma_s', 'expected_neurons', ...
    'use_global_zscore', 'global_mean', 'global_std', 'ax_limits', '-v7.3');
    
fprintf('Numerical results saved in:\n%s\n', results_folder);