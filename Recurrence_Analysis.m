%% DATA AND MODEL RECURRENCE ANALYSIS
% Analyzes experimental recordings or model conditions separately.
%
% For each condition:
%   1. Calculates the PCA trajectory of smoothed neural activity.
%   2. Detects future recurrence using a distance matrix.
%   3. Calculates the recurrence percentage using a moving window.
%   4. Detects the first sustained coalescence time.
%   5. Saves an independent figure and CSV table.
%
% Additionally, it displays a summary table in the command window.

%% DATA AND MODEL RECURRENCE ANALYSIS
% Analyzes experimental recordings or model conditions separately.
%
% For each condition:
%   1. Calculates the PCA trajectory of smoothed neural activity.
%   2. Detects future recurrence using a distance matrix.
%   3. Calculates the recurrence percentage using a moving window.
%   4. Detects the first sustained coalescence time.
%   5. Saves an independent figure and CSV table.

clear; 
clc; 
close all;

%% 1. PATHS AND DIRECTORIES
project_folder = pwd;
addpath(fullfile(project_folder, 'Functions'), '-end');

%% 2. SELECT DATA OR MODEL
% ================================================================
% KEEP ONLY ONE OF THESE TWO LINES ACTIVE:
mode_type = 'data';
% mode_type = 'model';
% ================================================================
data_path = "C:\Users\javie\Desktop\Data\Tesis_Data\Sep0622_serotonin.mat";
model_path =  fullfile(project_folder, 'Grid_Results', 'Models', 'Modelo_IZH_seed_009_perm_002.mat');

%% 3. ANALYSIS SETTINGS
language = 'EN';
n_conditions = 4;

% Neural activity settings
bin_s = 0.05;              % Time resolution: 20 points per second
smooth_sigma_s = 5;        % Gaussian kernel standard deviation [s]
num_components = 3;        % PCA components used for distances

% Recurrence and moving window settings
epsilon_percentile = 10;   % Percentile to define the recurrence threshold
exclusion_s = 5;           % Time to exclude before checking future recurrence [s]
window_s = 5;              % Moving window size [s]
step_s = 1;                % Moving window step [s]
coal_threshold = 0.90;     % Coalescence threshold (90%)
consecutive_windows = 3;   % Windows needed to confirm coalescence

% Export settings
save_figures = true;
save_results = true;
fig_width_cm = 15;
fig_height_cm = 8;

export_formats = {'pdf', 'png', 'svg'}; % Controlled by auxiliary function

%% 4. SOURCE CONFIGURATION
switch lower(mode_type)
    case 'data'
        input_path = data_path;
        min_time = 600;
        max_time = 900;
        condition_names = {'Recording 1'; 'Recording 2'; 'Recording 3'; 'Recording 4'};
        file_prefix = 'data';
        output_folder = fullfile(project_folder, 'Recurrence_Results', 'Data');
    case 'model'
        input_path = model_path;
        min_time = 10;
        max_time = 310;
        [~, model_name] = fileparts(model_path);
        model_name = regexprep(model_name, '[^a-zA-Z0-9_-]', '_');
        file_prefix = model_name;
        output_folder = fullfile(project_folder, 'Recurrence_Results', 'Model', model_name);
    otherwise
        error("Mode must be 'data' or 'model'.");
end

if ~isfile(input_path)
    error('File not found: %s', input_path);
end
if (save_figures || save_results) && ~isfolder(output_folder)
    mkdir(output_folder);
end

data = load(input_path);
for k = 1:n_conditions
    if ~isfield(data, sprintf('spks_%d', k))
        error('File is missing variable spks_%d.', k);
    end
end

% Retrieve model condition names if applicable
if strcmpi(mode_type, 'model') && isfield(data, 'modelo') && isfield(data.modelo, 'combinaciones')
    comb = double(data.modelo.combinaciones);
    condition_names = cellstr(compose('C%d: a=%.4g, b=%.4g, d=%.4g', (1:n_conditions)', comb));
end

% Determine expected total neurons
num_neurons_expected = 0;
if isfield(data, 'x') && isnumeric(data.x), num_neurons_expected = numel(data.x); end

%% 5. INITIALIZE SUMMARY VARIABLES
coal_time_s = nan(n_conditions, 1);
coal_time_rel_s = nan(n_conditions, 1);
epsilon_array = nan(n_conditions, 1);
num_neurons_array = zeros(n_conditions, 1);
num_active_array = zeros(n_conditions, 1);
num_spikes_array = zeros(n_conditions, 1);
num_windows_array = zeros(n_conditions, 1);

%% 6. PROCESS EACH CONDITION
for cond_idx = 1:n_conditions
    fprintf('\nAnalyzing condition %d of %d...\n', cond_idx, n_conditions);
    
    [PS, spks_raw] = Plot_Settings(cond_idx, data, language);
    main_color = validatecolor(PS.color_raster);
    
    % USE AUXILIARY FUNCTION: format_spikes
    [valid_spikes, n_neurons] = format_spikes(spks_raw, min_time, max_time, num_neurons_expected);
    if num_neurons_expected == 0, num_neurons_expected = n_neurons; end
    
    if isempty(valid_spikes)
        warning('No spikes in window for condition %d.', cond_idx); 
        continue; 
    end
    
    % Convert spikes to smoothed rates
    n_bins = floor((max_time - min_time) / bin_s);
    edges = min_time + (0:n_bins) * bin_s;
    time_vector = edges(1:end-1).' + bin_s / 2;
    
    counts = zeros(n_bins, n_neurons);
    for id = 1:n_neurons
        counts(:, id) = histcounts(valid_spikes(valid_spikes(:, 1) == id, 2), edges).';
    end
    
    n_active = nnz(sum(counts, 1) > 0);
    rates_hz = counts / bin_s;
    
    % Apply Gaussian smoothing
    if smooth_sigma_s > 0
        sigma_bins = smooth_sigma_s / bin_s;
        radius = max(1, ceil(4 * sigma_bins));
        kernel = exp(-0.5 * ((-radius:radius) / sigma_bins).^2);
        kernel = kernel / sum(kernel);
        activity = conv2(rates_hz, kernel(:), 'same');
    else
        activity = rates_hz;
    end
    
    % PCA coordinates extraction
    act_centered = activity - mean(activity, 1);
    [U, S_pca, ~] = svd(act_centered, 'econ');
    k_pca = min([num_components, size(U, 2), size(activity, 2)]);
    coords = U(:, 1:k_pca) * S_pca(1:k_pca, 1:k_pca);
    n_points = size(coords, 1);
    
    % Distance Matrix (Vectorized)
    norms = sum(coords.^2, 2);
    dist2 = norms + norms.' - 2 * (coords * coords.');
    distances = sqrt(max(dist2, 0));
    distances(1:n_points+1:end) = 0; % Force diagonal to 0
    
    % Calculate Epsilon via simple array sorting
    d_valid = distances(isfinite(distances));
    d_sort = sort(d_valid(:));
    epsilon = d_sort(max(1, round(epsilon_percentile / 100 * numel(d_sort))));
    
    % Check future recurrence
    fs = 1 / bin_s;
    excl_points = round(exclusion_s * fs);
    is_recurrent = false(n_points, 1);
    
    for i = 1:n_points
        first_future = i + excl_points;
        if first_future <= n_points
            is_recurrent(i) = any(distances(i, first_future:end) < epsilon);
        end
    end
    
    % Moving window for recurrence percentages
    win_points = max(1, round(window_s * fs));
    step_points = max(1, round(step_s * fs));
    starts = (1:step_points:(n_points - win_points + 1)).';
    
    rec_rates = zeros(numel(starts), 1);
    win_times = zeros(numel(starts), 1);
    
    for w_idx = 1:numel(starts)
        idx_start = starts(w_idx);
        idx_end = idx_start + win_points - 1;
        rec_rates(w_idx) = mean(is_recurrent(idx_start:idx_end));
        win_times(w_idx) = time_vector(idx_start + floor((win_points - 1) / 2));
    end
    
    % Coalescence detection
    over_threshold = rec_rates >= coal_threshold;
    idx_coal = [];
    if numel(over_threshold) >= consecutive_windows
        sum_consec = conv(double(over_threshold), ones(consecutive_windows, 1), 'valid');
        idx_coal = find(sum_consec == consecutive_windows, 1, 'first');
    end
    
    if isempty(idx_coal)
        t_coal = NaN;
        t_coal_rel = NaN;
    else
        t_coal = win_times(idx_coal);
        t_coal_rel = t_coal - min_time;
    end
    
    % Store data into arrays for final summary
    coal_time_s(cond_idx) = t_coal;
    coal_time_rel_s(cond_idx) = t_coal_rel;
    epsilon_array(cond_idx) = epsilon;
    num_neurons_array(cond_idx) = n_neurons;
    num_active_array(cond_idx) = n_active;
    num_spikes_array(cond_idx) = size(valid_spikes, 1);
    num_windows_array(cond_idx) = numel(rec_rates);
    
    %% 6.1 Condition Table & Figure
    cond_table = table(repmat(cond_idx, numel(win_times), 1), ...
        repmat(string(condition_names{cond_idx}), numel(win_times), 1), ...
        win_times, rec_rates * 100, over_threshold, ...
        'VariableNames', {'Condition', 'Name', 'Time_s', 'Recurrence_perc', 'Over_threshold'});
    
    fig = figure('Color', 'white', 'Name', sprintf('Recurrence - Condition %d', cond_idx), 'Visible', 'off');
    ax = axes(fig); hold(ax, 'on');
    plot(ax, win_times, rec_rates * 100, 'Color', main_color, 'LineWidth', 2, 'DisplayName', 'Window Recurrence');
    yline(ax, 100 * coal_threshold, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Threshold: %.0f %%', 100 * coal_threshold));
    
    if isfinite(t_coal)
        xline(ax, t_coal, '-.', 'Color', [0 0.55 0], 'LineWidth', 2, 'DisplayName', sprintf('Coalescence: %.2f s', t_coal));
    end
    
    xlim(ax, [min_time max_time]); ylim(ax, [0 105]);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Recurrence (%)');
    box(ax, 'on'); set(ax, 'FontSize', 11, 'LineWidth', 1);
    hold(ax, 'off');
    
    %% 6.2 Export files
    base_name = fullfile(output_folder, regexprep(sprintf('Recurrence_%s_cond_%d', file_prefix, cond_idx), '[^a-zA-Z0-9_-]', '_'));
    
    if save_figures
        % USE AUXILIARY FUNCTION: export_figure
        export_figure(fig, base_name, fig_width_cm, fig_height_cm, export_formats, 600);
        fprintf('Figure saved: %s\n', base_name);
    end
    close(fig);
    
    if save_results
        writetable(cond_table, [base_name '.csv']);
    end
end

%% 7. SUMMARY TABLE
summary_table = table((1:n_conditions)', string(condition_names(:)), ...
    coal_time_s, coal_time_rel_s, isfinite(coal_time_s), ...
    epsilon_array, num_neurons_array, num_active_array, ...
    num_spikes_array, num_windows_array, ...
    'VariableNames', {'Condition', 'Name', 'Coal_Time_s', 'Coal_Time_Rel_s', ...
    'Coalescence_Detected', 'Epsilon', 'Total_Neurons', 'Active_Neurons', ...
    'Total_Spikes', 'Num_Windows'});

disp(' ');
disp('--- Recurrence Analysis Summary ---');
disp(summary_table);