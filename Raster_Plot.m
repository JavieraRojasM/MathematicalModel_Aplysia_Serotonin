%% RATE-SORTED RASTER PLOT (DATA AND MODEL)
% Generates raster plots sorted by firing frequency for the four
% experimental recordings or the four conditions of an Izhikevich model.

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
% Experimental file
data_path = "C:\Users\javie\Desktop\Data\Tesis_Data\Sep0622_serotonin.mat";

% Model file
model_path =  fullfile(project_folder, 'Grid_Results', 'Models', 'Modelo_IZH_seed_009_perm_002.mat');
% ================================================================

%% 3. COMMON CONFIGURATION
language = 'EN';
dt = 1e-3;                 % Time step [s]
n_conditions = 4;

% Sorting behavior:
% 'ascend': lowest frequency at the bottom, highest at the top.
% 'descend': highest frequency at the bottom, lowest at the top.
sort_direction = 'descend'; 

spike_height = 0.8;
spike_thickness = 0.8;
fig_width_cm = 15;
fig_height_cm = 10;

% Figure export settings controlled by auxiliary function
save_png = true;
save_svg = true;
export_formats = {'pdf'};
if save_png, export_formats{end+1} = 'png'; end
if save_svg, export_formats{end+1} = 'svg'; end

close_figures = true;

%% 4. SOURCE CONFIGURATION
switch lower(mode_type)
    case 'data'
        input_path = data_path;
        min_plot_time = 605;
        max_plot_time = 1005;
        condition_names = {'No serotonin'; 'Serotonin injection'; ...
                           'Maintained serotonin'; 'Washout, No serotonin'};
        file_prefix = 'data';
        output_folder = fullfile(project_folder, 'Rate_Sorted_Rasters', 'Data');
    case 'model'
        input_path = model_path;
        min_plot_time = 15;
        max_plot_time = 415;
        [~, model_name] = fileparts(model_path);
        model_name = regexprep(model_name, '[^a-zA-Z0-9_-]', '_');
        file_prefix = model_name;
        output_folder = fullfile(project_folder, 'Rate_Sorted_Rasters', 'Model', model_name);
    otherwise
        error("mode_type must be 'data' or 'model'.");
end

if ~isfile(input_path), error('File not found: %s', input_path); end
if ~isfolder(output_folder), mkdir(output_folder); end

%% 5. LOAD AND VALIDATE DATA
data = load(input_path);
for k = 1:n_conditions
    if ~isfield(data, sprintf('spks_%d', k))
        error('The file does not contain the variable spks_%d.', k);
    end
end

% Retrieve condition parameters (model) or setup defaults (data)
combinations = nan(n_conditions, 3);
if strcmpi(mode_type, 'model')
    if isfield(data, 'modelo') && isfield(data.modelo, 'combinaciones')
        combinations = double(data.modelo.combinaciones);
        condition_names = cellstr(compose('C%d: a=%.4g, b=%.4g, d=%.4g', (1:n_conditions)', combinations));
    else
        condition_names = cellstr(compose('Recording %d', (1:n_conditions)'));
    end
end

% Extract stimulus time (converted to seconds)
stim_time_s = NaN;
if isfield(data, 'stim_time') && ~isempty(data.stim_time)
    if isnumeric(data.stim_time)
        stim_time_s = double(data.stim_time(1)) * 60;
    else
        val = sscanf(char(data.stim_time), '%f', 1);
        if ~isempty(val), stim_time_s = val * 60; end
    end
end

analysis_duration = max_plot_time - min_plot_time;
if analysis_duration <= 0
    error('max_plot_time must be greater than min_plot_time.');
end

% Setup sorting multiplier
sign_order = 1;
if strcmpi(sort_direction, 'descend'), sign_order = -1; end

% Determine expected number of neurons from coordinates (if available)
expected_neurons = 0;
if isfield(data, 'x') && ~isempty(data.x), expected_neurons = numel(data.x); end

%% 6. GENERATE RASTERS
freq_sort_data = repmat(struct('recording', [], 'name', '', ...
    'sorted_IDs', [], 'sorted_frequencies_Hz', [], 'frequencies_by_ID_Hz', []), n_conditions, 1);

for cond_idx = 1:n_conditions
    fprintf('Processing condition %d...\n', cond_idx);
    
    [PS, spks_raw] = Plot_Settings(cond_idx, data, language);
    
    % USE AUXILIARY FUNCTION: format_spikes
    [spks_filtered, n_neurons] = format_spikes(spks_raw, min_plot_time, max_plot_time, expected_neurons);
    if expected_neurons == 0, expected_neurons = n_neurons; end
    
    original_IDs = (1:n_neurons)';
    
    % Count spikes and calculate frequencies directly
    if isempty(spks_filtered)
        spike_counts = zeros(n_neurons, 1);
    else
        spike_counts = accumarray(spks_filtered(:,1), 1, [n_neurons, 1], @sum, 0);
    end
    freq_Hz = spike_counts / analysis_duration;
    
    % Sort neurons: primarily by firing rate, secondarily by ID (to handle ties)
    [~, sort_idx] = sortrows([sign_order * freq_Hz, original_IDs], [1, 2]);
    sorted_IDs = original_IDs(sort_idx);
    sorted_freqs = freq_Hz(sort_idx);
    
    % Map original IDs to new Y-axis positions based on sorting
    new_y_positions = zeros(n_neurons, 1);
    new_y_positions(sorted_IDs) = 1:n_neurons;
    
    % Store data for this condition
    freq_sort_data(cond_idx).recording = cond_idx;
    freq_sort_data(cond_idx).name = condition_names{cond_idx};
    freq_sort_data(cond_idx).sorted_IDs = sorted_IDs;
    freq_sort_data(cond_idx).sorted_frequencies_Hz = sorted_freqs;
    freq_sort_data(cond_idx).frequencies_by_ID_Hz = freq_Hz;
    
    %% --- PLOTTING ---
    fig_name = sprintf('Rate-Sorted Raster - %s %d', mode_type, cond_idx);
    fig = figure('Color', 'w', 'Name', fig_name, 'Visible', 'off');
    ax = axes(fig); hold(ax, 'on');
    
    % Fast raster plotting using vectorized NaNs separation 
    if ~isempty(spks_filtered)
        spike_times = spks_filtered(:, 2)';
        y_centers = new_y_positions(spks_filtered(:, 1))';
        
        x_lines = [spike_times; spike_times; nan(1, length(spike_times))];
        y_lines = [y_centers - spike_height/2; y_centers + spike_height/2; nan(1, length(spike_times))];
        
        plot(ax, x_lines(:), y_lines(:), 'Color', PS.color_raster, ...
            'LineWidth', spike_thickness, 'HandleVisibility', 'off');
    end
    
    % Axis formatting
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Neurons');
    xlim(ax, [min_plot_time max_plot_time]);
    ylim(ax, [0.5 n_neurons+0.5]);
    yticks(ax, []); % Hide Y-ticks as they don't represent raw IDs anymore
    
    % Draw Stimulus Line if applicable
    if isfinite(stim_time_s) && stim_time_s >= min_plot_time && stim_time_s <= max_plot_time
        xline(ax, stim_time_s, 'r-', 'LineWidth', 2, 'HandleVisibility', 'off');
        text(ax, stim_time_s, n_neurons + 0.3, 'Stimulus', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'Color', 'red');
    end
    
    set(ax, 'Layer', 'top', 'TickDir', 'out', 'Box', 'on');
    hold(ax, 'off');
    
    %% --- SAVING EXPORTS ---
    base_name = fullfile(output_folder, sprintf('Rate_Sorted_Raster_%s_cond_%d', file_prefix, cond_idx));
    
    % USE AUXILIARY FUNCTION: export_figure
    export_figure(fig, base_name, fig_width_cm, fig_height_cm, export_formats, 600);
    
    if close_figures, close(fig); else, set(fig, 'Visible', 'on'); end
end

%% 7. SAVE SORTING DATA TO MAT FILE
mat_file_path = fullfile(output_folder, 'Firing_Rate_Order.mat');
save(mat_file_path, 'freq_sort_data', 'mode_type', 'input_path', ...
    'min_plot_time', 'max_plot_time', 'dt', 'sort_direction', ...
    'combinations', 'expected_neurons');

fprintf('\nProcess completed. Results saved in:\n%s\n', output_folder);