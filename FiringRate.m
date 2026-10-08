%% FIRING RATE ANALYSIS (DATA AND MODEL)
% Analyzes the firing rate of four experimental conditions or four
% simulations contained in a Model_IZH_seed_XXX.mat file.

%% FIRING RATE ANALYSIS (DATA AND MODEL)
% Analyzes the firing rate of four experimental conditions or four
% simulations contained in a Model_IZH_seed_XXX.mat file.
%
% Generates:
%   1. A temporal population firing rate curve per condition.
%   2. The mean firing frequency for each neuron.
%   3. A summary of median and 25-75th percentiles.
%   4. Pairwise comparisons using Wilcoxon and Holm correction.
%   5. CSV tables and a MAT file with all results.

clear; 
clc; 
close all;

%% 1. PATHS AND DIRECTORIES
project_folder = pwd;
addpath(fullfile(project_folder, 'Functions'), '-end');

%% 2. SELECT DATA OR MODEL
% ================================================================
% LEAVE ONLY ONE OF THE FOLLOWING TWO LINES ACTIVE:
mode_type = 'data';
% mode_type = 'model';
% ================================================================
% Experimental file
data_path = "C:\Users\javie\Desktop\Data\Tesis_Data\Sep0622_serotonin.mat";

% Model file
model_path =  fullfile(project_folder, 'Grid_Results', 'Models', 'Modelo_IZH_seed_009_perm_002.mat');
% ================================================================

%% 3. COMMON CONFIGURATION
language = 'EN';
bin_s = 0.5;               % Resolution for temporal rate curve [s]
smooth_window_s = 10;      % Moving average window [s]
n_conditions = 4;

fig_width_temporal_cm = 15;
fig_height_temporal_cm = 10;

fig_width_summary_cm = 15;
fig_height_summary_cm = 10;

% Figure export settings controlled by auxiliary function
save_png = true;
save_svg = true;
export_formats = {'pdf'};
if save_png, export_formats{end+1} = 'png'; end
if save_svg, export_formats{end+1} = 'svg'; end

%% 4. SOURCE CONFIGURATION
switch lower(mode_type)
    case 'data'
        input_path = data_path;
        min_time = 605;
        max_time = 1005;
        y_limits = [0 1.5];
        condition_names = {'Recording 1'; 'Recording 2'; 'Recording 3'; 'Recording 4'};
        file_prefix = 'data';
        results_folder = fullfile(project_folder, 'Firing_Rate_Results', 'Data');
    case 'model'
        input_path = model_path;
        min_time = 15;
        max_time = 415;
        y_limits = [0 1]; 
        [~, model_name] = fileparts(model_path);
        model_name = regexprep(model_name, '[^a-zA-Z0-9_-]', '_');
        file_prefix = model_name;
        results_folder = fullfile(project_folder, 'Firing_Rate_Results', 'Model', model_name);
    otherwise
        error("mode_type must be 'data' or 'model'.");
end

if ~isfile(input_path), error('File not found: %s', input_path); end
if ~isfolder(results_folder), mkdir(results_folder); end

data = load(input_path);
for k = 1:n_conditions
    if ~isfield(data, sprintf('spks_%d', k))
        error('File is missing variable spks_%d.', k);
    end
end

% Retrieve condition combinations if model
combinations = nan(n_conditions, 3);
if strcmpi(mode_type, 'model')
    if isfield(data, 'modelo') && isfield(data.modelo, 'combinaciones')
        combinations = double(data.modelo.combinaciones);
        condition_names = cell(n_conditions, 1);
        for k = 1:n_conditions
            condition_names{k} = sprintf('C%d: a=%.4g, b=%.4g, d=%.4g', ...
                k, combinations(k,1), combinations(k,2), combinations(k,3));
        end
    else
        condition_names = cellstr(compose('Condition %d', (1:n_conditions)'));
    end
end

% Determine the total number of neurons across all conditions
num_neurons = 0;
if isfield(data, 'x') && isnumeric(data.x), num_neurons = numel(data.x); end

%% 5. INITIALIZE RESULT VARIABLES
freq_per_neuron = cell(n_conditions, 1);
median_freq_hz = nan(n_conditions, 1);
prctile_25_hz = nan(n_conditions, 1);
prctile_75_hz = nan(n_conditions, 1);
valid_neurons_count = zeros(n_conditions, 1);
cond_colors = zeros(n_conditions, 3);
analysis_duration_s = max_time - min_time;

%% 6. PROCESS EACH CONDITION
for cond = 1:n_conditions
    [PS, spks_raw] = Plot_Settings(cond, data, language);
    
    main_color = validatecolor(PS.color_raster);
    light_color = 0.5 * main_color + 0.5 * [1 1 1];
    cond_colors(cond, :) = main_color;
    
    % USE AUXILIARY FUNCTION: format_spikes
    [valid_spikes, current_n_neurons] = format_spikes(spks_raw, min_time, max_time, num_neurons);
    if num_neurons == 0, num_neurons = current_n_neurons; end
    
    %% 6.1 Mean Firing Rate per Neuron
    if isempty(valid_spikes)
        spike_counts = zeros(num_neurons, 1);
    else
        spike_counts = accumarray(valid_spikes(:,1), 1, [num_neurons, 1], @sum, 0);
    end
    
    neuron_freqs = spike_counts / analysis_duration_s;
    neuron_freqs(~isfinite(neuron_freqs)) = NaN;
    freq_per_neuron{cond} = neuron_freqs;
    
    valid_freqs = neuron_freqs(isfinite(neuron_freqs));
    valid_neurons_count(cond) = numel(valid_freqs);
    
    if ~isempty(valid_freqs)
        median_freq_hz(cond) = median(valid_freqs);
        prctile_25_hz(cond) = prctile(valid_freqs, 25);
        prctile_75_hz(cond) = prctile(valid_freqs, 75);
    end
    
    %% 6.2 Population Rate in Bins
    bin_edges = min_time:bin_s:max_time;
    if length(bin_edges) < 2
        warning('Condition %d does not have enough bins.', cond);
        continue;
    end
    
    time_vector = bin_edges(1:end-1) + bin_s / 2;
    spikes_per_bin = histcounts(valid_spikes(:,2), bin_edges);
    
    % Convert to Hz per neuron
    pop_rate_real = spikes_per_bin / (num_neurons * bin_s);
    
    % Moving average smoothing
    smooth_pts = max(1, round(smooth_window_s / bin_s));
    pop_rate_smooth = movmean(pop_rate_real, smooth_pts, 'omitnan');
    
    %% 6.3 Temporal Figure
    fig = figure('Color', 'w', 'Visible', 'off');
    ax = axes(fig); hold(ax, 'on');
    
    plot(ax, time_vector, pop_rate_real, 'Color', light_color, 'LineWidth', 1, ...
         'DisplayName', sprintf('Measured every %.1f s', bin_s));
    plot(ax, time_vector, pop_rate_smooth, 'Color', main_color, 'LineWidth', 2, ...
         'DisplayName', sprintf('%d s moving average', smooth_window_s));
         
    yline(ax, median_freq_hz(cond), '--', 'Color', main_color, 'LineWidth', 1.5);
    
    xlim(ax, [min_time max_time]);
    if ~isempty(y_limits), ylim(ax, y_limits); end
    box(ax, 'on'); hold(ax, 'off');
    
    %% 6.4 Save Temporal Figure using export_figure
    base_name = fullfile(results_folder, sprintf('Freq_%s_cond_%d', file_prefix, cond));
    export_figure(fig, base_name, fig_width_temporal_cm, fig_height_temporal_cm, export_formats, 600);
    close(fig);
end

%% 7. PAIRWISE COMPARISONS
comparison_pairs = [1 2; 2 3; 3 4; 1 4];
n_comps = size(comparison_pairs, 1);
p_original = nan(n_comps, 1);
paired_neurons_count = zeros(n_comps, 1);

for comp = 1:n_comps
    c1 = comparison_pairs(comp, 1);
    c2 = comparison_pairs(comp, 2);
    
    v1 = freq_per_neuron{c1};
    v2 = freq_per_neuron{c2};
    
    if isempty(v1) || isempty(v2), continue; end
    if numel(v1) ~= numel(v2)
        warning('Conditions %d and %d have different neuron counts.', c1, c2);
        continue;
    end
    
    valid_pairs = isfinite(v1) & isfinite(v2);
    paired_neurons_count(comp) = sum(valid_pairs);
    
    if paired_neurons_count(comp) >= 3
        x1 = v1(valid_pairs);
        x2 = v2(valid_pairs);
        if all(x1 == x2), p_original(comp) = 1;
        else, p_original(comp) = signrank(x1, x2); end
    end
end

%% 8. HOLM CORRECTION AND SIGNIFICANCE (Using apply_holm_correction)
[p_adjusted, significance] = apply_holm_correction(p_original);

%% 9. TABLES
freq_table = table((1:n_conditions)', string(condition_names(:)), ...
    median_freq_hz, prctile_25_hz, prctile_75_hz, valid_neurons_count, ...
    'VariableNames', {'Condition', 'Name', 'Median_Freq_Hz', ...
                      'Percentile_25_Hz', 'Percentile_75_Hz', 'Num_Neurons'});

name_c1 = string(condition_names(comparison_pairs(:, 1)));
name_c2 = string(condition_names(comparison_pairs(:, 2)));

comparisons_table = table(comparison_pairs(:,1), comparison_pairs(:,2), ...
    name_c1, name_c2, paired_neurons_count, p_original, p_adjusted, significance, ...
    'VariableNames', {'Cond_1', 'Cond_2', 'Name_1', 'Name_2', ...
                      'Paired_Neurons', 'p_Original', 'p_Holm_Adj', 'Significance'});

fprintf('\n--- Firing Rate Summary ---\n'); disp(freq_table);
fprintf('\n--- Pairwise Comparisons ---\n'); disp(comparisons_table);

%% 10. SUMMARY FIGURE: MEDIAN AND IQR
fig_summary = figure('Color', 'w', 'Visible', 'off');
ax_res = axes(fig_summary); hold(ax_res, 'on');

bars = bar(ax_res, 1:n_conditions, median_freq_hz, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1);
bars.CData = cond_colors;

error_low = median_freq_hz - prctile_25_hz;
error_high = prctile_75_hz - median_freq_hz;
errorbar(ax_res, 1:n_conditions, median_freq_hz, error_low, error_high, ...
         'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

set(ax_res, 'XTick', 1:n_conditions, 'XTickLabel', condition_names, 'FontSize', 11, 'TickLabelInterpreter', 'none');
ylabel(ax_res, 'Median Firing Rate per Neuron (Hz)');
xlim(ax_res, [0.5 n_conditions+0.5]);
box(ax_res, 'on');

%% 11. SIGNIFICANCE BRACKETS
base_height = max(prctile_75_hz, [], 'omitnan');
if isempty(base_height) || ~isfinite(base_height) || base_height <= 0, base_height = 1; end
step_h = 0.12 * base_height;

for comp = 1:n_comps
    c1 = comparison_pairs(comp, 1);
    c2 = comparison_pairs(comp, 2);
    h = base_height + comp * step_h;
    
    plot(ax_res, [c1 c1 c2 c2], [h-step_h/4 h h h-step_h/4], 'k-', 'LineWidth', 1.2);
    text(ax_res, mean([c1 c2]), h + step_h/10, significance(comp), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', 13, 'FontWeight', 'bold');
end

ylim(ax_res, [0, base_height + (n_comps + 1.2) * step_h]);
hold(ax_res, 'off');

%% 12. SAVE OVERALL RESULTS USING export_figure
summary_base_name = fullfile(results_folder, 'Summary_Median_Freq');
export_figure(fig_summary, summary_base_name, fig_width_summary_cm, fig_height_summary_cm, export_formats, 600);
close(fig_summary);

save(fullfile(results_folder, 'Complete_Firing_Rate_Results.mat'), ...
    'mode_type', 'input_path', 'min_time', 'max_time', 'freq_per_neuron', ...
    'median_freq_hz', 'prctile_25_hz', 'prctile_75_hz', 'freq_table', ...
    'comparisons_table', 'combinations', '-v7.3');
    
fprintf('\nResults saved in:\n%s\n', results_folder);