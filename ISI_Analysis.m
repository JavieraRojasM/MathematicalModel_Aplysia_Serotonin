%% ISI ANALYSIS (DATA AND MODEL)
%
% Analyzes the Inter-Spike Intervals (ISI) of four conditions from
% experimental data or an Izhikevich model.
%
% Generates:
%   1. An independent histogram per condition.
%   2. A summary of the neuronal median and 25-75th percentiles.
%   3. Pairwise comparisons using Wilcoxon (signrank) and Holm correction.
%   4. CSV tables and a MAT file with all results.

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
data_path = "C:\Users\javie\Desktop\Data\Tesis_Data\Sep0622_serotonin.mat";
model_path =  fullfile(project_folder, 'Grid_Results', 'Models', 'Modelo_IZH_seed_009_perm_002.mat');

%% 3. COMMON CONFIGURATION
language = 'EN';
n_conditions = 4;
num_bins = 40;
display_percentile = 99;

fig_width_cm = 15;
fig_height_cm = 10;

% Figure export settings controlled by auxiliary function
save_png = true;
save_svg = true;
export_formats = {'pdf'};
if save_png, export_formats{end+1} = 'png'; end
if save_svg, export_formats{end+1} = 'svg'; end

% Pairwise comparisons to perform (Row = one pair of conditions)
comparison_pairs = [1 2; 2 3; 3 4; 1 4];

%% 4. SOURCE CONFIGURATION
switch lower(mode_type)
    case 'data'
        input_path = data_path;
        min_time = 605; max_time = 1005;
        condition_names = {'Recording 1'; 'Recording 2'; 'Recording 3'; 'Recording 4'};
        source_name = 'Experimental Data'; file_prefix = 'data';
        y_limits_hist = [0 0.3];
        combinations = nan(n_conditions, 3);
        results_folder = fullfile(project_folder, 'ISI_Results', 'Data');
    case 'model'
        input_path = model_path;
        min_time = 30; max_time = 270;
        condition_names = cellstr(compose('Simulation %d', (1:n_conditions)'));
        [~, model_name] = fileparts(model_path);
        model_name = regexprep(model_name, '[^a-zA-Z0-9_-]', '_');
        source_name = strrep(model_name, '_', ' '); file_prefix = model_name;
        y_limits_hist = []; % Auto limits to avoid cutting model histograms
        results_folder = fullfile(project_folder, 'ISI_Results', 'Model', model_name);
end

if ~isfile(input_path), error('File not found: %s', input_path); end
if ~isfolder(results_folder), mkdir(results_folder); end

data = load(input_path);
for k = 1:n_conditions
    if ~isfield(data, sprintf('spks_%d', k)), error('File is missing spks_%d.', k); end
end

% Infer number of neurons
num_neurons = 0;
if isfield(data, 'x') && isnumeric(data.x), num_neurons = numel(data.x); end

if strcmpi(mode_type, 'model') && isfield(data, 'modelo') && isfield(data.modelo, 'combinaciones')
    combinations = double(data.modelo.combinaciones);
end

fprintf('\nISI Analysis: %s\n', source_name);
fprintf('Analyzed Window: %.0f-%.0f s\n\n', min_time, max_time);

%% 5. CALCULATE ISI (Optimized)
all_isi_by_cond = cell(n_conditions, 1);
median_isi_per_neuron = cell(n_conditions, 1);
cond_colors = zeros(n_conditions, 3);
median_of_medians_s = nan(n_conditions, 1);

for cond = 1:n_conditions
    [PS, spks_raw] = Plot_Settings(cond, data, language);
    cond_colors(cond, :) = validatecolor(PS.color_raster);
    
    % USE AUXILIARY FUNCTION: format_spikes
    [valid_spikes, current_n_neurons] = format_spikes(spks_raw, min_time, max_time, num_neurons);
    if num_neurons == 0, num_neurons = current_n_neurons; end
    
    neuron_medians = nan(num_neurons, 1);
    current_cond_isi = [];
    
    if isempty(valid_spikes)
        warning('Condition %d contains no spikes in window.', cond);
    else
        % Fast calculation of ISI per neuron
        for neuron_id = 1:num_neurons
            neuron_times = sort(valid_spikes(valid_spikes(:,1) == neuron_id, 2));
            if length(neuron_times) >= 2
                isi_vals = diff(neuron_times);
                isi_vals = isi_vals(isi_vals > 0); % Strictly positive ISIs
                if ~isempty(isi_vals)
                    current_cond_isi = [current_cond_isi; isi_vals]; %#ok<AGROW>
                    neuron_medians(neuron_id) = median(isi_vals);
                end
            end
        end
    end
    
    all_isi_by_cond{cond} = current_cond_isi;
    median_isi_per_neuron{cond} = neuron_medians;
    
    valid_medians = neuron_medians(isfinite(neuron_medians));
    if ~isempty(valid_medians)
        median_of_medians_s(cond) = median(valid_medians);
    end
end

%% 6. CREATE COMPARABLE HISTOGRAMS
merged_isi = vertcat(all_isi_by_cond{:});
if isempty(merged_isi), error('No valid ISIs found in any condition.'); end

limit_isi = prctile(merged_isi, display_percentile);
if ~isfinite(limit_isi) || limit_isi <= 0, error('Invalid limit for histograms.'); end
isi_edges = linspace(0, limit_isi, num_bins + 1);

for cond = 1:n_conditions
    fig_hist = figure('Color', 'w', 'Visible', 'off'); ax = axes(fig_hist); hold(ax, 'on');
    
    if isempty(all_isi_by_cond{cond})
        text(ax, 0.5, 0.5, 'No valid ISIs found', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 12);
    else
        histogram(ax, all_isi_by_cond{cond}, 'BinEdges', isi_edges, 'Normalization', 'probability', ...
            'FaceColor', cond_colors(cond,:), 'FaceAlpha', 0.65, 'EdgeColor', cond_colors(cond,:), 'LineWidth', 1.1);
        xline(ax, median_of_medians_s(cond), '--', 'Color', 'k', 'LineWidth', 1, ...
            'Label', sprintf('Neuronal Median = %.2f s', median_of_medians_s(cond)), ...
            'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right');
    end
    
    xlim(ax, [-0.5 limit_isi]); yticks(ax, 0:0.1:0.3); xticks(ax, 0:2:10);
    if ~isempty(y_limits_hist), ylim(ax, y_limits_hist); end
    box(ax, 'on'); hold(ax, 'off');
    
    % USE AUXILIARY FUNCTION: export_figure
    base_name = fullfile(results_folder, sprintf('ISI_Hist_%s_cond_%d', file_prefix, cond));
    export_figure(fig_hist, base_name, fig_width_cm, fig_height_cm, export_formats, 600);
    close(fig_hist);
end

%% 7. CONDITION SUMMARY
prctile_25_s = nan(n_conditions, 1);
prctile_75_s = nan(n_conditions, 1);
valid_neurons_count = zeros(n_conditions, 1);

for cond = 1:n_conditions
    vals = median_isi_per_neuron{cond};
    vals = vals(isfinite(vals));
    valid_neurons_count(cond) = numel(vals);
    if ~isempty(vals)
        median_of_medians_s(cond) = median(vals);
        prctile_25_s(cond) = prctile(vals, 25);
        prctile_75_s(cond) = prctile(vals, 75);
    end
end

summary_table = table((1:n_conditions)', string(condition_names(:)), ...
    combinations(:,1), combinations(:,2), combinations(:,3), ...
    median_of_medians_s, prctile_25_s, prctile_75_s, valid_neurons_count, ...
    'VariableNames', {'Condition', 'Condition_Name', 'a_E', 'b_E', 'd_E', ...
                      'Neuronal_Median_ISI_s', 'Percentile_25_s', 'Percentile_75_s', 'Num_Neurons'});

disp('--- ISI Summary Table ---'); disp(summary_table);

%% 8. PAIRWISE COMPARISONS (Wilcoxon + Holm-Bonferroni)
num_comps = size(comparison_pairs, 1);
p_original = nan(num_comps, 1);
paired_neurons = zeros(num_comps, 1);

for comp = 1:num_comps
    c1 = comparison_pairs(comp, 1); c2 = comparison_pairs(comp, 2);
    v1 = median_isi_per_neuron{c1}; v2 = median_isi_per_neuron{c2};
    
    valid_pair = isfinite(v1) & isfinite(v2);
    paired_neurons(comp) = sum(valid_pair);
    
    if paired_neurons(comp) >= 3
        p_original(comp) = signrank(v1(valid_pair), v2(valid_pair));
    end
end

% USE AUXILIARY FUNCTION: apply_holm_correction
[p_adjusted, significance] = apply_holm_correction(p_original);

comparisons_table = table(comparison_pairs(:,1), comparison_pairs(:,2), ...
    string(condition_names(comparison_pairs(:, 1))), string(condition_names(comparison_pairs(:, 2))), ...
    paired_neurons, p_original, p_adjusted, significance, ...
    'VariableNames', {'Cond_1', 'Cond_2', 'Name_1', 'Name_2', 'Paired_Neurons', 'p_Original', 'p_Holm_Adj', 'Significance'});

fprintf('\n--- Pairwise Comparisons ---\n'); disp(comparisons_table);

%% 9. MEDIAN, ERROR, AND SIGNIFICANCE PLOT
fig_summary = figure('Color', 'w', 'Visible', 'off'); ax = axes(fig_summary); hold(ax, 'on');

bars = bar(ax, 1:n_conditions, median_of_medians_s, 'FaceColor', 'flat', 'EdgeColor', 'k');
bars.CData = cond_colors;

errorbar(ax, 1:n_conditions, median_of_medians_s, median_of_medians_s - prctile_25_s, prctile_75_s - median_of_medians_s, ...
         'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

set(ax, 'XTick', 1:n_conditions, 'XTickLabel', condition_names);
ylabel(ax, 'Neuronal Median ISI (s)'); box(ax, 'on');

base_height = max([1, max(prctile_75_s, [], 'omitnan')]);
step_h = 0.12 * base_height;

for comp = 1:num_comps
    c1 = comparison_pairs(comp, 1); c2 = comparison_pairs(comp, 2);
    h = base_height + comp * step_h;
    plot(ax, [c1 c1 c2 c2], [h-step_h/4 h h h-step_h/4], 'k-', 'LineWidth', 1.2);
    text(ax, mean([c1 c2]), h + step_h/10, significance(comp), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 13, 'FontWeight', 'bold');
end
ylim(ax, [0 base_height + (num_comps + 1) * step_h]); hold(ax, 'off');

% USE AUXILIARY FUNCTION: export_figure
base_name_summary = fullfile(results_folder, 'ISI_Median_Summary');
export_figure(fig_summary, base_name_summary, fig_width_cm, fig_height_cm, export_formats, 600);
close(fig_summary);

%% 10. SAVE RESULTS
writetable(summary_table, fullfile(results_folder, 'ISI_Summary_Table.csv'));
writetable(comparisons_table, fullfile(results_folder, 'ISI_Comparisons_Table.csv'));

save(fullfile(results_folder, 'ISI_Complete_Results.mat'), 'mode_type', 'input_path', 'source_name', ...
    'condition_names', 'combinations', 'num_neurons', 'min_time', 'max_time', 'all_isi_by_cond', ...
    'median_isi_per_neuron', 'median_of_medians_s', 'prctile_25_s', 'prctile_75_s', ...
    'summary_table', 'comparisons_table', '-v7.3');
    
fprintf('\nAnalysis finished. Results saved in:\n%s\n', results_folder);