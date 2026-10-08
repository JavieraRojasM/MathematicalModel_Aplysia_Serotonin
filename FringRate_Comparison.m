%% FIRING RATE COMPARISON (DATA VS MODELS)
%
% This script compares the firing rates of four experimental recordings against
% the firing rates obtained from a selectable set of optimized Izhikevich models.
%
% It first automatically filters the model directory to find permutations that
% successfully reached the optimization targets. Then, it calculates the per-neuron
% firing rates (Hz) for both the experimental data and the valid models.

clear; 
clc; 
close all;

%% 1. PATHS AND DIRECTORIES
project_folder = pwd;
addpath(fullfile(project_folder, 'Functions'), '-end');

%% 2. FILES AND MODELS TO ANALYZE
data_file = "C:\Users\javie\Desktop\Data\Tesis_Data\Sep0622_serotonin.mat";
models_folder = fullfile(pwd, 'Grid_Results', 'Models');
heterogeneity_seed = 9; 

%% 3. AUTOMATIC FILTERING (BY INDIVIDUAL CONDITION)
fprintf('\n====================================================\n');
fprintf('FILTERING SUCCESSFUL MODELS\n');
fprintf('====================================================\n');

search_pattern = sprintf('Modelo_IZH_seed_%03d_perm_*.mat', heterogeneity_seed);
found_files = dir(fullfile(models_folder, search_pattern));
found_files = found_files(~[found_files.isdir]);

if isempty(found_files)
    error('No files found in:\n%s', models_folder);
end

model_perm_numbers = [];
success_map = false(numel(found_files), 4); 

for k = 1:numel(found_files)
    mat_path = fullfile(found_files(k).folder, found_files(k).name);
    
    token = regexp(found_files(k).name, 'perm_(\d+)\.mat$', 'tokens', 'once');
    if isempty(token), continue; end
    current_perm = str2double(token{1});
    
    try
        S = load(mat_path, 'T_mejores');
        if isfield(S, 'T_mejores')
            current_success = S.T_mejores.Objetivo_alcanzado == true;
            
            if any(current_success) 
                model_perm_numbers(end+1) = current_perm;
                success_map(numel(model_perm_numbers), :) = current_success;
                fprintf(' [OK] Permutation %03d included (Success in G: %s).\n', ...
                    current_perm, mat2str(find(current_success)));
            else
                fprintf(' [X]  Permutation %03d discarded (failed all recordings).\n', current_perm);
            end
        else
            fprintf(' [?]  Permutation %03d missing valid T_mejores.\n', current_perm);
        end
    catch
        fprintf(' [!]  Error reading permutation %03d.\n', current_perm);
    end
end

if isempty(model_perm_numbers)
    error('No permutations achieved the target in any recording.');
end

% Sort mathematically
[model_perm_numbers, sort_idx] = sort(model_perm_numbers);
success_map = success_map(sort_idx, :);
model_file_pattern = sprintf('Modelo_IZH_seed_%03d_perm_%%03d.mat', heterogeneity_seed);
spk_vars = compose("spks_%d", 1:4);
condition_names = ["No serotonin", "Serotonin injection", "Maintained serotonin", "Washout, no serotonin"];

%% 4. TIME WINDOWS
data_start_s = 605;  data_end_s = 1005;
model_start_s = 15;  model_end_s = 415;

%% 5. OUTPUT CONFIGURATION & COLORS
output_folder = fullfile(project_folder,'Grid_Results', 'Opt_Results_Freq_Comparison');
if ~isfolder(output_folder), mkdir(output_folder); end

fig_width_cm = 15; fig_height_cm = 10;
export_formats = {'pdf', 'png', 'svg', 'fig'};

% Colors using auxiliary or inline definitions
color_blue = [0.20 0.60 0.80];
color_orange = [0.80 0.40 0.20];
point_colors = [color_blue; color_orange; color_orange; color_blue];
segment_colors = [color_blue; color_orange; color_blue];

%% 6. VALIDATIONS
if ~isfile(data_file), error('Experimental file not found:\n%s', data_file); end
if ~isfolder(models_folder), error('Models folder not found:\n%s', models_folder); end
n_conditions = numel(spk_vars);
n_models = numel(model_perm_numbers);

%% 7. EXPERIMENTAL DATA FREQUENCY
fprintf('\n===============================================\n');
fprintf('EXPERIMENTAL DATA\n');
fprintf('===============================================\n');

[data_median_hz, data_q1_hz, data_q3_hz, data_freq_per_neuron, data_n_neurons, data_status] = ...
    analyze_firing_rates(data_file, spk_vars, data_start_s, data_end_s);

for cond = 1:n_conditions
    fprintf('Condition %d: Median = %.4f Hz; Q1-Q3 = [%.4f, %.4f] Hz; N = %d.\n', ...
        cond, data_median_hz(cond), data_q1_hz(cond), data_q3_hz(cond), data_n_neurons(cond));
end

%% 8. MODEL MEDIANS (PER MUTATION)
fprintf('\n===============================================\n');
fprintf('MODELS (Permutations)\n');
fprintf('===============================================\n');

model_medians_hz = nan(n_models, n_conditions);
model_q1_hz = nan(n_models, n_conditions);
model_q3_hz = nan(n_models, n_conditions);
model_n_neurons = zeros(n_models, n_conditions);
model_freq_per_neuron = cell(n_models, n_conditions);
model_status = strings(n_models, n_conditions);
model_files_list = strings(n_models, 1);

for row = 1:n_models
    perm_num = model_perm_numbers(row);
    mat_name = sprintf(model_file_pattern, perm_num);
    mat_path = string(fullfile(models_folder, mat_name));
    model_files_list(row) = mat_path;
    
    if ~isfile(mat_path)
        model_status(row, :) = "File not found"; continue;
    end
    
    [model_medians_hz(row,:), model_q1_hz(row,:), model_q3_hz(row,:), ...
     model_freq_per_neuron(row,:), model_n_neurons(row,:), model_status(row,:)] = ...
        analyze_firing_rates(mat_path, spk_vars, model_start_s, model_end_s);
end

%% 9. MEDIAN OF MODEL MEDIANS
median_of_medians_hz = nan(1, n_conditions);
q1_across_models_hz = nan(1, n_conditions);
q3_across_models_hz = nan(1, n_conditions);
valid_models_count = zeros(1, n_conditions);

for cond = 1:n_conditions
    vals = model_medians_hz(:, cond);
    vals = vals(isfinite(vals));
    valid_models_count(cond) = numel(vals);
    
    if ~isempty(vals)
        median_of_medians_hz(cond) = median(vals);
        q1_across_models_hz(cond) = prctile(vals, 25);
        q3_across_models_hz(cond) = prctile(vals, 75);
    end
end

fprintf('\n===============================================\n');
fprintf('SUMMARY: MEDIAN OF MEDIANS [Q1, Q3]\n');
fprintf('===============================================\n');
for cond = 1:n_conditions
    fprintf('Condition %d: %.4f [%.4f, %.4f] Hz; Valid models = %d/%d.\n', ...
        cond, median_of_medians_hz(cond), q1_across_models_hz(cond), ...
        q3_across_models_hz(cond), valid_models_count(cond), n_models);
end

%% 10. TREND SIMILARITY & STATISTICS
fprintf('\n===============================================\n');
fprintf('SIMILARITY & STATISTICS\n');
fprintf('===============================================\n');

% 10.1 Cosine Similarity (Macroscopic Trend)
v_data = data_median_hz(:);
v_model = median_of_medians_hz(:);
valid_idx = isfinite(v_data) & isfinite(v_model);

cos_sim = NaN; rmse_trend = NaN;
if any(valid_idx)
    cos_sim = dot(v_data(valid_idx), v_model(valid_idx)) / (norm(v_data(valid_idx)) * norm(v_model(valid_idx)));
    rmse_trend = sqrt(mean((v_data(valid_idx) - v_model(valid_idx)).^2));
end
fprintf('Cosine Similarity (Trend): %.4f\n', cos_sim);
fprintf('RMSE (Magnitude difference): %.4f Hz\n\n', rmse_trend);

% 10.2 Wilcoxon Rank-Sum Test (Neuronal Distributions)
p_values_wilcoxon = nan(1, n_conditions);
fprintf('Wilcoxon Rank-Sum Test:\n');

for c = 1:n_conditions
    freq_exp = data_freq_per_neuron{c};
    freq_exp = freq_exp(isfinite(freq_exp));
    
    freq_mod_pool = [];
    for m = 1:n_models
        if isfinite(model_medians_hz(m, c))
            f_curr = model_freq_per_neuron{m, c};
            freq_mod_pool = [freq_mod_pool; f_curr(isfinite(f_curr))]; %#ok<AGROW>
        end
    end
    
    try
        if ~isempty(freq_exp) && ~isempty(freq_mod_pool)
            p_values_wilcoxon(c) = ranksum(freq_exp, freq_mod_pool);
        end
        if p_values_wilcoxon(c) > 0.05
            verdict = "No significant difference (Replicates biology)";
        else
            verdict = "Significantly different (Deterministic/scale bias)";
        end
        fprintf('  Condition %d: p-value = %.4f -> %s\n', c, p_values_wilcoxon(c), verdict);
    catch
        fprintf('  Condition %d: Test failed (Missing toolbox or empty data)\n', c);
    end
end

%% 11. TABLES
T_data = table((1:n_conditions)', condition_names(:), spk_vars(:), ...
    repmat(data_start_s, n_conditions, 1), repmat(data_end_s, n_conditions, 1), ...
    data_n_neurons(:), data_median_hz(:), data_q1_hz(:), data_q3_hz(:), data_status(:), ...
    'VariableNames', {'Condition', 'Name', 'Spk_Var', 'Start_s', 'End_s', ...
    'Num_Neurons', 'Median_Freq_Hz', 'Q1_Hz', 'Q3_Hz', 'Status'});

T_summary_models = table((1:n_conditions)', condition_names(:), ...
    median_of_medians_hz(:), q1_across_models_hz(:), q3_across_models_hz(:), ...
    valid_models_count(:), p_values_wilcoxon(:), ...
    'VariableNames', {'Condition', 'Name', 'Median_of_Medians_Hz', ...
    'Q1_Across_Models_Hz', 'Q3_Across_Models_Hz', 'Valid_Models', 'p_Value_Wilcoxon'});

%% 12. COMPARATIVE PLOT
fig = figure('Color', 'white', 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');

x_coords = 1:n_conditions;

% Draw IQR Bands for Model
for i = 1:(n_conditions-1)
    if all(isfinite([q1_across_models_hz(i), q1_across_models_hz(i+1), q3_across_models_hz(i), q3_across_models_hz(i+1)]))
        patch(ax, [x_coords(i) x_coords(i+1) x_coords(i+1) x_coords(i)], ...
            [q1_across_models_hz(i) q1_across_models_hz(i+1) q3_across_models_hz(i+1) q3_across_models_hz(i)], ...
            segment_colors(i,:), 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

% Draw Lines
for i = 1:(n_conditions-1)
    if all(isfinite([data_median_hz(i), data_median_hz(i+1)]))
        plot(ax, x_coords(i:i+1), data_median_hz(i:i+1), '-', 'Color', segment_colors(i,:), 'LineWidth', 2.2, 'HandleVisibility', 'off');
    end
    if all(isfinite([median_of_medians_hz(i), median_of_medians_hz(i+1)]))
        plot(ax, x_coords(i:i+1), median_of_medians_hz(i:i+1), '--', 'Color', segment_colors(i,:), 'LineWidth', 2.2, 'HandleVisibility', 'off');
    end
end

% Draw Scatters
for c = 1:n_conditions
    scatter(ax, x_coords(c), data_median_hz(c), 70, 'o', 'MarkerFaceColor', point_colors(c,:), ...
        'MarkerEdgeColor', point_colors(c,:), 'LineWidth', 1.0, 'HandleVisibility', 'off');
    scatter(ax, x_coords(c), median_of_medians_hz(c), 82, 's', 'MarkerFaceColor', point_colors(c,:), ...
        'MarkerEdgeColor', point_colors(c,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

xlim(ax, [0.65 n_conditions+0.35]);
xticks(ax, x_coords);
yticks(ax, 0:0.1:0.5);
xticklabels(ax, {'Condition 1', 'Condition 2', 'Condition 3', 'Condition 4'});
ylabel(ax, 'Mean Firing Rate per Neuron (Hz)');

all_vals = [data_median_hz(:); median_of_medians_hz(:); q3_across_models_hz(:)];
all_vals = all_vals(isfinite(all_vals));
if isempty(all_vals) || max(all_vals) <= 0
    ylim(ax, [0 1]);
else
    ylim(ax, [0 0.5]);
end

box(ax, 'on');
set(ax, 'FontSize', 11, 'LineWidth', 1, 'Layer', 'top');
hold(ax, 'off');

%% 13. SAVE RESULTS USING export_figure
results_file = fullfile(output_folder, 'Freq_Comparison_Results.mat');
save(results_file, 'T_data', 'T_summary_models', 'data_freq_per_neuron', 'model_freq_per_neuron', ...
    'data_median_hz', 'data_q1_hz', 'data_q3_hz', 'model_medians_hz', 'model_q1_hz', 'model_q3_hz', ...
    'median_of_medians_hz', 'q1_across_models_hz', 'q3_across_models_hz', 'cos_sim', 'rmse_trend', ...
    'p_values_wilcoxon', 'data_n_neurons', 'model_n_neurons', 'valid_models_count', 'model_perm_numbers', ...
    'model_files_list', 'condition_names', 'spk_vars', 'data_start_s', 'data_end_s', 'model_start_s', 'model_end_s');

base_name = fullfile(output_folder, 'Comparison_Median_Freq');
export_figure(fig, base_name, fig_width_cm, fig_height_cm, export_formats, 300);
close(fig);

fprintf('\nResults saved in:\n%s\n', output_folder);

%% LOCAL HELPER FUNCTION (Using format_spikes)
function [medians, q1, q3, freqs_per_neuron, num_neurons, statuses] = analyze_firing_rates(file_path, spk_vars, start_s, end_s)
    n_conds = numel(spk_vars);
    medians = nan(1, n_conds); q1 = nan(1, n_conds); q3 = nan(1, n_conds);
    freqs_per_neuron = cell(1, n_conds);
    num_neurons = zeros(1, n_conds);
    statuses = strings(1, n_conds);
    duration_s = end_s - start_s;
    
    info = whos('-file', char(file_path));
    avail_vars = string({info.name});
    vars_to_load = spk_vars(ismember(spk_vars, avail_vars));
    
    if ismember("x", avail_vars), vars_to_load(end+1) = "x"; end
    if ismember("p", avail_vars), vars_to_load(end+1) = "p"; end
    
    if isempty(vars_to_load)
        statuses(:) = "No spks_N variables found"; return;
    end
    
    vars_cell = cellstr(vars_to_load);
    S = load(char(file_path), vars_cell{:});
    
    % Determine total neurons
    if isfield(S, 'x') && ~isempty(S.x)
        n_total = numel(S.x);
    elseif isfield(S, 'p') && isfield(S.p, 'N')
        n_total = double(S.p.N);
    else
        n_total = 0;
        for c = 1:n_conds
            var_name = char(spk_vars(c));
            if isfield(S, var_name) && ~isempty(S.(var_name))
                [~, curr_n] = format_spikes(S.(var_name), [], [], 0);
                n_total = max(n_total, curr_n);
            end
        end
    end
    if n_total == 0, n_total = 1; end
    
    for c = 1:n_conds
        var_name = char(spk_vars(c));
        if ~isfield(S, var_name)
            statuses(c) = "Variable missing"; continue;
        end
        
        try
            % USO DE LA FUNCION AUXILIAR: format_spikes
            [valid_spks, ~] = format_spikes(S.(var_name), start_s, end_s, n_total);
            
            if isempty(valid_spks)
                freqs = zeros(n_total, 1);
            else
                freqs = accumarray(valid_spks(:, 1), 1, [n_total 1], @sum, 0) / duration_s;
            end
            
            freqs_per_neuron{c} = freqs;
            num_neurons(c) = numel(freqs);
            
            valid_freqs = freqs(isfinite(freqs));
            if isempty(valid_freqs)
                statuses(c) = "No valid frequencies"; continue;
            end
            
            medians(c) = median(valid_freqs);
            q1(c) = prctile(valid_freqs, 25);
            q3(c) = prctile(valid_freqs, 75);
            statuses(c) = "Success";
        catch ME
            statuses(c) = "Error: " + string(ME.message);
        end
    end
end