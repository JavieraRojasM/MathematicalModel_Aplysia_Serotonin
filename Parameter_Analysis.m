%% PARAMETER TRAJECTORIES AND DELTAS ANALYSIS (PARTIAL FILTERING)
%
% This script extracts the optimal parameters of the models, keeping only the 
% permutations that reached the optimization target in at least one recording. 
% It computes paired changes (deltas) between consecutive conditions and 
% generates an independent trajectory plot for each parameter (J_E, J_I, a_E, b_E, d_E).
% The plots include individual permutation paths, the median trajectory, and the 
% interquartile range (IQR) shaded band, using distinct colors for each condition 
% transition (Blue -> Orange -> Blue).

clc; 
clear; 
close all;

%% 1. CONFIGURATION AND PATHS
project_folder = pwd;
addpath(fullfile(project_folder, 'Functions'), '-end');

models_folder = fullfile(pwd, 'Grid_Results', 'Models');
out_folder = fullfile('Grid_Results', 'Parameter_Analysis');
if ~isfolder(out_folder), mkdir(out_folder); end

% Parameters to analyze
params = ["J_E", "J_I", "a_E", "b_E", "d_E"];
param_labels = {'J_E', 'J_I', 'a_E', 'b_E', 'd_E'};

% Comparisons: [Source_Recording, Target_Recording]
comparisons = [1 2; 2 3; 3 4; 1 4];
trans_names = ["G1_to_G2", "G2_to_G3", "G3_to_G4", "G1_to_G4"];

cond_names = {'No serotonin', 'Serotonin application', 'Maintained serotonin', 'Washout'};
cond_labels = {'Condition 1', 'Condition 2', 'Condition 3', 'Condition 4'};

%% 2. FIGURE CONFIGURATION
fig_w_cm = 15; fig_h_cm = 10;

% Figure export settings controlled by auxiliary function
save_png = true;
save_svg = true;
save_fig = true;
export_formats = {'pdf'};
if save_png, export_formats{end+1} = 'png'; end
if save_svg, export_formats{end+1} = 'svg'; end
if save_fig, export_formats{end+1} = 'fig'; end

% Colors
c_blue = [0.20 0.60 0.80];
c_orange = [0.80 0.40 0.20];
c_points = [c_blue; c_orange; c_orange; c_blue];
c_segments = [c_blue; c_orange; c_blue];

%% 3. SEARCH FOR MODEL FILES
files = dir(fullfile(models_folder, 'Modelo_IZH_seed_*_perm_*.mat'));
if isempty(files), error('No models found in:\n%s', models_folder); end

fprintf('\n===============================================\n');
fprintf('FILES FOUND\n');
fprintf('===============================================\n');
fprintf('Found %d files in:\n%s\n', numel(files), models_folder);

%% 4. EXTRACT AND FILTER SUCCESSFUL PERMUTATIONS
T_absolute = table();
perm_numbers = [];
success_map = false(0, 4);

req_cols = ["Grabacion", "Objetivo_alcanzado", params];

for i = 1:numel(files)
    file_path = fullfile(files(i).folder, files(i).name);
    try
        S = load(file_path, 'T_mejores', 'modelo');
        if ~isfield(S, 'T_mejores') || ~istable(S.T_mejores)
            fprintf(' [?]  File %s missing valid T_mejores.\n', files(i).name); continue;
        end
        
        T_temp = S.T_mejores;
        missing_cols = setdiff(req_cols, string(T_temp.Properties.VariableNames));
        if ~isempty(missing_cols)
            fprintf(' [?]  File %s missing columns: %s.\n', files(i).name, strjoin(missing_cols, ', ')); continue;
        end
        
        % Extract permutation seed
        if isfield(S, 'modelo') && isfield(S.modelo, 'semilla_permutacion')
            curr_perm = double(S.modelo.semilla_permutacion);
        else
            token = regexp(files(i).name, 'perm_(\d+)\.mat$', 'tokens', 'once');
            if isempty(token), continue; end
            curr_perm = str2double(token{1});
        end
        
        % Check success per condition (1 to 4)
        curr_success = false(1, 4);
        for g = 1:4
            rows_g = T_temp.Grabacion == g;
            curr_success(g) = any(T_temp.Objetivo_alcanzado(rows_g) == true);
        end
        
        if any(curr_success)
            perm_numbers(end+1, 1) = curr_perm; 
            success_map(end+1, :) = curr_success;
            T_temp.Permutation_Seed = repmat(curr_perm, height(T_temp), 1);
            T_absolute = [T_absolute; T_temp]; %#ok<AGROW>
            fprintf(' [OK] Permutation %03d included (Success in G: %s).\n', curr_perm, mat2str(find(curr_success)));
        else
            fprintf(' [X]  Permutation %03d discarded (failed all recordings).\n', curr_perm);
        end
    catch ME
        fprintf(' [!]  Error reading permutation from file %s: %s\n', files(i).name, ME.message);
    end
end

if isempty(perm_numbers), error('No permutation achieved the target in any recording.'); end

% Sort permutations and maintain alignment
[perm_numbers, sort_idx] = sort(perm_numbers);
success_map = success_map(sort_idx, :);
T_absolute = sortrows(T_absolute, {'Permutation_Seed', 'Grabacion'});
unique_perms = perm_numbers;

T_success_map = table(perm_numbers, success_map(:,1), success_map(:,2), success_map(:,3), success_map(:,4), ...
    'VariableNames', {'Permutation_Seed', 'G1', 'G2', 'G3', 'G4'});

fprintf('\n===============================================\n');
fprintf('FILTERING SUMMARY\n');
fprintf('===============================================\n');
fprintf('Permutations with at least one success: %d\n', numel(unique_perms));

%% 5. CALCULATE PAIRED DELTAS (Delta = Target - Source)
delta_col_names = ["Permutation_Seed", "Source_Rec", "Target_Rec", "Transition", params];
delta_col_types = ["double", "double", "double", "string", repmat("double", 1, numel(params))];

T_deltas = table('Size', [0 numel(delta_col_names)], 'VariableTypes', cellstr(delta_col_types), 'VariableNames', cellstr(delta_col_names));

for i = 1:numel(unique_perms)
    perm = unique_perms(i);
    perm_data = T_absolute(T_absolute.Permutation_Seed == perm, :);
    
    for c = 1:size(comparisons, 1)
        src_rec = comparisons(c, 1);
        tgt_rec = comparisons(c, 2);
        
        src_row = perm_data(perm_data.Grabacion == src_rec & perm_data.Objetivo_alcanzado == true, :);
        tgt_row = perm_data(perm_data.Grabacion == tgt_rec & perm_data.Objetivo_alcanzado == true, :);
        
        if isempty(src_row) || isempty(tgt_row), continue; end
        
        src_row = src_row(1, :); tgt_row = tgt_row(1, :);
        
        delta_row = table();
        delta_row.Permutation_Seed = perm;
        delta_row.Source_Rec = src_rec;
        delta_row.Target_Rec = tgt_rec;
        delta_row.Transition = trans_names(c);
        
        for p = 1:numel(params)
            p_name = char(params(p));
            delta_row.(p_name) = tgt_row.(p_name) - src_row.(p_name);
        end
        T_deltas = [T_deltas; delta_row]; %#ok<AGROW>
    end
end

%% 6. DELTA SUMMARY TABLE
T_summary_deltas = table();
for p = 1:numel(params)
    p_name = char(params(p));
    for c = 1:numel(trans_names)
        vals = T_deltas.(p_name)(T_deltas.Transition == trans_names(c));
        vals = vals(isfinite(vals));
        
        res_row = table();
        res_row.Parameter = string(p_name);
        res_row.Transition = trans_names(c);
        res_row.N = numel(vals);
        
        if isempty(vals)
            res_row.Median_Delta = NaN; res_row.Q1_Delta = NaN; res_row.Q3_Delta = NaN;
        else
            res_row.Median_Delta = median(vals);
            res_row.Q1_Delta = prctile(vals, 25);
            res_row.Q3_Delta = prctile(vals, 75);
        end
        T_summary_deltas = [T_summary_deltas; res_row]; %#ok<AGROW>
    end
end

%% 8. PLOTTING PARAMETER TRAJECTORIES
x_coords = 1:4;

for p = 1:numel(params)
    p_name = char(params(p));
    val_matrix = nan(numel(unique_perms), 4);
    
    for i = 1:numel(unique_perms)
        perm = unique_perms(i);
        perm_data = T_absolute(T_absolute.Permutation_Seed == perm, :);
        for g = 1:4
            row = perm_data(perm_data.Grabacion == g & perm_data.Objetivo_alcanzado == true, :);
            if ~isempty(row), val_matrix(i, g) = row.(p_name)(1); end
        end
    end
    
    % Median and IQR
    param_median = nan(1, 4); param_q1 = nan(1, 4); param_q3 = nan(1, 4);
    for g = 1:4
        vals_g = val_matrix(:, g);
        vals_g = vals_g(isfinite(vals_g));
        if ~isempty(vals_g)
            param_median(g) = median(vals_g);
            param_q1(g) = prctile(vals_g, 25);
            param_q3(g) = prctile(vals_g, 75);
        end
    end
    
    %% Setup Figure
    fig = figure('Name', sprintf('Trajectory of %s', p_name), 'Color', 'white', 'Visible', 'off');
    ax = axes(fig); hold(ax, 'on');
    
    %% Draw IQR Band (Inlined)
    for seg = 1:3
        if all(isfinite([param_q1(seg) param_q1(seg+1) param_q3(seg) param_q3(seg+1)]))
            patch(ax, [x_coords(seg) x_coords(seg+1) x_coords(seg+1) x_coords(seg)], ...
                [param_q1(seg) param_q1(seg+1) param_q3(seg+1) param_q3(seg)], ...
                c_segments(seg,:), 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    
    %% Draw Median Segments (Inlined)
    for seg = 1:3
        if all(isfinite([param_median(seg), param_median(seg+1)]))
            plot(ax, x_coords(seg:seg+1), param_median(seg:seg+1), ':', ...
                'Color', c_segments(seg,:), 'LineWidth', 2.8, 'HandleVisibility', 'off');
        end
    end
    
    %% Draw Scatter Points
    for g = 1:4
        if isfinite(param_median(g))
            scatter(ax, x_coords(g), param_median(g), 82, 's', 'MarkerFaceColor', c_points(g,:), ...
                'MarkerEdgeColor', c_points(g,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    %% Axis Formatting
    lims = val_matrix(isfinite(val_matrix));
    if isempty(lims), ylim(ax, [0 1]);
    else
        min_y = min(lims); max_y = max(lims);
        range_y = max_y - min_y;
        margin_y = max(0.12 * range_y, 0.10);
        ylim(ax, [min_y - margin_y, max_y + margin_y]);
    end
    
    xlim(ax, [0.65 4.35]);
    xticks(ax, x_coords);
    xticklabels(ax, cond_labels);
    ylabel(ax, param_labels{p}, 'Interpreter', 'tex');
    box(ax, 'on');
    set(ax, 'FontSize', 11, 'LineWidth', 1, 'Layer', 'top', 'TickLabelInterpreter', 'none');
    hold(ax, 'off');
    
    %% Save Figure using export_figure auxiliary function
    base_name = fullfile(out_folder, sprintf('Trajectory_param_%s', p_name));
    export_figure(fig, base_name, fig_w_cm, fig_h_cm, export_formats, 600);
    close(fig);
    
    fprintf('Figure saved for parameter %s.\n', p_name);
end

%% 9. SAVE RESULTS TO MAT
results_file = fullfile(out_folder, 'Parameter_Trajectories_Results.mat');
save(results_file, 'T_absolute', 'T_deltas', 'T_summary_deltas', 'T_success_map', ...
    'unique_perms', 'perm_numbers', 'success_map', 'params', 'param_labels', ...
    'comparisons', 'trans_names', 'cond_names', 'cond_labels', 'models_folder', ...
    'out_folder', 'c_points', 'c_segments', '-v7.3');
    
fprintf('\n===============================================\n');
fprintf('ANALYSIS COMPLETED\n');
fprintf('===============================================\n');
fprintf('Results saved in:\n%s\n', out_folder);