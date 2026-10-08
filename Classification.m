%% CLASSIFICATION ANALYSIS (DATA AND MODEL) - UNIFIED
%
% Neuronal classification based on firing patterns.
%
% Data mode : Reduces the experimental recordings, runs the Bruno Toolbox,
%             computes percentages, and saves the reference Class_data_Bruno.mat.
% Model mode: Aggregates all successful permutations of a model seed, computes
%             Median and IQR, and compares them with the experimental
%             reference (L1, Cosine, Chi-square). Saves a consolidated plot.

clear; 
clc; 
close all;

%% 1. PATHS AND DIRECTORIES
project_folder = pwd;
addpath(fullfile(project_folder, 'Functions'), '-end');
test_data_folder = "C:\Users\javie\Desktop\Data\TestData"; % Required by Bruno
if isfolder(test_data_folder), addpath(test_data_folder, '-end'); end

%% 2. SELECT DATA OR MODEL
% ================================================================
% LEAVE ONLY ONE OF THE FOLLOWING TWO LINES ACTIVE:
% mode_type = 'data';
mode_type = 'model';
% ================================================================
% Experimental Data File
data_path = "C:\Users\javie\Desktop\Data\Tesis_Data\Sep0622_serotonin.mat";

% Model File
model_path = fullfile(project_folder, 'Grid_Results', 'Models', 'Modelo_IZH_seed_009_perm_002.mat');
% ================================================================

%% 3. COMMON CONFIGURATION
n_conditions = 4;
n_total_neurons = 49; % Required for Chi-square expected counts
fig_width_cm = 20;
fig_height_cm = 15;

% Figure export settings controlled by auxiliary function
save_png = true;
save_svg = true;
save_fig = true;
export_formats = {'pdf'};
if save_png, export_formats{end+1} = 'png'; end
if save_svg, export_formats{end+1} = 'svg'; end
if save_fig, export_formats{end+1} = 'fig'; end

types = ["Non_oscillatory"; "Oscillator"; "Burster"; "Pauser"];
n_types = numel(types);
cond_labels = ["Condition 1", "Condition 2", "Condition 3", "Condition 4"];

%% 4. EXECUTE SELECTED MODE
switch lower(mode_type)
    case 'data'
        %% ====================================================================
        %                          DATA MODE
        % ====================================================================
        fprintf('\n====================================================\n');
        fprintf('RUNNING CLASSIFICATION ON EXPERIMENTAL DATA\n');
        fprintf('====================================================\n');
        
        min_time = 605;
        max_time = 1005;
        
        bruno_out_folder = fullfile(project_folder, 'Bruno_Analysis');
        if ~isfolder(bruno_out_folder), mkdir(bruno_out_folder); end
        addpath(bruno_out_folder, '-end');
        
        if exist('ReduceDataFx', 'file') ~= 2 || exist('RunBrunoAnalysis', 'file') ~= 2
            error('Bruno Toolbox functions not found. Check TestData and Functions paths.');
        end
        
        % --- 4.1 REDUCE RECORDINGS ---
        fprintf('Reducing data for %d conditions...\n', n_conditions);
        DataList = struct('spikes', cell(1, n_conditions));
        for g = 1:n_conditions
            fprintf('  Processing condition %d...\n', g);
            DataList(g).spikes = ReduceDataFx(char(data_path), min_time, max_time, g, char(test_data_folder));
        end
        
        % --- 4.2 RUN BRUNO TOOLBOX ---
        fprintf('\nRunning Bruno Toolbox (this may take a moment)...\n');
        bruno_filename = 'Class_data_Bruno.mat';
        
        old_dir = pwd; clean_dir = onCleanup(@() cd(old_dir));
        cd(bruno_out_folder);
        Classification_Bruno = RunBrunoAnalysis(DataList, char(test_data_folder), char(bruno_out_folder), 1, bruno_filename);
        cd(old_dir); clear clean_dir;
        
        if ~isstruct(Classification_Bruno) || ~isfield(Classification_Bruno, 'Summary')
            error('RunBrunoAnalysis did not return a valid Summary.');
        end
        
        % --- 4.3 CALCULATE PERCENTAGES (INLINED) ---
        Summary_counts = Classification_Bruno.Summary;
        pct_table = Summary_counts; % Copy to retain structure (Neuron_type)
        
        if istable(Summary_counts)
            is_num_col = varfun(@isnumeric, Summary_counts, 'OutputFormat', 'uniform');
            num_cols_idx = find(is_num_col);
            
            counts_mat = double(table2array(Summary_counts(:, num_cols_idx)));
            totals = sum(counts_mat, 1, 'omitnan');
            pct_mat = 100 * (counts_mat ./ totals);
            pct_mat(:, totals == 0) = NaN; % Handle empty recordings
            
            for k = 1:length(num_cols_idx)
                col_name = Summary_counts.Properties.VariableNames{num_cols_idx(k)};
                pct_table.(col_name) = pct_mat(:, k);
            end
        else
            counts_mat = double(Summary_counts);
            totals = sum(counts_mat, 1, 'omitnan');
            pct_mat = 100 * (counts_mat ./ totals);
            pct_mat(:, totals == 0) = NaN;
            pct_table = pct_mat;
        end
        
        Classification_Bruno.Summary_percentages = pct_table;
        
        % --- 4.4 PLOT & SAVE (NATIVE GENERATION) ---
        results_file = fullfile(bruno_out_folder, bruno_filename);
        save(results_file, 'Classification_Bruno', 'Summary_counts', 'pct_table', 'DataList', 'data_path', '-v7.3');
        fprintf('\nData classification complete. Saved reference to:\n%s\n', results_file);
        
        cond_names_data = {'Recording 1', 'Recording 2', 'Recording 3', 'Recording 4'};
        x_axis = 1:n_conditions;

        % Plot Total Counts
        fig_counts = figure('Color', 'w', 'Renderer', 'painters', 'Visible', 'off');
        ax_counts = axes(fig_counts); hold(ax_counts, 'on');
        
        counts_data = table2array(Summary_counts(:, 2:end));
        for t = 1:n_types
            c = get_project_colors(types(t));
            plot(ax_counts, x_axis, counts_data(t, :), '-o', 'LineWidth', 2.0, ...
                'MarkerSize', 8, 'MarkerFaceColor', c, 'Color', c, 'DisplayName', types(t));
        end
        xticks(ax_counts, x_axis); xticklabels(ax_counts, cond_names_data);
        ylabel(ax_counts, 'Number of Neurons'); xlabel(ax_counts, 'Experimental Condition');
        legend(ax_counts, 'Location', 'northeast'); box(ax_counts, 'on');
        ylim(ax_counts, [0, max(counts_data(:)) + 7]); xlim(ax_counts, [0.8, 4.2]);
        hold(ax_counts, 'off');
        
        export_figure(fig_counts, fullfile(bruno_out_folder, 'Class_Counts'), fig_width_cm, fig_height_cm, export_formats, 300);
        close(fig_counts);
        
        % Plot Percentages
        fig_pct = figure('Color', 'w', 'Renderer', 'painters', 'Visible', 'off');
        ax_pct = axes(fig_pct); hold(ax_pct, 'on');
        
        pct_data = table2array(pct_table(:, 2:end));
        for t = 1:n_types
            c = get_project_colors(types(t));
            plot(ax_pct, x_axis, pct_data(t, :), '-o', 'LineWidth', 2.0, ...
                'MarkerSize', 8, 'MarkerFaceColor', c, 'Color', c, 'DisplayName', types(t));
        end
        xticks(ax_pct, x_axis); xticklabels(ax_pct, cond_names_data);
        ylabel(ax_pct, 'Classified Neurons (%)'); xlabel(ax_pct, 'Experimental Condition');
        legend(ax_pct, 'Location', 'northeast'); box(ax_pct, 'on');
        ylim(ax_pct, [0, 105]); xlim(ax_pct, [0.8, 4.2]);
        hold(ax_pct, 'off');
        
        export_figure(fig_pct, fullfile(bruno_out_folder, 'Class_Percentages'), fig_width_cm, fig_height_cm, export_formats, 300);
        close(fig_pct);
        
    case 'model'
        %% ====================================================================
        %                          MODEL MODE
        % ====================================================================
        [model_dir, model_name, ~] = fileparts(char(model_path));
        output_folder = fullfile(fileparts(model_dir), 'Classification_Permutations_Results');
        if ~isfolder(output_folder), mkdir(output_folder); end
        
        token = regexp(model_name, 'Modelo_IZH_seed_(\d+)_perm_', 'tokens', 'once');
        if isempty(token), error('Could not determine base seed from model_path.'); end
        seed_het = str2double(token{1});
        
        fprintf('\n====================================================\n');
        fprintf('EXTRACTING MODEL CLASSIFICATION PERCENTAGES\n');
        fprintf('Catalog Seed: %d\n', seed_het);
        fprintf('====================================================\n');
        
        files = dir(fullfile(model_dir, sprintf('Modelo_IZH_seed_%03d_perm_*.mat', seed_het)));
        files = files(~[files.isdir]);
        if isempty(files), error('No permutation files found in %s', model_dir); end
        
        % --- 4.1 AGGREGATE PERMUTATIONS ---
        all_pct = nan(n_types, numel(files) * n_conditions);
        cond_ids = nan(1, numel(files) * n_conditions);
        count_valid = 0;
        perms_with_success = 0;
        
        for i = 1:numel(files)
            mat_file = fullfile(files(i).folder, files(i).name);
            S = load(mat_file, 'T_mejores', 'porcentajes_modelo');
            
            if ~isfield(S, 'T_mejores') || ~isfield(S, 'porcentajes_modelo'), continue; end
            
            success = S.T_mejores.Objetivo_alcanzado == true;
            if ~any(success), continue; end
            perms_with_success = perms_with_success + 1;
            
            for g = 1:n_conditions
                if success(g)
                    count_valid = count_valid + 1;
                    all_pct(:, count_valid) = S.porcentajes_modelo(:, g);
                    cond_ids(count_valid) = g;
                end
            end
        end
        
        if count_valid == 0, error('No successful permutations found.'); end
        all_pct = all_pct(:, 1:count_valid);
        cond_ids = cond_ids(1:count_valid);
        
        % Calculate Median and IQR across permutations
        med_pct_model = nan(n_types, n_conditions);
        q1_pct_model = nan(n_types, n_conditions);
        q3_pct_model = nan(n_types, n_conditions);
        
        for g = 1:n_conditions
            cols = cond_ids == g;
            if any(cols)
                cond_data = all_pct(:, cols);
                med_pct_model(:, g) = median(cond_data, 2, 'omitnan');
                q1_pct_model(:, g)  = prctile(cond_data, 25, 2);
                q3_pct_model(:, g)  = prctile(cond_data, 75, 2);
            end
        end
        
        % --- 4.2 EXPERIMENTAL COMPARISON & STATISTICS (INLINED) ---
        ref_file = fullfile(project_folder, 'Bruno_Analysis', 'Class_data_Bruno.mat');
        exp_pct = []; exp_types = types;
        
        if isfile(ref_file)
            try
                ref_data = load(ref_file);
                if isfield(ref_data, 'pct_table')
                    exp_table = ref_data.pct_table;
                    if istable(exp_table)
                        first_col = exp_table{:,1};
                        if iscellstr(first_col) || isstring(first_col)
                            exp_types = string(first_col);
                            exp_table(:,1) = [];
                        end
                        exp_pct = double(table2array(exp_table));
                    else
                        exp_pct = exp_table;
                    end
                elseif isfield(ref_data, 'Classification_Bruno')
                    exp_table = ref_data.Classification_Bruno.Summary;
                    if istable(exp_table)
                        first_col = exp_table{:,1};
                        if iscellstr(first_col) || isstring(first_col), exp_table(:,1) = []; end
                        counts = double(table2array(exp_table));
                    else
                        counts = double(exp_table);
                    end
                    exp_pct = 100 * (counts ./ sum(counts, 1, 'omitnan'));
                end
                
                % Align Experimental Data Rows to standard 'types'
                exp_aligned = nan(n_types, n_conditions);
                for t = 1:n_types
                    norm_target = regexprep(lower(char(types(t))), '[^a-z]', '');
                    for e = 1:length(exp_types)
                        norm_exp = regexprep(lower(char(exp_types(e))), '[^a-z]', '');
                        if strcmp(norm_target, norm_exp)
                            exp_aligned(t, 1:n_conditions) = exp_pct(e, 1:n_conditions);
                            break;
                        end
                    end
                end
                exp_pct = exp_aligned;
                
                fprintf('\n--- CLASSIFICATION SIMILARITY (BY CONDITION) ---\n');
                for g = 1:n_conditions
                    p_real = exp_pct(:, g) / 100;
                    p_sim  = med_pct_model(:, g) / 100;
                    
                    if all(isfinite(p_real)) && all(isfinite(p_sim))
                        L1_err = sum(abs(p_sim - p_real)) * 100;
                        cosine_sim = dot(p_sim, p_real) / (norm(p_sim) * norm(p_real));
                        
                        obs = (p_sim * n_total_neurons) + 1;
                        exp = (p_real * n_total_neurons) + 1;
                        chi2_stat = sum(((obs - exp).^2) ./ exp);
                        try p_val = 1 - chi2cdf(chi2_stat, length(p_real) - 1); catch, p_val = NaN; end
                        
                        verdict = "Significant diff. (Model DEVIATES)";
                        if p_val > 0.05, verdict = "No significant diff. (Model REPLICATES)"; end
                        
                        fprintf('>> Condition %d: L1 Error: %.2f pp | Cosine: %.4f | Chi2 p-val: %.4f -> %s\n', ...
                                g, L1_err, cosine_sim, p_val, verdict);
                    end
                end
            catch ME
                warning('Could not calculate similarity: %s', ME.message);
            end
        else
            fprintf('\n[!] Experimental reference not found. Skipping statistics.\n');
        end
        
        % --- 4.3 INLINED GLOBAL PLOT (MEDIAN + IQR) ---
        fig_global = figure('Color', 'w', 'Visible', 'off', 'Name', 'Model vs Experimental Classification');
        ax_global = axes('Parent', fig_global);
        hold(ax_global, 'on');
        x_coords = 1:n_conditions;
        max_y = 0;
        
        for t = 1:n_types
            c = get_project_colors(types(t));
            
            % Model Data
            m_med = med_pct_model(t, :);
            m_q1  = max(0, q1_pct_model(t, :)); 
            m_q3  = q3_pct_model(t, :);
            
            % Draw IQR Patch
            patch(ax_global, [x_coords, fliplr(x_coords)], [m_q1, fliplr(m_q3)], ...
                c, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                
            % Draw Model Line
            plot(ax_global, x_coords, m_med, '--s', 'Color', c, 'LineWidth', 2.2, ...
                'MarkerSize', 7, 'MarkerFaceColor', c, 'DisplayName', sprintf('Model - %s', types(t)));
                
            % Draw Experimental Data if available
            if ~isempty(exp_pct) && any(isfinite(exp_pct(t, :)))
                exp_val = exp_pct(t, :);
                plot(ax_global, x_coords, exp_val(1:n_conditions), '-o', 'Color', c, ...
                    'LineWidth', 2.2, 'MarkerSize', 8, 'MarkerFaceColor', c, ...
                    'DisplayName', sprintf('Data - %s', types(t)));
                max_y = max([max_y, max(m_q3), max(exp_val)]);
            else
                max_y = max([max_y, max(m_q3)]);
            end
        end
        
        % Formatting Axes
        xticks(ax_global, x_coords); xticklabels(ax_global, cond_labels);
        xlim(ax_global, [0.75, n_conditions + 0.25]);
        if max_y <= 0 || isnan(max_y), max_y = 100; end
        ylim(ax_global, [0, min(100, max_y + 10)]); 
        
        ylabel(ax_global, 'Classified Neurons (%)');
        title(ax_global, sprintf('Percentage Classification | Median & IQR (%d Perms, Seed %d)', perms_with_success, seed_het));
        legend(ax_global, 'Location', 'eastoutside', 'Interpreter', 'none');
        ax_global.FontSize = 10;
        box(ax_global, 'on'); hold(ax_global, 'off');
        
        % Save Plot using export_figure auxiliary function
        base_name = fullfile(output_folder, 'Classification_Percentages_Median_IQR');
        export_figure(fig_global, base_name, fig_width_cm, fig_height_cm, export_formats, 600);
        close(fig_global);
        
        fprintf('\nModel extraction complete. Global plot saved in:\n%s\n', output_folder);
end