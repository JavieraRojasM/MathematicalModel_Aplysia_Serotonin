%% GRID OPTIMIZATION
%
% This script performs a multi-objective grid search to optimize Izhikevich 
% network parameters (J_E, J_I, a_E, b_E, d_E) across four experimental 
% conditions. In Phase 1, it evaluates parameter combinations in parallel, 
% filtering models first by population firing rate, and then classifying 
% bursting behavior to compute an L1 error score (while automatically saving 
% checkpoints). In Phase 2, it selects the best-performing parameter set 
% for each condition, re-runs the full simulation to capture all detailed 
% network data (spikes, connectivity matrices, and stats), and exports the 
% final compiled model into a .mat file.
%

clc; 
clear; 
close all;

%% 1. USER CONFIGURATION
project_folder = fileparts(mfilename('fullpath'));
if isempty(project_folder), project_folder = pwd; end

% Add required paths
addpath(project_folder, '-begin');
if isfolder(fullfile(project_folder, 'Functions')), addpath(fullfile(project_folder, 'Functions'), '-end'); end
if isfolder(fullfile(project_folder, 'Functions2')), addpath(fullfile(project_folder, 'Functions2'), '-end'); end

heterogeneity_seeds = 9;
permutation_seeds = 1:62;
data_folder = fullfile(project_folder);
positions_file = ""; % Leave empty to auto-detect from experimental data

% Initialize output folders
output_folder = fullfile(project_folder, 'Grid_Results');
models_folder = fullfile(output_folder, 'Models');
checkpoints_folder = fullfile(output_folder, 'Checkpoints');
if ~isfolder(output_folder), mkdir(output_folder); end
if ~isfolder(models_folder), mkdir(models_folder); end
if ~isfolder(checkpoints_folder), mkdir(checkpoints_folder); end

final_files_list = strings(0, 1);
catalog_seed_record = [];
perm_seed_record = [];

for seed_2 = 1:numel(permutation_seeds)
    current_perm_seed = permutation_seeds(seed_2);
    fprintf('\n====================================================\n');
    fprintf('CATALOG SEED %d (PERMUTATION %d)\n', heterogeneity_seeds, current_perm_seed);
    fprintf('====================================================\n');
    
    resume_checkpoint = true;
    
    %% 2. GRID SPACE AND LIMITS DEFINITION
    cfg = struct();
    cfg.names = ["J_E", "J_I", "a_E", "b_E", "d_E"];
    
    % --- PHASE 1 (GRID) CONFIGURATION ---
    cfg.grid.J_E = 25:1:35;   
    cfg.grid.J_I = 10:1:18;   
    cfg.grid.a_E = 0.01:0.01:0.040; 
    cfg.grid.b_E = 0.20:0.01:0.250; 
    cfg.grid.d_E = 0.50:0.5:3;   
    
    cfg.search.parfor_batch_size = 15; 
    cfg.search.max_grid_evaluations = 500; % Exploration limit
    
    %% 3. SEQUENTIAL FILTER TOLERANCES (FREQUENCY ONLY)
    cfg.filters.rel_freq_tolerance = 0.35;
    cfg.filters.abs_freq_tolerance_hz = 0.35;
    cfg.filters.class_L1_tolerance_pp = 30;
    
    %% 4. WINDOWS AND FIXED MODEL PARAMETERS
    cfg.model.T_prefilters_ms = 410000;
    cfg.model.T_ms = 500000;
    cfg.model.metrics_window_s = [15 415];
    cfg.model.max_total_spikes = 30000;
    cfg.model.I_base_E = 3.0;
    cfg.model.sigma_noise = 0.55;
    
    cfg.model.seed_connectivity = 91;
    cfg.model.seed_dynamics = 12;
    cfg.model.seed_permutation = current_perm_seed;
    
    cfg.classification.duration_s = 400;
    cfg.classification.which_spks = 5;
    cfg.classification.seed_bruno = 700000 + (1:4);
    cfg.classification.types = ["Non_oscillatory"; "Oscillator"; "Burster"; "Pauser"];
        
%% 5. LOAD EXPERIMENTAL TARGETS AND POSITIONS
    % Load Firing Rate Targets
    freq_file_path = fullfile(data_folder, 'Firing_Rate_Results', 'Data', 'Complete_Firing_Rate_Results.mat');
    freq_data = load(freq_file_path);
    
    % Load Classification Targets
    class_file_path = fullfile(data_folder, 'Bruno_Analysis', 'Class_data_Bruno.mat');
    class_data = load(class_file_path);
    
    targets.freq_hz = double(freq_data.median_freq_hz(:).'); 

    numeric_cols = varfun(@isnumeric, class_data.pct_table, 'OutputFormat', 'uniform');
    targets.percentages = double(table2array(class_data.pct_table(:, numeric_cols)));

    if strlength(positions_file) == 0 
        if isfield(freq_data, 'input_path')
            positions_file = string(freq_data.input_path);
        else
            error('positions_file is empty and input_path was not found in freq_data.');
        end
    end
    
    % Load Positions
    pos_data = load(positions_file, 'x', 'xstd', 'y', 'ystd');
    x = double(pos_data.x(:)); xstd = double(pos_data.xstd(:)); 
    y = double(pos_data.y(:)); ystd = double(pos_data.ystd(:));
    pos_original = [x, y];
    pos_model = pos_original / max(max(pos_original, [], 1)); % Normalize
    
    heterogeneity_seeds = unique(double(heterogeneity_seeds(:).'), 'stable');
        
    %% --- PHASE 1: MULTI-OBJECTIVE GRID SEARCH ---
    [J, I, a, b, d] = ndgrid(cfg.grid.J_E, cfg.grid.J_I, cfg.grid.a_E, cfg.grid.b_E, cfg.grid.d_E);
    full_grid = [J(:), I(:), a(:), b(:), d(:)];
    n_combinations = size(full_grid, 1);
    
    % Initialize empty history table format
    empty_history = table('Size', [0 23], ...
        'VariableTypes', [repmat({'double'},1,9), repmat({'logical'},1,2), repmat({'double'},1,8), ...
                          repmat({'logical'},1,2), {'string','string'}], ...
        'VariableNames', {'J_E','J_I','a_E','b_E','d_E','Freq_Hz','Error_Freq_Norm','Error_Activity', ...
                          'Failed_Filters','Passes_Filters','Classified','Start_Class_s','End_Class_s', ...
                          'Error_L1_pp','Pct_NonOsc','Pct_Oscillator','Pct_Burster','Pct_Pauser', ...
                          'Score','Target_Reached','Valid_Simulation','Status','Message'});
    
    for iSeed = 1:numel(heterogeneity_seeds)
        seed = heterogeneity_seeds(iSeed);
        seed_folder = fullfile(checkpoints_folder, sprintf('Catalog_%03d_perm_%03d', seed, current_perm_seed));
        if ~isfolder(seed_folder), mkdir(seed_folder); end
        
        histories = cell(1, 4);
        pending_conditions = [];
        
        % Load Checkpoints
        for g = 1:4
            chk_file = fullfile(seed_folder, sprintf('Checkpoint_Opt_G%d.mat', g));
            histories{g} = empty_history;
            if resume_checkpoint && isfile(chk_file)
                chk_data = load(chk_file, 'historial');
                histories{g} = chk_data.historial;
            end
            if ~any(histories{g}.Target_Reached)
                pending_conditions(end+1) = g;
            end
        end
        
        fprintf('\n----------------------------------------------------\n');
        fprintf('PHASE 1: MULTI-OBJECTIVE GRID EXPLORATION\n');
        fprintf('----------------------------------------------------\n');
        
        if isempty(pending_conditions)
            fprintf('All recordings (G1-G4) have already reached their target.\n');
        else
            % Shuffle grid predictably
            grid_stream = RandStream('mt19937ar', 'Seed', 98000 + 100 * seed);
            shuffled_grid = full_grid(randperm(grid_stream, n_combinations), :);
            
            % Remove already evaluated combinations
            if ~isempty(histories{pending_conditions(1)})
                evaluated_params = [histories{pending_conditions(1)}.J_E, histories{pending_conditions(1)}.J_I, ...
                                    histories{pending_conditions(1)}.a_E, histories{pending_conditions(1)}.b_E, ...
                                    histories{pending_conditions(1)}.d_E];
                is_evaluated = ismembertol(shuffled_grid, evaluated_params, 1e-6, 'ByRows', true);
                shuffled_grid(is_evaluated, :) = [];
            end
            
            batch_size = cfg.search.parfor_batch_size;
            n_remaining = size(shuffled_grid, 1);
            evals_done = height(histories{pending_conditions(1)});
            
            for i_batch = 1:batch_size:n_remaining
                if isempty(pending_conditions) || evals_done >= cfg.search.max_grid_evaluations
                    break;
                end
                
                batch_end = min(i_batch + batch_size - 1, n_remaining);
                current_batch = shuffled_grid(i_batch:batch_end, :);
                n_batch = size(current_batch, 1);
                batch_results = cell(n_batch, 1);
                
                conds_to_evaluate = pending_conditions; 
                
                parfor in = 1:n_batch
                    theta_v = current_batch(in, :);
                    worker_stream = RandStream('mrg32k3a', 'Seed', 98000 + 100 * seed + i_batch + in);
                    RandStream.setGlobalStream(worker_stream);
                    
                    batch_results{in} = evaluate_grid_point(theta_v, seed, conds_to_evaluate, ...
                                                            pos_model, pos_original, targets, cfg, false);
                end
                
                evals_done = evals_done + n_batch;
                
                % Unpack batch results
                for in = 1:n_batch
                    theta_v = current_batch(in, :);
                    res_multi = batch_results{in};
                    
                    for g = conds_to_evaluate
                        if ismember(g, pending_conditions)
                            res = res_multi{g};
                            new_row = {theta_v(1), theta_v(2), theta_v(3), theta_v(4), theta_v(5), ...
                                res.Freq_Hz, res.Error_Freq_Norm, res.Error_Activity, res.Failed_Filters, ...
                                res.Passes_Filters, res.Classified, res.Start_Class_s, res.End_Class_s, ...
                                res.Error_L1_pp, res.Pct_Class(1), res.Pct_Class(2), res.Pct_Class(3), res.Pct_Class(4), ...
                                res.Score, res.Target_Reached, res.Valid_Simulation, res.Status, res.Message};
                                
                            histories{g} = [histories{g}; new_row];
                            
                            if res.Target_Reached
                                fprintf('*** GRID: G%d reached L1 target with score %.4f ***\n', g, res.Score);
                                pending_conditions(pending_conditions == g) = []; 
                            end
                        end
                    end
                end
                
                % Save checkpoints
                for g = conds_to_evaluate
                    chk_file = fullfile(seed_folder, sprintf('Checkpoint_Opt_G%d.mat', g));
                    historial = histories{g}; 
                    save(chk_file, 'historial', 'cfg', 'seed', 'g', '-v7.3');
                end
                
                fprintf('Grid Progress: %d / %d max evaluations.\n', evals_done, cfg.search.max_grid_evaluations);
            end
        end
        
        %% --- PHASE 2: RE-RUN AND EXPORT THE BEST POINT OF THE GRID ---
        fprintf('\n----------------------------------------------------\n');
        fprintf('PHASE 2: FINAL RE-EXECUTION AND EXPORT\n');
        fprintf('----------------------------------------------------\n');
        
        best_rows = cell(1, 4);
        best_T = table();
        
        % Final Data Containers
        spks_conds = cell(1,4); R_conds = cell(1,4); p_conds = cell(1,4);
        counts_conds = cell(1,4); freqs_conds = cell(1,4);
        class_counts = nan(4,4); class_pcts = nan(4,4);
        final_status = strings(1,4);
        
        for g = 1:4
            % Select best row
            T = histories{g};
            cands = find(T.Passes_Filters & T.Classified & isfinite(T.Error_L1_pp));
            if isempty(cands), cands = find(T.Passes_Filters & T.Valid_Simulation); end
            if isempty(cands), cands = find(T.Valid_Simulation); end
            if isempty(cands), cands = (1:height(T)).'; end
            [~, best_idx] = min(T.Score(cands));
            best_rows{g} = T(cands(best_idx), :);
            row = best_rows{g};
            
            fprintf('Best G%d: J_E=%.5g, J_I=%.5g, a=%.5g, b=%.5g, d=%.5g | %s | score=%.4f\n', ...
                g, row.J_E, row.J_I, row.a_E, row.b_E, row.d_E, row.Status, row.Score);
                
            % Re-execute with details
            theta = [row.J_E, row.J_I, row.a_E, row.b_E, row.d_E];
            fprintf('Simulating final model for G%d...\n', g);
            
            res = evaluate_grid_point(theta, seed, g, pos_model, pos_original, targets, cfg, true);
            res = res{g}; % Unpack single condition
            
            if res.Valid_Simulation
                spks_conds{g} = res.spks; R_conds{g} = res.R; p_conds{g} = res.p;
                counts_conds{g} = res.neuron_counts; freqs_conds{g} = res.neuron_freqs;
                class_counts(:,g) = res.Counts_Class; class_pcts(:,g) = res.Pct_Class;
                final_status(g) = res.Status;
            else
                spks_conds{g} = zeros(0,2); R_conds{g} = struct(); p_conds{g} = struct();
                counts_conds{g} = []; freqs_conds{g} = [];
                final_status(g) = "re_execution_error";
            end
            
            % Compile final table row
            new_row = {seed, g, theta(1), theta(2), theta(3), theta(4), theta(5), ...
                       res.Freq_Hz, res.Error_Freq_Norm, res.Error_Activity, res.Failed_Filters, ...
                       res.Passes_Filters, res.Classified, res.Start_Class_s, res.End_Class_s, ...
                       res.Error_L1_pp, res.Pct_Class(1), res.Pct_Class(2), res.Pct_Class(3), res.Pct_Class(4), ...
                       res.Score, res.Target_Reached, res.Valid_Simulation, res.Status, res.Message};
                       
            if isempty(best_T)
                best_T = cell2table(new_row, 'VariableNames', ['Semilla_catalogo', 'Grabacion', empty_history.Properties.VariableNames]);
            else
                best_T = [best_T; new_row];
            end
        end
        
        % Unpack structures for saving backward compatibility
        spks_1 = spks_conds{1}; spks_2 = spks_conds{2}; spks_3 = spks_conds{3}; spks_4 = spks_conds{4}; 
        p_1 = p_conds{1}; p_2 = p_conds{2}; p_3 = p_conds{3}; p_4 = p_conds{4}; 
        
        % Extract connectivity matrices from first valid R
        valid_Rs = find(~cellfun(@isempty, R_conds));
        if ~isempty(valid_Rs)
            base_R = R_conds{valid_Rs(1)};
            A = base_R.A; tipoNeurona = base_R.tipoNeurona; idxE = base_R.idxE; idxI = base_R.idxI;
            esE = base_R.esE; esI = base_R.esI; tipoCatalogo = base_R.tipoCatalogo;
            idxE_catalogo = base_R.idxE_catalogo; idxI_catalogo = base_R.idxI_catalogo;
            catalogoEnPosicion = base_R.catalogoEnPosicion; posicionPorCatalogo = base_R.posicionPorCatalogo;
        else
            A=[]; tipoNeurona=[]; idxE=[]; idxI=[]; esE=[]; esI=[]; tipoCatalogo=[]; 
            idxE_catalogo=[]; idxI_catalogo=[]; catalogoEnPosicion=[]; posicionPorCatalogo=[];
        end
        
        W_condiciones = cell(1,4); WE_condiciones = cell(1,4); WI_condiciones = cell(1,4);
        for g = 1:4
            if isfield(R_conds{g}, 'W'), W_condiciones{g} = R_conds{g}.W; end
            if isfield(R_conds{g}, 'WE'), WE_condiciones{g} = R_conds{g}.WE; end
            if isfield(R_conds{g}, 'WI'), WI_condiciones{g} = R_conds{g}.WI; end
        end
        W_1 = W_condiciones{1}; W_2 = W_condiciones{2}; W_3 = W_condiciones{3}; W_4 = W_condiciones{4};
        WE_1 = WE_condiciones{1}; WE_2 = WE_condiciones{2}; WE_3 = WE_condiciones{3}; WE_4 = WE_condiciones{4};
        WI_1 = WI_condiciones{1}; WI_2 = WI_condiciones{2}; WI_3 = WI_condiciones{3}; WI_4 = WI_condiciones{4};
        
        condicion_pesos_principal = valid_Rs(1);
        W = W_condiciones{condicion_pesos_principal}; 
        WE = WE_condiciones{condicion_pesos_principal}; 
        WI = WI_condiciones{condicion_pesos_principal};
        
        % Structure arrays for MAT file
        combinaciones_completas = [best_T.J_E, best_T.J_I, best_T.a_E, best_T.b_E, best_T.d_E]; 
        combinaciones = combinaciones_completas(:, 3:5); 
        nombres_parametros = cfg.names; 
        
        T_combinaciones = array2table(combinaciones_completas, 'VariableNames', cellstr(nombres_parametros)); 
        T_combinaciones = addvars(T_combinaciones, (1:4).', 'Before', 1, 'NewVariableNames', 'Grabacion');
            
        modelo = struct('version', cfg.version, 'funcion_modelo', 'red_izh_catalogo_permutado', ...
            'tipo_seed_optimizada', 'catalogo_heterogeneidad', 'semilla_identidad', seed, ...
            'semilla_catalogo', seed, 'semilla_permutacion', cfg.model.seed_permutation, ...
            'semilla_conectividad', cfg.model.seed_connectivity, 'semilla_dinamica', cfg.model.seed_dynamics, ...
            'combinaciones', combinaciones, 'combinaciones_completas', combinaciones_completas, ...
            'T_combinaciones', T_combinaciones, 'nombres_parametros', nombres_parametros, ...
            'parametros', {p_conds}, 'estado', final_status, 'duracion_s', cfg.model.T_ms / 1000, ...
            'ventana_metricas_s', cfg.model.metrics_window_s, 'duracion_clasificacion_s', cfg.classification.duration_s, ...
            'condicion_pesos_principal', condicion_pesos_principal); 
            
        % Formatting variables specifically named for downward compatibility
        file_length = sprintf('%g min', cfg.model.T_ms / 60000); 
        stim_time = sprintf('%.6g min', 10000 / 60000); 
        ganglion_imaged = 'simulated network'; p9_stim = 'model'; 
        objetivos_experimentales = targets; rutas_datos_experimentales = []; 
        semilla_catalogo = seed; semilla_permutacion = cfg.model.seed_permutation; semilla_identidad = seed; 
        fecha_optimizacion = datetime('now'); fuente_posiciones = char(positions_file); 
        version_optimizador = cfg.version; tipos = cfg.classification.types;
        
        % Rename for output matching 
        T_mejores = best_T; 
        estado_final = final_status; 
        R_condiciones = R_conds; 
        historiales = histories;
        conteos_clasificacion = class_counts; 
        porcentajes_modelo = class_pcts;
        conteos_neuronales_condiciones = counts_conds; 
        frecuencia_neuronal_hz_condiciones = freqs_conds;
        
        final_file = fullfile(models_folder, sprintf('Modelo_IZH_seed_%03d_perm_%03d.mat', seed, current_perm_seed));
            
        save(final_file, 'spks_1', 'spks_2', 'spks_3', 'spks_4', 'x', 'xstd', 'y', 'ystd', ...
            'file_length', 'stim_time', 'ganglion_imaged', 'p9_stim', 'A', 'W', 'WE', 'WI', ...
            'W_1', 'W_2', 'W_3', 'W_4', 'WE_1', 'WE_2', 'WE_3', 'WE_4', 'WI_1', 'WI_2', 'WI_3', 'WI_4', ...
            'W_condiciones', 'WE_condiciones', 'WI_condiciones', 'tipoNeurona', 'idxE', 'idxI', ...
            'esE', 'esI', 'tipoCatalogo', 'idxE_catalogo', 'idxI_catalogo', 'catalogoEnPosicion', ...
            'posicionPorCatalogo', 'p_1', 'p_2', 'p_3', 'p_4', 'R_condiciones', ...
            'conteos_neuronales_condiciones', 'frecuencia_neuronal_hz_condiciones', 'historiales', ...
            'T_mejores', 'conteos_clasificacion', 'porcentajes_modelo', 'objetivos_experimentales', ...
            'rutas_datos_experimentales', 'estado_final', 'cfg', 'semilla_identidad', 'semilla_catalogo', ...
            'semilla_permutacion', 'fecha_optimizacion', 'fuente_posiciones', 'version_optimizador', ...
            'tipos', 'modelo', 'combinaciones', 'combinaciones_completas', 'T_combinaciones', ...
            'nombres_parametros', 'condicion_pesos_principal', '-v7.3');
            
        final_files_list(end+1, 1) = string(final_file); 
        catalog_seed_record(end+1, 1) = seed;
        perm_seed_record(end+1, 1) = current_perm_seed;
        
        fprintf('\nSeed %d (Perm %d) final optimization saved.\n', seed, current_perm_seed);
    end
end

fprintf('\nFINAL EXPORTED FILES\n');
disp(final_files_list);

T_archivos = table(catalog_seed_record, perm_seed_record, final_files_list, ...
    'VariableNames', {'Semilla_catalogo', 'Semilla_permutacion', 'Archivo_MAT'}); 
index_file = fullfile(output_folder, 'Optimized_Models_Index.mat');
save(index_file, 'T_archivos', 'final_files_list', 'heterogeneity_seeds', 'permutation_seeds', 'cfg', '-v7.3');
