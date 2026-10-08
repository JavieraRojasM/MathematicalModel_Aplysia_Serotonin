function Classification_Bruno = RunBrunoAnalysis(DataList, datapath, outputFolder, analysis, outputName)
    % Internal setup for the toolbox
    addpath(outputFolder);
    nfiles = size(DataList,2);
    
    if analysis == 1
        disp('Step: Spike Train Properties...');
        Analyse_Spike_Train_Properties_MOD;
    
    else
        
    end
    %% Clasificación
    % --- Paso 1: Carga de Datos ---
    data = load(fullfile(outputFolder, 'Analyses_Neurons_and_Ensembles.mat'), 'neurondata');
    neurons_all = data.neurondata;
    
    % --- Paso 2: Clasificación Individual (Criterio de Bins Consecutivos) ---
    for i = 1:numel(neurons_all)
        % Limpieza de vacíos
        if isempty(neurons_all(i).blnOscN); neurons_all(i).blnOscN = 2; end
        if isempty(neurons_all(i).blnOscP); neurons_all(i).blnOscP = 2; end
        
        % Classification
        if neurons_all(i).blnOscP == 0 && neurons_all(i).blnOscN == 0
            neurons_all(i).Type = 1; % Non-oscillatory
        elseif neurons_all(i).blnOscP == 1 && neurons_all(i).blnOscN == 1
            neurons_all(i).Type = 2; % Oscillator
        elseif neurons_all(i).blnOscP == 1 && neurons_all(i).blnOscN == 0
            neurons_all(i).Type = 3; % Burster
        elseif neurons_all(i).blnOscP == 0 && neurons_all(i).blnOscN == 1
            neurons_all(i).Type = 4; 
        else
            neurons_all(i).Type = 5; % Unknown/Empty (los "2")
        end 
    end
  
    % --- Paso 3: Distribución por Grabación y Categoría ---
    categories = {'Non-oscillatory', 'Oscillator', 'Burster', 'Pauser'};
    type_codes = [1, 2, 3, 4];
    Classification_Bruno = struct();
    
    % Extraemos vectores para filtrar rápidamente (acelera el proceso)
    all_recs = [neurons_all.Recording];
    all_types = [neurons_all.Type];
    all_ids = [neurons_all.ID];
    
    unique_recs = unique(all_recs);
    num_recs = numel(unique_recs);
    summary_counts = zeros(4, num_recs);
    
    % Ciclo para procesar cada grabación
    for r = 1:num_recs
        rec_id = unique_recs(r);
        
        % Filtramos los índices de las neuronas que pertenecen a esta grabación
        idx_this_rec = (all_recs == rec_id);
        
        ids_col = cell(4, 1); 
        counts_col = zeros(4, 1);
        
        for c = 1:4
            % Buscamos dentro de la grabación las que coinciden con el tipo
            match_type = (all_types(idx_this_rec) == type_codes(c));
            temp_ids = all_ids(idx_this_rec);
            temp_list = temp_ids(match_type); 
            
            ids_col{c} = temp_list; % Guardamos lista de IDs
            counts_col(c) = length(temp_list);
            summary_counts(c, r) = length(temp_list); 
        end
        
        % Guardamos tabla por grabación
        field_name = sprintf('spks_%d', rec_id);
        Classification_Bruno.(field_name) = table(categories', ids_col, counts_col, ...
            'VariableNames', {'Neuron_type', 'IDs', 'Total_Count'});
    end
    
    % --- Paso 4: Tabla Resumen Final ---
    summary_varnames = cellfun(@(x) sprintf('Count_spks%d', x), num2cell(unique_recs), 'UniformOutput', false);
    
    SummaryTable = array2table(summary_counts, 'VariableNames', summary_varnames);
    SummaryTable = [table(categories', 'VariableNames', {'Neuron_type'}), SummaryTable];
    
    Classification_Bruno.Summary = SummaryTable;
    
    % --- Paso 5: Guardado de Resultados ---
    output_file = fullfile(outputFolder, outputName);
    save(output_file, 'Classification_Bruno');
    
end