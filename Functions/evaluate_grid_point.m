function results = evaluate_grid_point(theta, seed, conds_to_evaluate, pos, pos_original, targets, cfg, save_details)
    results = cell(1, 4);
    
    % Initialize all conditions with empty templates
    empty_res = struct('Valid_Simulation', false, 'Freq_Hz', NaN, 'Error_Freq_Norm', NaN, 'Error_Activity', NaN, ...
                       'Failed_Filters', 1, 'Passes_Filters', false, 'Classified', false, 'Start_Class_s', NaN, ...
                       'End_Class_s', NaN, 'Error_L1_pp', NaN, 'Score', 9999, 'Target_Reached', false, ...
                       'Status', "pending", 'Message', "", 'Counts_Class', nan(4,1), 'Pct_Class', nan(4,1), ...
                       'spks', zeros(0,2), 'R', struct(), 'p', struct(), 'neuron_counts', [], 'neuron_freqs', []);
    for g = 1:4, results{g} = empty_res; end
    
    try
        % Build base parameters
        p = get_base_params(pos, pos_original, size(pos, 1), seed, cfg);
        p.J_E = theta(1); p.J_I = theta(2); p.a_E = theta(3); p.b_E = theta(4); p.d_E = theta(5);
        p.T = cfg.model.T_prefilters_ms;
        if save_details, p.T = cfg.model.T_ms; end
        
        % Run base simulation
        R = red_izh_catalogo_permutado(p);
        spks_raw = [double(R.spikeNeuron(:)), double(R.spikeTimes(:)) / 1000];
        
        % --- USO DE LA FUNCION AUXILIAR: format_spikes ---
        window = cfg.model.metrics_window_s;
        [valid_spikes, ~] = format_spikes(spks_raw, window(1), window(2), R.N);
        
        if isempty(valid_spikes)
            neuron_counts = zeros(R.N, 1);
        else
            neuron_counts = accumarray(valid_spikes(:, 1), 1, [R.N 1], @sum, 0);
        end
        neuron_freqs = neuron_counts / diff(window);
        freq_hz_val = median(neuron_freqs);
        
        any_passed = false;
        
        for g = conds_to_evaluate
            res = results{g};
            res.Valid_Simulation = true;
            res.Freq_Hz = freq_hz_val;
            res.neuron_counts = neuron_counts;
            res.neuron_freqs = neuron_freqs;
            
            % Check Frequency Filters
            tol_f = max(cfg.filters.abs_freq_tolerance_hz, abs(targets.freq_hz(g)) * cfg.filters.rel_freq_tolerance);
            res.Error_Freq_Norm = abs(res.Freq_Hz - targets.freq_hz(g)) / tol_f;
            res.Error_Activity = res.Error_Freq_Norm;
            res.Failed_Filters = sum(res.Error_Freq_Norm > 1);
            res.Passes_Filters = (res.Failed_Filters == 0);
            
            if res.Passes_Filters, any_passed = true; end
            results{g} = res;
        end
        
        % If any condition passed the frequency filter, classify spikes
        if any_passed
            start_c = cfg.model.metrics_window_s(1);
            end_c = start_c + cfg.classification.duration_s;
            
            if save_details
                R_ext = R;
                % Para los detalles completos usamos format_spikes sin recortar ventana de tiempo inicial
                [spks_class, ~] = format_spikes(spks_raw, [], [], R.N);
            else
                % Re-run full simulation if it passed
                p.T = cfg.model.T_ms;
                R_ext = red_izh_catalogo_permutado(p);
                spks_full = [double(R_ext.spikeNeuron(:)), double(R_ext.spikeTimes(:)) / 1000];
                [spks_class, ~] = format_spikes(spks_full, [], [], R_ext.N);
            end
            
            for g = conds_to_evaluate
                res = results{g};
                if res.Passes_Filters
                    res.Start_Class_s = start_c;
                    res.End_Class_s = end_c;
                    try
                        [count, pct] = CLASIFICAR_SPIKES_BRUNO_MULTIFILTRO(spks_class, theta, start_c, end_c, ...
                            cfg.classification.which_spks, seed, g, cfg.classification.seed_bruno(g), ...
                            cfg.classification.types, sprintf('Catalog_%03d_G%d', seed, g));
                            
                        res.Counts_Class = count(:);
                        res.Pct_Class = pct(:);
                        res.Error_L1_pp = sum(abs(pct(:) - targets.percentages(:,g)));
                        res.Classified = true;
                        res.Score = res.Error_L1_pp;
                        res.Target_Reached = res.Error_L1_pp <= cfg.filters.class_L1_tolerance_pp;
                        
                        if res.Target_Reached, res.Status = "target_reached"; else, res.Status = "classified"; end
                        if save_details, res.spks = spks_class; res.R = R_ext; res.p = p; end
                    catch ME_class
                        res.Score = 900 + min(res.Error_Activity, 99);
                        res.Status = "classification_error";
                        res.Message = string(ME_class.message);
                    end
                else
                    res.Score = 1000 + 100 * res.Failed_Filters + min(res.Error_Activity, 99);
                    res.Status = "failed_filters";
                    [res.spks, ~] = format_spikes(spks_raw, [], [], R.N);
                    res.R = R; res.p = p;
                end
                results{g} = res;
            end
        else
            for g = conds_to_evaluate
                res = results{g};
                res.Score = 1000 + 100 * res.Failed_Filters + min(res.Error_Activity, 99);
                res.Status = "failed_filters";
                [res.spks, ~] = format_spikes(spks_raw, [], [], R.N);
                res.R = R; res.p = p;
                results{g} = res;
            end
        end
        
    catch ME
        for g = conds_to_evaluate
            res = results{g};
            res.Valid_Simulation = false;
            res.Score = 9999;
            res.Status = "simulation_error";
            res.Message = string(ME.message);
            results{g} = res;
        end
    end
end