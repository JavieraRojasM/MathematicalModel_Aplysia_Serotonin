%% SPATIAL MAPS (DATA AND MODEL)
% Generates a spatial map for each of the four conditions.
%
% EXPERIMENTAL DATA:
%   - Gray: Active neuron (without serotonin)
%   - Orange: Active neuron (with serotonin)
%   - White: Inactive/Absent neuron
%
% IZHIKEVICH MODEL:
%   - Blue: Inhibitory neuron or connection
%   - Orange: Excitatory neuron or connection
%   - Solid/Intense color: Active neuron (fired in the condition)
%   - Light/Faded color: Inactive neuron (no spikes in the condition)

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
% ================================================================

%% 3. FIGURE CONFIGURATION
n_conditions = 4;
marker_size = 320;
show_labels = true;
invert_y_axis = true;

% false: considers a neuron present if it fires at any point.
% true: only considers spikes within the specified time window.
use_time_window = false;

% false: shows full model topology.
% true: shows only connections where both pre- and post-synaptic neurons fired.
only_active_connections = false;

fig_width_cm = 15;
fig_height_cm = 15;

% Figure export settings controlled by auxiliary function
save_png = true;
save_svg = true;
save_fig = false;

export_formats = {'pdf'};
if save_png, export_formats{end+1} = 'png'; end
if save_svg, export_formats{end+1} = 'svg'; end
if save_fig, export_formats{end+1} = 'fig'; end

% Common limits (Use [] for automatic limits)
x_limits = [0 45];
y_limits = [-1 50];

% USE AUXILIARY FUNCTION: get_project_colors
color_inhibitory = get_project_colors('inhibitory');
color_excitatory = get_project_colors('excitatory');
color_with_5HT = get_project_colors('serotonin');
color_no_5HT = [0.70 0.70 0.70];
color_absent = [1.00 1.00 1.00];
color_edge = [0.10 0.10 0.10];
conds_with_5HT = [2 3]; % Conditions 2 and 3 have serotonin

% Light versions for inactive model neurons (blended with 72% white)
color_inhibitory_light = (1 - 0.72) * color_inhibitory + 0.72 * [1 1 1];
color_excitatory_light = (1 - 0.72) * color_excitatory + 0.72 * [1 1 1];

connection_thickness = 0.65;

%% 4. SOURCE CONFIGURATION
switch lower(mode_type)
    case 'data'
        input_path = data_path;
        start_time_s = 620;
        end_time_s = 860;
        condition_names = {'No serotonin'; 'Serotonin injection'; ...
                           'Maintained serotonin'; 'Washout, no serotonin'};
        output_folder = fullfile(project_folder, 'Spatial_Maps_Results', 'Data');
    case 'model'
        input_path = model_path;
        start_time_s = 30;
        end_time_s = 270;
        [~, model_name] = fileparts(model_path);
        model_name = regexprep(model_name, '[^a-zA-Z0-9_-]', '_');
        output_folder = fullfile(project_folder, 'Spatial_Maps_Results', 'Model', model_name);
    otherwise
        error("mode_type must be 'data' or 'model'.");
end

if ~isfile(input_path), error('File not found: %s', input_path); end
if ~isfolder(output_folder), mkdir(output_folder); end

%% 5. LOAD AND VALIDATE DATA
data = load(input_path);
spike_fields = cell(n_conditions, 1);
for cond = 1:n_conditions
    spike_fields{cond} = sprintf('spks_%d', cond);
    if ~isfield(data, spike_fields{cond})
        error('The file is missing the variable %s.', spike_fields{cond});
    end
end

if ~isfield(data, 'x') || ~isfield(data, 'y')
    error('The file must contain variables x and y.');
end

x = double(data.x(:));
y = double(data.y(:));
n_neurons = numel(x);
IDs = (1:n_neurons)';
valid_position = isfinite(x) & isfinite(y);

% Load topology if model
combinations = nan(n_conditions, 3);
if strcmpi(mode_type, 'model')
    if isfield(data, 'A'), A = logical(data.A); 
    elseif isfield(data, 'W'), A = (data.W ~= 0); 
    else, error('Model must contain connectivity matrix A or W.'); end
    A(1:n_neurons+1:end) = false; % Remove self-connections
    
    if isfield(data, 'idxE') && isfield(data, 'idxI')
        idxE = double(data.idxE(:));
        idxI = double(data.idxI(:));
    elseif isfield(data, 'tipoNeurona')
        idxE = find(startsWith(upper(string(data.tipoNeurona)), "E"));
        idxI = find(startsWith(upper(string(data.tipoNeurona)), "I"));
    else
        error('Model must contain idxE/idxI or tipoNeurona.');
    end
    
    neuron_type = repmat("E", n_neurons, 1);
    neuron_type(idxI) = "I";
    
    if isfield(data, 'modelo') && isfield(data.modelo, 'combinaciones')
        combinations = double(data.modelo.combinaciones);
        condition_names = cellstr(compose('C%d: a=%.4g, b=%.4g, d=%.4g', (1:n_conditions)', combinations));
    else
        condition_names = cellstr(compose('Condition %d', (1:n_conditions)'));
    end
else
    A = false(n_neurons);
    idxE = []; idxI = [];
    neuron_type = repmat("N/A", n_neurons, 1);
end

%% 6. GENERATE SPATIAL MAPS
spatial_tables = cell(n_conditions, 1);

for cond = 1:n_conditions
    spikes_raw = data.(spike_fields{cond});
    
    % USE AUXILIARY FUNCTION: format_spikes
    if use_time_window
        [valid_spikes, ~] = format_spikes(spikes_raw, start_time_s, end_time_s, n_neurons);
    else
        [valid_spikes, ~] = format_spikes(spikes_raw, [], [], n_neurons);
    end
    
    if isempty(valid_spikes)
        spike_IDs = [];
        spike_counts = zeros(n_neurons, 1);
    else
        spike_IDs = valid_spikes(:, 1);
        spike_counts = accumarray(spike_IDs, 1, [n_neurons, 1], @sum, 0);
    end
    
    active_IDs = unique(spike_IDs);
    is_present = ismember(IDs, active_IDs);

    %% 6.1 Spatial Table
    current_table = table(IDs, x, y, neuron_type, valid_position, is_present, spike_counts, ...
        'VariableNames', {'Neuron', 'X', 'Y', 'NeuronType', 'ValidPosition', 'Present', 'SpikeCount'});
    
    sort_matrix = [-y, x];
    if invert_y_axis, sort_matrix = [y, x]; end
    sort_matrix(~isfinite(sort_matrix)) = Inf;
    [~, spatial_order] = sortrows(sort_matrix, [1 2]);
    
    current_table = current_table(spatial_order, :);
    current_table.SpatialOrder = (1:height(current_table))';
    current_table = movevars(current_table, 'SpatialOrder', 'Before', 'Neuron');
    spatial_tables{cond} = current_table;

    %% 6.2 Plot Figure
    fig = figure('Color', 'w', 'Name', sprintf('Spatial Map - Condition %d', cond), 'Visible', 'off');
    ax = axes(fig); hold(ax, 'on');
    
    % --- PLOT CONNECTIONS (Model Only) ---
    if strcmpi(mode_type, 'model')
        [post, pre] = find(A);
        if only_active_connections
            valid_conn = is_present(post) & is_present(pre);
            post = post(valid_conn); pre = pre(valid_conn);
        end
        
        is_E = ismember(pre, idxE);
        pre_E = pre(is_E); post_E = post(is_E);
        pre_I = pre(~is_E); post_I = post(~is_E);
        
        if ~isempty(pre_E)
            plot(ax, [x(pre_E)'; x(post_E)'; nan(1, length(pre_E))], ...
                     [y(pre_E)'; y(post_E)'; nan(1, length(pre_E))], ...
                 'Color', color_excitatory, 'LineWidth', connection_thickness, 'HandleVisibility', 'off');
        end
        if ~isempty(pre_I)
            plot(ax, [x(pre_I)'; x(post_I)'; nan(1, length(pre_I))], ...
                     [y(pre_I)'; y(post_I)'; nan(1, length(pre_I))], ...
                 'Color', color_inhibitory, 'LineWidth', connection_thickness, 'HandleVisibility', 'off');
        end
    end

    % --- PLOT NEURONS (Vectorized Scatter) ---
    node_colors = repmat(color_absent, n_neurons, 1);
    text_colors = repmat([0.1 0.1 0.1], n_neurons, 1);
    
    if strcmpi(mode_type, 'model')
        node_colors(idxE, :) = repmat(color_excitatory_light, length(idxE), 1);
        node_colors(idxI, :) = repmat(color_inhibitory_light, length(idxI), 1);
        
        active_E = intersect(idxE, find(is_present));
        active_I = intersect(idxI, find(is_present));
        
        node_colors(active_E, :) = repmat(color_excitatory, length(active_E), 1);
        node_colors(active_I, :) = repmat(color_inhibitory, length(active_I), 1);
        
        text_colors(is_present, :) = repmat([1 1 1], sum(is_present), 1);
    else
        if ismember(cond, conds_with_5HT)
            node_colors(is_present, :) = repmat(color_with_5HT, sum(is_present), 1);
        else
            node_colors(is_present, :) = repmat(color_no_5HT, sum(is_present), 1);
        end
        text_colors(is_present, :) = repmat([1 1 1], sum(is_present), 1);
    end
    
    scatter(ax, x(valid_position), y(valid_position), marker_size, ...
            node_colors(valid_position, :), 'filled', 'MarkerEdgeColor', color_edge, 'LineWidth', 1.2);
            
    % --- PLOT LABELS ---
    if show_labels
        for n = 1:n_neurons
            if valid_position(n)
                text(ax, x(n), y(n), num2str(n), 'FontSize', 8, 'FontWeight', 'bold', ...
                     'Color', text_colors(n, :), 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'Clipping', 'on');
            end
        end
    end
    
    % --- FORMATTING ---
    axis(ax, 'equal'); axis(ax, 'tight');
    if ~isempty(x_limits), xlim(ax, x_limits); end
    if ~isempty(y_limits), ylim(ax, y_limits); end
    if invert_y_axis, set(ax, 'YDir', 'reverse'); end
    box(ax, 'on'); yticks(ax, []); xticks(ax, []); hold(ax, 'off');
    
    % USE AUXILIARY FUNCTION: export_figure
    base_name = fullfile(output_folder, sprintf('Spatial_Map_Condition_%d', cond));
    export_figure(fig, base_name, fig_width_cm, fig_height_cm, export_formats, 600);
    close(fig);
    
    fprintf('Condition %d: %d of %d neurons present.\n', cond, nnz(is_present), n_neurons);
end

%% 7. SAVE SUMMARY DATA
save(fullfile(output_folder, 'Spatial_Maps_Results.mat'), ...
    'spatial_tables', 'mode_type', 'input_path', 'spike_fields', ...
    'condition_names', 'combinations', 'use_time_window', ...
    'start_time_s', 'end_time_s', 'only_active_connections', ...
    'idxE', 'idxI', 'neuron_type', 'color_inhibitory', 'color_excitatory');

fprintf('\nSpatial Maps saved in:\n%s\n', output_folder);