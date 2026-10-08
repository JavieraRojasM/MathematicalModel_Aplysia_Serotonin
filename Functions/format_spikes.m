function [spks_filtered, n_neurons] = format_spikes(spks_raw, min_time_s, max_time_s, expected_neurons)
% FORMAT_SPIKES Standardizes spike matrix to [ID, Time] and filters by time window.

    if isempty(spks_raw)
        spks_filtered = zeros(0, 2);
        n_neurons = max(1, expected_neurons);
        return;
    end

    spks = double(spks_raw);
    
    % Ensure N x 2 format
    if size(spks, 2) ~= 2 && size(spks, 1) == 2
        spks = spks.';
    end
    
    % Round IDs and handle 0-indexing
    spks(:, 1) = round(spks(:, 1));
    if min(spks(:,1)) == 0
        spks(:,1) = spks(:,1) + 1; 
    end
    
    % Filter by time window
    if nargin >= 3 && ~isempty(min_time_s) && ~isempty(max_time_s)
        in_window = spks(:, 2) >= min_time_s & spks(:, 2) <= max_time_s;
        spks_filtered = spks(in_window, :);
    else
        spks_filtered = spks;
    end
    
    % Determine total neurons
    if nargin < 4 || isempty(expected_neurons)
        expected_neurons = 0;
    end
    
    max_id = 0;
    if ~isempty(spks_filtered)
        max_id = max(spks_filtered(:, 1));
    end
    n_neurons = max(expected_neurons, max_id);
end