function c = get_project_colors(color_name)
% GET_PROJECT_COLORS Returns standard RGB vectors for consistent plotting.

    switch lower(color_name)
        case {'blue', 'control', 'g1', 'g4'}
            c = [0.20 0.60 0.80];
        case {'orange', 'serotonin', 'g2', 'g3'}
            c = [0.80 0.40 0.20];
        case {'non_oscillatory', 'non'}
            c = [46 204 113] / 255; % Green
        case 'oscillator'
            c = [52 152 219] / 255; % Blue
        case 'burster'
            c = [243 156 18] / 255; % Orange/Yellow
        case 'pauser'
            c = [142 68 173] / 255; % Purple
        case 'inhibitory'
            c = [0.0000 0.4470 0.7410];
        case 'excitatory'
            c = [0.8500 0.3250 0.0980];
        otherwise
            c = [0.5 0.5 0.5]; % Gray default
    end
end