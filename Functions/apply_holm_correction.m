function [p_adjusted, significance] = apply_holm_correction(p_original)
% APPLY_HOLM_CORRECTION Applies Holm-Bonferroni correction and assigns significance stars.

    p_adjusted = nan(size(p_original));
    valid_idx = find(isfinite(p_original));
    
    if ~isempty(valid_idx)
        [p_sorted, sort_idx] = sort(p_original(valid_idx));
        m = length(p_sorted);
        p_holm = p_sorted .* (m:-1:1)';
        p_holm = min(1, cummax(p_holm)); % cummax ensures monotonicity
        p_adjusted(valid_idx(sort_idx)) = p_holm;
    end
    
    % Assign significance stars
    significance = strings(length(p_original), 1);
    for comp = 1:length(p_original)
        p = p_adjusted(comp);
        if isnan(p)
            significance(comp) = "NA";
        elseif p < 0.001
            significance(comp) = "***";
        elseif p < 0.01
            significance(comp) = "**";
        elseif p < 0.05
            significance(comp) = "*";
        else
            significance(comp) = "ns";
        end
    end
end