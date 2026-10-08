function export_figure(fig, base_path, width_cm, height_cm, formats, dpi)
% EXPORT_FIGURE Exports a MATLAB figure to multiple formats with precise dimensions.

    if nargin < 6 || isempty(dpi)
        dpi = 300;
    end
    if nargin < 5 || isempty(formats)
        formats = {'pdf', 'png', 'svg'};
    end
    
    set(fig, 'Color', 'w');
    set(fig, 'Units', 'centimeters', 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [width_cm, height_cm], 'PaperPosition', [0, 0, width_cm, height_cm]);
    
    pos = get(fig, 'Position');
    set(fig, 'Position', [pos(1), pos(2), width_cm, height_cm]);
    drawnow;
    
    for i = 1:length(formats)
        fmt = lower(formats{i});
        file_name = sprintf('%s.%s', base_path, fmt);
        try
            switch fmt
                case 'pdf'
                    exportgraphics(fig, file_name, 'ContentType', 'vector', 'BackgroundColor', 'white');
                case 'svg'
                    exportgraphics(fig, file_name, 'ContentType', 'vector', 'BackgroundColor', 'white');
                case 'png'
                    exportgraphics(fig, file_name, 'Resolution', dpi, 'BackgroundColor', 'white');
                case 'fig'
                    savefig(fig, file_name);
            end
        catch ME
            warning('Could not save %s: %s', file_name, ME.message);
        end
    end
end