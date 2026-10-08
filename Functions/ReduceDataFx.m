%% Reduce Data

function of_name = ReduceDataFx(file_path, min_plottime, max_plottime, which_spks, TestData_Path)

%% -------------------- CONFIGURACIÓN --------------------
%file_path = fullfile('C:\Users\javie\OneDrive\Escritorio\Datos','Sep0622_serotonin.mat');
%file_path = fullfile('/home/javiera/Documents/MATLAB/Sep0622_serotonin.mat');

data = load(file_path);


if which_spks == 5
    pca_cmap_name = 'cool';
    title_extra = 'Izh';
    spks = data.spks; 
    spks(spks(:, 2) < min_plottime, :) = [];
    spks(spks(:, 2) > max_plottime, :) = [];
    
    
    of_name = sprintf('spks_%d_from%d_to%d', which_spks, min_plottime, max_plottime);
    
    
    %Se guarda la matriz en un archivo .mat
    save(fullfile(TestData_Path, of_name),'spks', "max_plottime", 'min_plottime', 'title_extra');

else
    
    if which_spks == 1
        pca_cmap_name = 'parula';
        title_extra = 'Sin Serotonina';
        if isfield(data,'spks_1')
            spks = data.spks_1; 
        else
            error('spks_1 no existe en el archivo.'); 
        end
    elseif which_spks == 2
        pca_cmap_name = 'hot';
        title_extra = 'Serotonina inyectada durante grabacion';
        if isfield(data,'spks_2')
            spks = data.spks_2; 
        else
            error('spks_2 no existe en el archivo.'); 
        end
    
    elseif which_spks == 3
        pca_cmap_name = 'spring';
        title_extra = 'Con Serotonina';
        if isfield(data,'spks_3')
            spks = data.spks_3; 
        else
            error('spks_3 no existe en el archivo.'); 
        end
    
    elseif which_spks == 4
        pca_cmap_name = 'cool';
        title_extra = 'Sin Serotonina, lavada antes de grabacion';
        if isfield(data,'spks_4')
            spks = data.spks_4; 
        else
            error('spks_4 no existe en el archivo.'); 
        end
    
    else
        error(fprintf('spks_%d no existe en el archivo.', which_spks)); 
    end
    
    
    data_serotonintime = "5 min";
    
    
    % Leer tiempo de estímulo y largo de archivo
    total_time   = sscanf(data.file_length, '%d')*60;
    stim_app   = sscanf(data.stim_time, '%d')*60;
    %ser_app = sscanf(data_serotonintime, '%d') * 60;    % Minutes to seconds conversion
    
    
    spks(spks(:, 2) < min_plottime, :) = [];
    spks(spks(:, 2) > max_plottime, :) = [];
    
    
    of_name = sprintf('spks_%d_from%d_to%d', which_spks, min_plottime, max_plottime);
    
    
    %Se guarda la matriz en un archivo .mat
    save(fullfile(TestData_Path, of_name), 'spks', "stim_app", "max_plottime", 'min_plottime', 'title_extra');

end 

end 