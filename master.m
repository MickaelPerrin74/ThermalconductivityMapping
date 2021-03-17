Params.N_workers = 19;
Params.power = 4;
Params.temperature = 317;
Params.name = '8um_S1_K12_post';
Params.absorption_file = '0.027';
Params.absorption_file_ref = '0.027';
Params.laser_spot_radius = 0.3;
Params.membrane_radius = 4.0;
Params.membrane_mask = 3.5;
spatial_averaging = 0;
Params.power_correction = 1;
Params.g = 28000000;
Params.convection = 29000;
distcomp.feature( 'LocalUseMpiexec', false )

%% clean up
system('rm tmp/status_*');
system('rm tmp/pid_*');
system('rm tmp/port_*');
system('rm tmp/logfile_*');
system('killall java');
system('rm tmp/worker_time*');

%%
mkdir('tmp');

%% parameters
T_max_error_conv = 5;    % K
T_mean_error_conv = 1;    % K
convergence_slope_kappa = 1.1;
convergence_slope_absorption = 0.0001;
addpath(genpath('functions'));
addpath(genpath('/opt/ud/comsol55/mli'));

%Params.g = 1.724e7;		%  Villaroman et al. Carbon 2017
Params.support_radius = 20.0;    % um
Params.grid_resolution = 1;
Params.display = 0;
Params.Timeout = 120;
Params.laser_spot_radius_tot = 2.5; % um
Params.export_data = 0;         % export each temperature map as dat file
Params.export_mph = 0;          % export each mph file

Kappa_min = 100;
Kappa_max = 4000;
Kappa_start = 1000;
Params.Kappa_support = 500;
Params.N_kappa = 50;
Params.N_absorption = 100;
Params.N_kappa_init = 100;

N_total = 200;

Absorption_min = 0.023;

%% set parameters
fid = fopen('logfile','w');
try
    load Data.mat
    fprintf('%s -- Resuming...', datetime('now'));
    fprintf(fid, '%s -- Resuming...', datetime('now'));
catch
    
    %% get T experimental
    fprintf('%s -- Getting experimental T...', datetime('now'));
    fprintf(fid, '%s -- Getting experimental T...', datetime('now'));
    load(sprintf('Input_data/Formatted data/%s/calculated_temperature_%1.2fmW_%1.0fK.mat', Params.name, Params.power, Params.temperature));
    X_array = X;
    Y_array = Y;
    
    fprintf('done\n');
    fprintf(fid, 'done\n');
    
    [N_y, N_x] = size(T_exp);
    
    %% load absorption
    fprintf('%s -- Getting absorption...', datetime('now'));
    fprintf(fid, '%s -- Getting absorption...', datetime('now'));
    
    Absorption = str2double(Params.absorption_file) * ones(size(T_exp));
    
    %% double absorption on support
    [X,Y] = meshgrid(X_array, Y_array);
    Absorption( sqrt(X.^2 + Y.^2) > Params.membrane_radius) = 2 * Absorption_min;
    
    Absorption_init = Absorption;
    
    fprintf('done\n');
    fprintf(fid, 'done\n');
    
    %% Intialize mask
    Kappa = zeros(size(X)) + Kappa_start;
    Kappa( sqrt(X.^2 + Y.^2) > Params.membrane_radius) = Params.Kappa_support;
    mask = zeros(size(Kappa));
    mask( sqrt(X.^2 + Y.^2) < Params.membrane_mask) = 1;
    Params.mask = mask;
    
    %% Initialize arrays
    Kappa_init = Kappa;
    
    dX = mean(diff(X_array));
    dY = mean(diff(Y_array));
       
    T_mean_error_min = 20000;
    T_max_error = 10000;
    T_mean_error = 10000;
    
    convergence_counter = 1;
    TEMP = cell(1,1);
    KAPPA = cell(1,1);
    ABSORPTION = cell(1,1);
    DIFF_KAPPA = cell(1,1);
    DIFF_ABSORPTION = cell(1,1);
    DIFF_KAPPA_SMOOTH = cell(1,1);
    DIFF_ABSORPTION_SMOOTH = cell(1,1);
    T_Max_Error = zeros(1);
    T_Mean_Error = zeros(1);
    
end

%% clear persistent variable
clear comsol_model

%% parrallel execution
pool = gcp('nocreate');
if isempty(pool)
    parpool(Params.N_workers);
    pool = gcp('nocreate');
end
if pool.NumWorkers ~= Params.N_workers
    delete(gcp('nocreate'));
    parpool(Params.N_workers);
end

%% generate status files
for i = 1:Params.N_workers
    dlmwrite(sprintf('tmp/status_%02d', i), 1);
end

%% Generate adjust
Adjust = vertcat(repmat({'Kappa'},[Params.N_kappa_init 1]),repmat(vertcat(repmat({'Absorption'},[Params.N_absorption 1]), repmat({'Kappa'},[Params.N_kappa 1])), [1 1]));

%% Run Fit Kappa
while T_max_error > T_max_error_conv && T_mean_error > T_mean_error_conv && convergence_counter <= N_total
    
    %% get T for given Kappa & Absorption
    fprintf('%s -- iteration = %01d, ', datetime('now'), convergence_counter );
    fprintf(fid,'%s -- iteration = %01d ,', datetime('now'), convergence_counter);
    
    T_new = get_temp2D_parfor_gaussian(X_array, Y_array, Kappa, Absorption, Params);
    
    KAPPA{convergence_counter} = Kappa;
    TEMP{convergence_counter} = T_new;
    ABSORPTION{convergence_counter} = Absorption;
    
    %% get new error
    tmp = abs(T_exp - T_new);
    tmp(mask==0) = 0;
    tmp = tmp(:);
    tmp(tmp == 0) = [];
    
    T_Max_Error(convergence_counter) = max(tmp);
    T_Mean_Error(convergence_counter) = mean(tmp);
    
    %% display error
    fprintf('Max error = %1.2d, Mean error = %1.2d\n', T_Max_Error(convergence_counter), T_Mean_Error(convergence_counter));
    fprintf(fid,'Max error = %1.2d, Mean error = %1.2d\n', T_Max_Error(convergence_counter), T_Mean_Error(convergence_counter));
    
    %% adjust kappa
    if strcmp(Adjust{convergence_counter}, 'Kappa')
        
        Diff = convergence_slope_kappa * log10(abs(T_new - T_exp));
        Diff(mask==0) = 0;
        Diff(isinf(Diff)) = 0;
        Diff = sign(T_new - T_exp) .* 10.^Diff;
        
        Diff_av = zeros(size(Diff));
        for k = spatial_averaging+1:N_y-spatial_averaging
            for l = spatial_averaging+1:N_x-spatial_averaging
                if Diff(k,l)~=0
                    tmp = Diff(k-spatial_averaging:k+spatial_averaging,l-spatial_averaging:l+spatial_averaging);
                    tmp(tmp==0) = [];
                    Diff_av(k,l) = mean(tmp(:));
                end
            end
        end
        Diff_av(isnan(Diff_av)) = 0;
        
    else
        Diff = zeros(N_y, N_x);
        Diff_av = zeros(N_y, N_x);
    end
    
    Diff_av(mask == 0) = 0;
    Kappa_new = Kappa + Diff_av;
    Kappa_new(Kappa_new < Kappa_min) = Kappa_min;
    Kappa_new(Kappa_new > Kappa_max) = Kappa_max;
    
    DIFF_KAPPA{convergence_counter} = Diff;
    DIFF_KAPPA_SMOOTH{convergence_counter} = Diff_av;
    
    %% adjust absorption
    if strcmp(Adjust{convergence_counter}, 'Absorption')
        Diff = convergence_slope_absorption * (T_new - T_exp);
        Diff(mask==0) = 0;
        
        Diff_av = zeros(size(Diff));
        
        for k = spatial_averaging+1:N_y-spatial_averaging
            for l = spatial_averaging+1:N_x-spatial_averaging
                if Diff(k,l)~=0
                    tmp = Diff(k-spatial_averaging:k+spatial_averaging,l-spatial_averaging:l+spatial_averaging);
                    tmp(tmp==0) = [];
                    Diff_av(k,l) = mean(tmp(:));
                end
            end
        end
        Diff_av(isnan(Diff_av)) = 0;
        
    else
        Diff = zeros(N_y, N_x);
        Diff_av = zeros(N_y, N_x);
    end
    
    Diff_av(mask == 0) = 0;
    Absorption_new = Absorption - Diff_av;
    Absorption_new(Absorption_new < Absorption_min) = Absorption_min;
    
    DIFF_ABSORPTION{convergence_counter} = Diff;
    DIFF_ABSORPTION_SMOOTH{convergence_counter} = Diff_av;
    
    %% save data
    save('Data.mat')
    
    %% next iteration
    Absorption = Absorption_new;
    Kappa = Kappa_new;
    
    convergence_counter = convergence_counter + 1;
    
end


