%% clear
clear
clc
close all

%% Settings
Params.N_workers = 19;                       % number of parallel workers
Params.power = 4;                           % laser power [mW]
Params.temperature = 317;                   % bath temperature [K]
Params.name = '8um_S1_K12_post';            % sample name
Params.absorption_value = 0.027;            % absorption
Params.laser_spot_radius = 0.185;           % laser spot radius [um]
Params.membrane_radius = 4.0;               % membrane radius [um]
Params.membrane_mask = 3.5;                 % mask radius for ignoring edge effects [um]
spatial_averaging = 0;                      % spatial averaging when adjusting Kappa and absoprtion maps 
Params.power_correction = 1;                % optional power correction for accounting for power loss in system
Params.g = 28000000;                        % Interfacial [W/(K*m^2)]
Params.convection = 29000;                  % convection [W/(K*m^2)]

T_max_error_conv = 5;                       % optional convergence parameter for max temperature error [K]
T_mean_error_conv = 1;                      % optional convergence parameter for mean temperature error [K]
convergence_slope_kappa = 1.1;              % kappa adjustment (Beta_K)
convergence_slope_absorption = 0.0001;      % absorption adjustment (Beta_A)

Params.support_radius = 20.0;               % size of support (outer domain) [um]
Params.grid_resolution = 1;                 % comsol mesh resolution 1-10. 1- best, 10 worst
Params.display = 1;                         % set verbose on
Params.Timeout = 120;                       % timeout for connecting to comsol server
Params.laser_spot_radius_tot = 2.5;         % cut-off for laser Gaussian beam [um]

Kappa_min = 100;                            % lower bound for kappa [W/m/K]
Kappa_max = 4000;                           % upper bound for kappa [W/m/K]
Kappa_start = 1000;                         % initial guess kappa [W/m/K]
Params.Kappa_support = 500;                 % kappa value on support [W/m/K]

Absorption_min = 0.023;                     % lower bound for absorption 

Params.N_kappa = 50;                        % number of cycle to adjust kappa
Params.N_absorption = 100;                  % number of cycle to adjust absorption
Params.N_kappa_init = 100;                  % number of cycle to adjust kappa at start
N_total = 200;                              % Total number of iterations

%% add paths
addpath(genpath('functions'));
addpath(genpath('/opt/ud/comsol56/mli'));   % Matlab LiveLink path

%% disable mpiexec
distcomp.feature( 'LocalUseMpiexec', false )

%% clean up temporary files if available
system('rm tmp/ -r');
mkdir('tmp');

%% kill existing comsol servers
system('killall java');

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
    load(sprintf('Input_data/Formatted_data/%s/calculated_temperature_%1.2fmW_%1.0fK.mat', Params.name, Params.power, Params.temperature));
    X_array = X;
    Y_array = Y;
    
    fprintf('done\n');
    fprintf(fid, 'done\n');
    
    [N_y, N_x] = size(T_exp);
    
    %% load absorption
    fprintf('%s -- Getting absorption...', datetime('now'));
    fprintf(fid, '%s -- Getting absorption...', datetime('now'));
    
    Absorption = Params.absorption_value * ones(size(T_exp));
    
    %% double absorption on support
    [X,Y] = meshgrid(X_array, Y_array);
    Absorption( sqrt(X.^2 + Y.^2) > Params.membrane_radius) = 2 * Absorption_min;
    
    Absorption_init = Absorption;
    
    fprintf('done\n');
    fprintf(fid, 'done\n');
    
    %% Intialize mask for adjustement
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
        
    %% save data
    save('Data.mat')
    
    %% next iteration
    Absorption = Absorption_new;
    Kappa = Kappa_new;
    
    convergence_counter = convergence_counter + 1;
    
end