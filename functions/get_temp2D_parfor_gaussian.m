%% clear
function T = get_temp2D_parfor_gaussian(X_pos, Y_pos, Kappa, Absorption, Params)

import com.comsol.model.*
import com.comsol.model.util.*

N_X = length(X_pos);
N_Y = length(Y_pos);

%% write kappa to file
[X,Y] = meshgrid(X_pos, Y_pos);
X = X(:);
Y = Y(:);
Z = Kappa(:);
Mask = Params.mask(:);

fid = fopen([pwd '/tmp/kappa.txt'],'w');
for i = 1:length(Z)
    fprintf(fid,'%1.3f %1.3f %1.3f\n',X(i), Y(i), Z(i));
end
fclose(fid);

%% write absorption to file
Z = Params.power * 1e-3 * Absorption(:) * Params.power_correction;
fid = fopen([pwd '/tmp/absorption.txt'],'w');
for i = 1:length(Z)
    fprintf(fid,'%1.3f %1.3f %1.3e\n',X(i), Y(i), Z(i));
end
fclose(fid);

%% remove points outside mask
T = zeros(N_Y * N_X, 1);
N = find(Mask==1);
N_calcs = length(N);
X = X(N);
Y = Y(N);

Temp = zeros(N_calcs, 1);

%% solve heat equation for all laser spot positions
Display = Params.display;

%% par for loop
parfor index1 = 1:N_calcs
    
    run = true;
    
    %% run single calculation
    while run
        
        %% get worker id
        t = getCurrentTask();
        
        %% logfile
        fid = fopen(sprintf('tmp/logfile_%02d', t.ID),'a');
        
        %% start server
        if dlmread(sprintf('tmp/status_%02d', t.ID)) == 1
            
            done = 0;
            pause(t.ID);

            while done == 0
                fprintf(fid, '%s -- Starting comsol server... \n', datetime('now'));
                fprintf('%s -- Starting comsol server... \n', datetime('now'));
                system( sprintf('comsol server -np 1 -portfile tmp/port_%02d & echo $! > tmp/pid_%02d', t.ID, t.ID));
                
                counter = 1;
                
                % wait for port file
                while counter < 60 && done == 0
                    if exist(sprintf('tmp/port_%02d', t.ID),'file') == 2
                        done = 1;
                        dlmwrite(sprintf('tmp/status_%02d', t.ID), 2);
                        Port = dlmread(sprintf('tmp/port_%02d', t.ID));
                        fprintf(fid, 'port %02d \n', Port );
                        fprintf('port %02d \n', Port );
                    else
                        pause(1)
                        counter = counter + 1;
                    end
                end
            end
        end
        
        %% connect to server
        if dlmread(sprintf('tmp/status_%02d', t.ID)) == 2
            
            try
                Port = dlmread(sprintf('tmp/port_%02d', t.ID));
                mphstart(Port);
                dlmwrite(sprintf('tmp/status_%02d', t.ID), 3);
                fprintf(fid, '%s -- Connected to server\n', datetime('now'));
                fprintf('%s -- Connected to server\n', datetime('now'));
            catch
                fprintf(fid, '%s -- Cannot connect to server\n', datetime('now'));
                fprintf('%s -- Cannot connect to server\n', datetime('now'));
                dlmwrite(sprintf('tmp/status_%02d', t.ID), 1);
                delete(sprintf('tmp/port_%02d', t.ID));
                system(sprintf('kill %01d', dlmread(sprintf('tmp/pid_%02d', t.ID))));
                delete(sprintf('tmp/pid_%02d', t.ID));
            end

        end
        
        %% run calculations
        if dlmread(sprintf('tmp/status_%02d', t.ID)) == 3
            
            if  Display
                fprintf(fid, sprintf('%s - %04d/%04d\n', datetime('now'), index1, N_calcs));
                fprintf(sprintf('%s - %04d/%04d\n', datetime('now'), index1, N_calcs));
            end
            
            pos_x = X(index1);
            pos_y = Y(index1);
            
            % record start time
            fid2 = fopen(sprintf('tmp/worker_time_%02d', t.ID), 'w');
            fprintf(fid2,'%1.0f\n', posixtime(datetime));
            fclose(fid2);
            
            % run calculation
            try
                Temp(index1) = comsol_model_gaussian(pos_x, pos_y, Params);
                run = false;
            catch ME
                display(ME)
                fprintf(fid, '%s -- Calculation %01d failed\n', datetime('now'), index1);
                fprintf('%s -- Calculation %01d failed\n', datetime('now'), index1);
                dlmwrite(sprintf('tmp/status_%02d', t.ID), 2);
            end
        end
        
        fclose(fid);
        
    end
end

%% assemble data
T(Mask==0) = Params.temperature;
T(Mask==1) = Temp;

T = reshape(T, [N_Y N_X]);

end
