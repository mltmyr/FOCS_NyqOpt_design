%% FILE: evaluate_regulator.m
%% AUTHOR: Andreas Hanssen Moltumyr

%Function that takes in a path to a run folder and adds evaluation criteria to the data.mat files
function evaluate_regulator(results_path, regulator_name, run_num_or_name)
path_to_data = {};
switch nargin
    case 1
        path_to_data{1} = results_path;
        eval_single_reg(path_to_data);
    case 3
        if isa(run_num_or_name, 'char')
            path_to_data{1} = fullfile(results_path, regulator_name, run_num_or_name, 'data.mat');
        elseif isa(run_num_or_name, 'cell')
            for j = 1:length(run_num_or_name)
                path_to_data{j} = fullfile(results_path, regulator_name, run_num_or_name{j}, 'data.mat');
            end
        elseif isa(run_num_or_name, 'double') || isa(run_num_or_name, 'single')
            for j = 1:length(run_num_or_name)
                path_to_data{j} = fullfile(results_path, regulator_name, ['Run', num2str(run_num_or_name(j))], 'data.mat');
            end
        else
            error('Error in input arguments!');
        end
        
        for j = 1:length(path_to_data)
            %try
                eval_single_reg(path_to_data{j});
            %catch
                %warning(['Can not find/access/open ', path_to_data{j}]);
            %end
        end
    otherwise
        error('Error in input arguments!');
end
    
end

%% Worker function
function eval_single_reg(path_to_data)

    data = load(path_to_data);
    
    switch data.regulator_type
        case {'PPF', 'FO-PPF_1', 'FO-PPF_2', 'FO-PPF_3'}
            C_t = data.misc(1);
            C_d = data.misc(2);
            C = C_t*C_d;
            
            So = data.So;
            T = data.G*C*So;
        case {'PID', 'FO-PID'}
            C =data.misc(1);
            
            So = data.So;
            T = data.G*C*So;
        otherwise
            error('Regulator type not supported!');
    end
    
    %% Calculate Bandwidth
    mag_db_at_bandwidth = -3;
    
    w_low  = -5;
    w_high = 5;
    w_N    = 10000;
    
    w_array = logspace(w_low,w_high,w_N);
    T_frd = bode(T,w_array);
    MAG_T_dB = 20.*log10(abs(squeeze(T_frd.ResponseData)));
    W_T = w_array;
    
    j = 0;
    found_bandwidth_point = false;
    while ~found_bandwidth_point && j <= length(MAG_T_dB)-1
        j = j + 1;
        if MAG_T_dB(j) > mag_db_at_bandwidth && MAG_T_dB(j+1) < mag_db_at_bandwidth
            found_bandwidth_point = true;
        end
    end
    T_bandwidth = W_T(j);
    
    data.bandwidth = T_bandwidth;
    
    %% Calculate ||N(s)||_{\infty} max
    So_frd = bode(So,w_array); %[MAG_So,~]
    MAG_So_dB = 20.*log10(abs(squeeze(So_frd.ResponseData)));
    W_So = w_array;
    [~, max_index] = max(MAG_So_dB);
    
    w_So_zoom_scale = 0.5;
    w_So_low  = log10(W_So(max_index))+log10(w_So_zoom_scale);
    w_So_high = log10(W_So(max_index))-log10(w_So_zoom_scale);
    w_So_N    = 10000;
    w_array_zoomed = logspace(w_So_low,w_So_high,w_So_N);
    So_frd = bode(So,w_array_zoomed);
    MAG_So_dB = 20.*log10(abs(squeeze(So_frd.ResponseData)));
    So_max_val = max(MAG_So_dB);
    
    data.sensitivity_max_dB = So_max_val;
    
    %% Reference-Response Error of Step response
    simStepAmplitude         = data.simStepAmplitude;
    simStepSimulationTime    = data.simStepSimulationTime;
    
    simFractionalFilterOrder      = data.simFractionalFilterOrder;
    simFractionalOrderFrequencies = data.simFractionalOrderFrequencies;
    
    
    sim_options = simset('SrcWorkspace','current');
    
    G = data.G;
    
    switch data.regulator_type
        case {'PPF', 'FO-PPF_1', 'FO-PPF_2', 'FO-PPF_3'}
            C_t = data.misc(1);
            C_d = data.misc(2);
            sim('step_response_simulation_foppf.slx',[],sim_options);
        case {'PID', 'FO-PID'}
            C = data.misc(1);
            sim('step_response_simulation_fopid.slx',[],sim_options);
        otherwise
                error('Choose a supported regulator or implement support for new regulator');
    end

    t_sim  = sim_data.Time;
    y_sim  = sim_data.Data(:,1);
    y0_sim = sim_data.Data(:,2);
    
    for j = 1:length(t_sim)
        if t_sim(j) >= 0.01*simStepSimulationTime
            idx = j;
            break;
        end
    end
    
    e_y = y0_sim(1:end-1) - y_sim(1:end-1);
    dt  = t_sim(2:end) - t_sim(1:end-1);
    
    % Integral of Squared Error (ISE)
    data.ISE  = sum(e_y(idx:end).^2.*dt(idx:end));
    
    % Integral of Absolute Error (IAE)
    data.IAE  = sum(abs(e_y(idx:end)).*dt(idx:end));
    
    % Integral of Time-weighted Absolute Error (ITAE)
    data.ITAE = sum(t_sim(idx+1:end).*abs(e_y(idx:end)).*dt(idx:end));
    
    %% Save evaluation criteria
    save(path_to_data, '-struct', 'data');
end
