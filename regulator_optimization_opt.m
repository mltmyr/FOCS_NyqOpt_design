%% FILE: regulator_optimization_opt.m
%% AUTHOR: Andreas Hanssen Moltumyr

function y = regulator_optimization_opt(x, G, reg_type)
    %% Construct closed-loop system transfer function
    [~, T, ~] = construct_system_fotfs(G, reg_type, x);
    
    %% Calculate frequency responses
    w_low  = -1;
    w_high = 4;
    N      = 1000;
    W = logspace(w_low, w_high, N);
    
    Z_T   = freqresp(1i*W, T);
    T_MAG = 20.*log10(abs(Z_T));

    %% Objectivity choice
    objectivity_to_use = 2;
    
    %% Objectivity function calculation
    switch objectivity_to_use
        case 1
            T_vec = 1-abs(squeeze(T_MAG));
            objective = norm(T_vec, 2);
        case 2
            obj_weight = 0.1; % [0-1] Weight hight bandwidth against 
                              % flat closed-loop response up to bandwidth.
            bw_defined_at_mag_dB = -6;
            
            [bw, idx] = find_bandwidth(T_MAG, W, bw_defined_at_mag_dB);
            T_vec = 1-abs(squeeze(T_MAG(1:idx)));
            obj1 = norm(T_vec/idx, 2);
            obj2 = -bw;
            objective = (1-obj_weight)*obj1 + obj_weight*obj2;
        otherwise
            error(['Objectivity function nr. ', ...
                num2str(objectivity_to_use), 'does not exist']);
    end
    
    %% Return objective value
    y = objective;
end

function [bandwidth, index] = find_bandwidth(T_mag_dB, W, magnitude_dB_at_bandwidth)
    mag_dB_at_bw = magnitude_dB_at_bandwidth;
    
    bw  = W(1); idx = 1;
    for j = 1:length(T_mag_dB)-1
        if (T_mag_dB(j) <= abs(mag_dB_at_bw) && T_mag_dB(j) >= -abs(mag_dB_at_bw))
            bw  = W(j);
            idx = j;
        else
            break;
        end
    end
    
    bandwidth = bw;
    index     = idx;
end
