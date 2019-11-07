function y = regulator_optimization_opt(x, G, reg_type)
    % REGULATOR_OPTIMIZATION_OPT Function for defining and calculating
    % the objectivity used in the optimization problem.
    % 
    % y = regulator_optimization_opt(x, G, reg_type) - x are the control parameters to optimize,
    % G is the plant and reg_type is a regulator identifier that is used by the 
    % construct_system fotfs(...) function in order to build the correct regulator from the control
    % parameters x. y are the returned score used by the optimization algorithm to assess the solution x.
    %
    % If a custom objectivity function is needed, it can be defined below.
    %
    % AUTHOR: Andreas Hanssen Moltumyr

    %% Construct closed-loop system transfer function
    [~, T, So, ~] = construct_system_fotfs(G, reg_type, x);
    
    %% Calculate frequency responses
    w_low  = -4;
    w_high = 4;
    N      = 1000;
    W = logspace(w_low, w_high, N);
    
    Z_T   = freqresp(1i*W, T);
    T_MAG = 20.*log10(abs(Z_T));

    Z_So   = freqresp(1i*W, So);
    So_MAG = 20.*log10(abs(Z_So));

    %% Objectivity choice
    objectivity_to_use = 4;
    
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
        case 3
            obj_weight_1 = 4;
            obj_weight_2 = 1;
            obj_weight_3 = 4;
            bw_defined_at_mag_dB = -6;
            
            [bw, idx] = find_bandwidth(T_MAG, W, bw_defined_at_mag_dB);
            T_vec_below_bw = abs(squeeze(T_MAG(1:idx)));
            T_vec_above_bw = squeeze(T_MAG(idx:end));
            obj1 = norm(T_vec_below_bw, 2)/idx;
            obj2 = -bw;
            obj3 = sum(T_vec_above_bw)/(length(squeeze(T_MAG))-idx);
            
            obj_weight_sum = obj_weight_1 + obj_weight_2 + obj_weight_3;
            obj_weight_1 = obj_weight_1/obj_weight_sum;
            obj_weight_2 = obj_weight_2/obj_weight_sum;
            obj_weight_3 = obj_weight_3/obj_weight_sum;
            
            objective = obj_weight_1*obj1 + obj_weight_2*obj2 + obj_weight_3*obj3;
        case 4
            obj_w1 = 50;
            obj_w2 = 1;
            obj_w3 = 50;
            bw_defined_at_mag_dB = -6;
            
            [bwT, idxT]  = find_bw_T(T_MAG, W, bw_defined_at_mag_dB);
            [bwSo,idxSo] = find_bw_So(So_MAG, W, bw_defined_at_mag_dB);
            T_vec_below_bw  = abs(squeeze(T_MAG(1:idxT)));
            So_vec_above_bw = abs(squeeze(So_MAG(idxSo:length(So_MAG))));
            
            obj1 = norm(T_vec_below_bw, 2)/idxT;
            obj2 = -min([bwT,bwSo]);
            obj3 = norm(So_vec_above_bw, 2)/(length(So_MAG)-idxSo);
            
            obj_weight_sum = obj_w1 + obj_w2 + obj_w3;
            obj_w1 = obj_w1/obj_weight_sum;
            obj_w2 = obj_w2/obj_weight_sum;
            obj_w3 = obj_w3/obj_weight_sum;
            
            objective = obj_w1*obj1 + obj_w2*obj2 + obj_w3*obj3;
        case 5
            obj_w1 = 4;
            obj_w2 = 1;
            obj_w3 = 4;
            bw_defined_at_mag_dB = -6;
            
            [bwT,~]  = find_bw_T(T_MAG, W, bw_defined_at_mag_dB);
            [bwSo,~] = find_bw_So(So_MAG, W, bw_defined_at_mag_dB);
            T_vec  = squeeze(T_MAG);
            So_vec = squeeze(So_MAG);
            
            obj1 = max(T_vec);
            obj2 = -min([bwT,bwSo]);
            obj3 = max(So_vec);
            
            obj_weight_sum = obj_w1 + obj_w2 + obj_w3;
            obj_w1 = obj_w1/obj_weight_sum;
            obj_w2 = obj_w2/obj_weight_sum;
            obj_w3 = obj_w3/obj_weight_sum;
            
            objective = obj_w1*obj1 + obj_w2*obj2 + obj_w3*obj3;
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

function [bw, idx] = find_bw_T(MAG, W, mag_at_bw)
    bw = W(1); idx = 1;
    for j = 1:length(MAG)-1
        if (MAG(j) >= mag_at_bw)
            bw  = W(j);
            idx = j;
        else 
            break;
        end
    end
end

function [bw, idx] = find_bw_So(MAG, W, mag_at_bw)
    bw = W(1); idx = 1;
    for j = 1:length(MAG)-1
        if (MAG(j) <= mag_at_bw)
            bw  = W(j);
            idx = j;
        else 
            break;
        end
    end
end