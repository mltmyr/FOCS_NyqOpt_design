%% FILE: regulator_optimization_con.m
%% AUTHOR: Andreas Hanssen Moltumyr

function [c, ceq] = regulator_optimization_con(x, G, reg_type, Np)
    %% Construct open-loop system transfer function
    L = construct_system_fotfs(G, reg_type, x);
    
    %% Check regulator open loop stability so to adjust needed encirclement number Np
    % N_p: Needed encirclements in nyquist diagram for closed-loop stability
    switch reg_type
        case 'FO-PPF_2'
            zeta_d = x(3);
            beta_d = x(5);
            if zeta_d > -cos(pi/2*beta_d)
                N_p = Np;
            else
                N_p = Np + 2;
            end
        case 'FO-PPF_3'
            zeta_d = x(4);
            beta_d = x(6);
            if zeta_d > -cos(pi/2*beta_d)
                N_p = Np;
            else
                N_p = Np + 2;
            end
        otherwise
            N_p = Np;
    end
    
    %% Evaluate Stability (Nyquist stability criterion approach)
    w_low  = 10^-10;
    w_high = 10^5;
    [isStable,~,GMs_dB,~,~,~] = assess_stability(L,N_p,{w_low,w_high});
    
    %% Calculate Equality Conditions (Stability)
    c_eq1 = 100*(1 - isStable);
    
    ceq = [c_eq1];
    
    %% Calculate Inequality Conditions (Stability)
    gain_margin_db = 6;
    
    if size(GMs_dB,2) ~= 0
        c_ineq1 = gain_margin_db - GMs_dB(2,1);
        c_ineq2 = GMs_dB(1,1) + gain_margin_db;
    else
        c_ineq1 = -realmax;
        c_ineq2 =  realmax;
    end
        
    if c_ineq1 == inf
        c_ineq1 = realmax;
    elseif c_ineq1 == -inf
        c_ineq1 = -realmax;
    end
    if c_ineq2 == inf
        c_ineq2 = realmax;
    elseif c_ineq2 == -inf
        c_ineq2 = -realmax;
    end
    
    c = [c_ineq1, c_ineq2];
end
