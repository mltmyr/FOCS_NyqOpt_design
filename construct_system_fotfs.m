function [L, T, So, misc] = construct_system_fotfs(G_plant, reg_type, reg_param)
    % CONSTRUCT_SYSTEM_FOTFS Function where regulator parameters reg_param are
    % converted to a fotf object (fractional-order transfer function) based on
    % the identifier reg_type.
    %
    % [L, T, So, misc] = construct_system_fotfs(G_plant, reg_type, reg_param)
    %
    % Open-loop and Closed-loop transfer functions L(s),
    % T(s) (Complementary sensitivity function) and So(s) (Sensitivity function)
    % are also calculated and returned. Possible other interesting transfer functions
    % and data can be returned as MISC.
    %
    % If other controllers than the ones defined here should be used and optimized, 
    % necessary conversions and calculations should be added here.
    %
    % AUTHOR: Andreas Hanssen Moltumyr

    s = fotf('s');
    
    %% Input guard
    if isa(G_plant, 'tf')
        G_plant = fotf(G_plant);
    elseif ~isa(G_plant, 'fotf')
        error(['Unsupported type ', class(G_plant), ' of argument G']);
    end
    
    %% Extract Regulator Parameters from x based on regulator
    x = reg_param;
    
    switch reg_type
        case 'PPF'
            k_t = x(1);
            k_d = x(2); zeta_d  = x(3); omega_d = x(4);
        case 'FO-PPF_1'
            k_t = x(1); alpha_t = x(2);
            k_d = x(3); zeta_d  = x(4); omega_d = x(5);
        case 'FO-PPF_2'
            k_t = x(1);
            k_d = x(2); zeta_d  = x(3); omega_d = x(4); beta_d = x(5);
        case 'FO-PPF_3'
            k_t = x(1); alpha_t = x(2);
            k_d = x(3); zeta_d  = x(4); omega_d = x(5); beta_d = x(6);
        case 'FO-PID'
            k_p = x(1); k_i  = x(2); mu_i  = x(4);
            k_d = x(3); mu_d = x(5); tau_f = x(6);
        case 'PID'
            k_p = x(1); k_i = x(2); k_d = x(3); tau_f = x(4);
        otherwise
            error('Choose a supported regulator or implement support for new regulator');
    end
    
    %% Construct Regulators and Define System-Interconnection
    % G - Plant transfer function
    % C - Forward Controller transfer function
    % F - Backward controller transfer function
    % L - Open-loop transfer function
    switch reg_type
        case 'PPF'
            C_t = -k_t/s;
            C_d = -(k_d*omega_d^2)/(s^2 + 2*zeta_d*omega_d*s + omega_d^2);
            
            C = C_t*C_d;
            F = 1 + C_t^(-1);
            L = F*G_plant*C;
            
            misc = [C_t, C_d];
        case 'FO-PPF_1'
            C_t = -k_t/(s^alpha_t);
            C_d = -(k_d*omega_d^2)/(s^2 + 2*zeta_d*omega_d*s + omega_d^2);
            
            C = C_t*C_d;
            F = 1 + C_t^(-1);
            L = F*G_plant*C;
            
            misc = [C_t, C_d];
        case 'FO-PPF_2'
            C_t = -k_t/s;
            C_d = -(k_d*omega_d^2)/(s^(2*beta_d) + 2*zeta_d*omega_d*s^(beta_d) + omega_d^2);
            
            C = C_t*C_d;
            F = 1 + C_t^(-1);
            L = F*G_plant*C;
            
            misc = [C_t, C_d];
        case 'FO-PPF_3'
            C_t = -k_t/(s^alpha_t);
            C_d = -(k_d*omega_d^2)/(s^(2*beta_d) + 2*zeta_d*omega_d*s^(beta_d) + omega_d^2);
            
            C = C_t*C_d;
            F = 1 + C_t^(-1);
            L = F*G_plant*C;
            
            misc = [C_t, C_d];
        case 'FO-PID'
            C = k_p + k_i/(s^mu_i) + (k_d*s^(mu_d))/(1+tau_f*s^mu_d);
            
            L = G_plant*C;
            misc = [C];
        case 'PID'
            C = k_p + k_i/s + (k_d*s)/(1+tau_f*s);
            
            L = G_plant*C;
            misc = [C];
        otherwise
            error(['reg_type = ', reg_type, ' is not supported! If it should be supported then add it!']);
    end
    
    %% Calculate Sensitivity and Complementary Sensitivity functions
    switch nargout
        case 1
            % Allready finished so do nothing
        case {2,3,4}
            So = (1+L)^(-1);
            T  = G_plant*C*So;
        otherwise
            error('Number of output arguments not supported')
    end
        
end    
