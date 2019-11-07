function [re_out,im_out,mag_out,phase_out,w_out] = astep_fotf_freqresp(G, w_low, w_high)   
    % ASTEP_FOTF_FREQRESP Calculate frequency response of fractional-order transfer function G
    % with an adaptive stepping technique. re_out and im_out are vectors of real and imaginary
    % values calculated at frequencies in vector w_out. mag_out and phase_out are the
    % magnitude and phase values of the frequency response at frequencies w_out.
    % 
    % [re_out,im_out,mag_out,phase_out,w_out] = ASTEP_FOTF_FREQRESP(G) - Frequency response of G
    % is evaluated and suitable lower and upper frequencies where the calculation should
    % terminate are automatically found.
    %
    % [re_out,im_out,mag_out,phase_out,w_out] = ASTEP_FOTF_FREQRESP(G, w_low, w_high) - Frequency
    % response of G is evaluated and the method stops when it moves outside the given frequency
    % bounds w_low and w_high.
    %
    % Note: Used correctly, the adaptive stepping technique should make sure that the calculated
    % frequency response is smooth in a logarithmic amplitude polar diagram. Frequency response
    % values are calculated sporadically so there are no linear or logarithmic structure to
    % the frequencies where the frequency response have been calculated.
    %
    % AUTHOR: Andreas Hanssen Moltumyr
    
    %% Choose function behavior based on input arguments
    switch nargin
        case 1
            w_start = 0;
            w_stop  = inf;
            w_interval_given = false;
        case 3
            w_start = w_low;
            w_stop  = w_high;
            w_interval_given = true;
        otherwise
            error('Number of input arguments not supported');
    end
    
    %% Configuration parameters (Tuneable)
    storage_size = 5000;
    
    % Precission tolerance of slope change between calculated values in degrees
    d     = 5;
    
    % Precission tolerance of step length magnitude compared to point magnitude
    y     = 12;
    
    % Multiplicative factor when changing step length up and down
    alpha = 2;
    
    % Tolerance relaxation value. Used when the method gets stuck [0,180/d - 1]
    dr    = 1;
    
    % Check extra frequency above w at stopping condition to be relatively
    % sure that the algorithm is not stopped to early. Scale the stopping
    % w with this value and check phase of that value.
    extra_stopping_test_w_scale = 100;
    
    % Change startpoint of w. The algorithm is started two times from this
    % point. One calculates from startpoint_w and towards inf, while
    % the other one calculates from startpoint_w and towards zero.
    startpoint_w                = 0.1;
    
    % If the algorithm only shows a circle or line (less than what you
    % would expect from the given fractional-order transfer function), then
    % try to change startpoint of w with "startpoint_w".
    
    %% Configuration parameters (Static)
    w_delta_start = 1e-6;
    deltaAngleTolerance_inDeg_atMagZero     = 5;
    deltaAngleTolerance_inDeg_atMagInf      = 5;
    magnitudeLength_at_Zero_constantPhase   = 30;
    magnitudeLength_at_Inf_constantPhase    = 30;
    constantEndValue_re_tolerance           = 1e-6;
    constantEndValue_im_tolerance           = 1e-6;
    constantEndValue_phase_tolerance_in_deg = 0.1;
    
    %% Extract coefficients and powers/orders from fotf object G
    [a,an,b,bn,L] = fotfdata(G);
    
    %% Calculate magnitude and phase at w = 0 and w = inf
    if w_interval_given == false
    [~, an_li] = min(an);
    [~, bn_li] = min(bn);
    [~, an_hi] = max(an);
    [~, bn_hi] = max(bn);
    
    % Calculate value at w = 0
    w_zero  = 0;
    magnitude_at_w_zero = abs(b(bn_li)/a(an_li))*w_zero^(bn(bn_li) - an(an_li));
    phase_at_w_zero     = wrapToPi(pi*(sign(b(bn_li)) ~= sign(a(an_li)))...
                                   + pi/2*(bn(bn_li) - an(an_li)));
    % Calculate value at w = inf
    w_inf  = inf;
    magnitude_at_w_inf = abs(b(bn_hi)/a(an_hi))*w_inf^(bn(bn_hi) - an(an_hi));
    phase_at_w_inf     = wrapToPi(pi*(sign(b(bn_hi)) ~= sign(a(an_hi)))...
                                  + pi/2*(bn(bn_hi) - an(an_hi)));
    
    % Calculate values used during stop checks
    switch magnitude_at_w_zero
        case 0
            lowerPhase_at_wZero = wrapToPi(phase_at_w_zero - (pi/180)*deltaAngleTolerance_inDeg_atMagZero);
            upperPhase_at_wZero = wrapToPi(phase_at_w_zero + (pi/180)*deltaAngleTolerance_inDeg_atMagZero);
            constPhase_stopLength_in_dB_at_wZero = -magnitudeLength_at_Zero_constantPhase;
        case inf
            lowerPhase_at_wZero = wrapToPi(phase_at_w_zero - (pi/180)*deltaAngleTolerance_inDeg_atMagInf);
            upperPhase_at_wZero = wrapToPi(phase_at_w_zero + (pi/180)*deltaAngleTolerance_inDeg_atMagInf);
            constPhase_stopLength_in_dB_at_wZero = magnitudeLength_at_Inf_constantPhase;
        otherwise
            if 0 < magnitude_at_w_zero && magnitude_at_w_zero < inf
                re_at_wZero = magnitude_at_w_zero*cos(phase_at_w_zero);
                im_at_wZero = magnitude_at_w_zero*sin(phase_at_w_zero);
            else
                error('NaN encountered');
            end
    end
    switch magnitude_at_w_inf
        case 0
            lowerPhase_at_wInf = wrapToPi(phase_at_w_inf - (pi/180)*deltaAngleTolerance_inDeg_atMagZero);
            upperPhase_at_wInf = wrapToPi(phase_at_w_inf + (pi/180)*deltaAngleTolerance_inDeg_atMagZero);
            constPhase_stopLength_in_dB_at_wInf = -magnitudeLength_at_Zero_constantPhase;
        case inf
            lowerPhase_at_wInf = wrapToPi(phase_at_w_inf - (pi/180)*deltaAngleTolerance_inDeg_atMagInf);
            upperPhase_at_wInf = wrapToPi(phase_at_w_inf + (pi/180)*deltaAngleTolerance_inDeg_atMagInf);
            constPhase_stopLength_in_dB_at_wInf = magnitudeLength_at_Inf_constantPhase;
        otherwise
            if 0 < magnitude_at_w_inf && magnitude_at_w_inf < inf
                re_at_wInf = magnitude_at_w_inf*cos(phase_at_w_inf);
                im_at_wInf = magnitude_at_w_inf*sin(phase_at_w_inf);
            else
                error('NaN encountered');
            end
    end
    constantEndValue_phase_tolerance_in_rad = (pi/180)*constantEndValue_phase_tolerance_in_deg;
    end
    
    %% Storage allocation
    REAL_LOW         = zeros(1,storage_size);
    IMAG_LOW         = zeros(1,storage_size);
    MAG_LOW          = zeros(1,storage_size);
    PHASE_LOW        = zeros(1,storage_size);
    OMEGA_LOW        = zeros(1,storage_size);
    DELTA_OMEGA_LOW  = zeros(1,storage_size);
    DELTA_PHASE_LOW  = zeros(1,storage_size);
    DELTA_MAGN_LOW   = zeros(1,storage_size);
    
    REAL_HIGH        = zeros(1,storage_size);
    IMAG_HIGH        = zeros(1,storage_size);
    MAG_HIGH         = zeros(1,storage_size);
    PHASE_HIGH       = zeros(1,storage_size);
    OMEGA_HIGH       = zeros(1,storage_size);
    DELTA_OMEGA_HIGH = zeros(1,storage_size);
    DELTA_PHASE_HIGH = zeros(1,storage_size);
    DELTA_MAGN_HIGH  = zeros(1,storage_size);
    
    %% Setting some initial values
    dd1     = (pi/180)*d;
    dd2     = (pi/180)*(360-d);
    d1      = dd1;
    d2      = dd2;
    
    j_stuck = 0;
    
    %% Calculate frequency response for first interval, w = [startpoint_w,inf]
    j                   = 1;
    w                   = startpoint_w;
    w_delta             = w_delta_start;
    stop_flag           = false;
    check_end_condition = false;
    sum_change_in_mag   = 0;
    
    while ( ((w_interval_given == false) && (stop_flag ~= true)) || ...
            ((w_interval_given == true)  && (w < w_stop))          )
        % 1.) Calculate frequency response value at w
        [re,im] = calc_freq_resp(w,a,an,b,bn,L);
        
        % The first two iterations of the loop is calculated without
        %  i. change in step length for w
        % ii. stopping condition check
        if j <= 2
            REAL_HIGH(j)        = re;
            IMAG_HIGH(j)        = im;
            MAG_HIGH(j)         = 10*log10(re^2 + im^2);
            PHASE_HIGH(j)       = atan2(im,re);
            OMEGA_HIGH(j)       = w;
            DELTA_OMEGA_HIGH(j) = w_delta;
            if j == 2
                DELTA_PHASE_HIGH(j) = atan2(im - IMAG_HIGH(j-1), re - REAL_HIGH(j-1));
                DELTA_MAGN_HIGH(j)  = 10.*log10((re - REAL_HIGH(j-1))^2 + (im - IMAG_HIGH(j-1))^2);
            end
            j = j + 1;
            w_prev = w;
            w = w_prev + w_delta;
        else
            % 2.) Check phase difference of the last line segment and the line
            % drawn between the previous point and the newly calculated
            % point. Also checks the increase in magnitude with
            % regards to the magnitude value at the current step.
            % If the checked values satisfy the given conditions then
            % increase the step size for w and save the frequency response
            % value. If the conditions are not meet, then decrease step
            % size for w and try a new step point by starting from the top
            % of the while loop.
            DELTA_PHASE_HIGH(j) = atan2(im - IMAG_HIGH(j-1), re - REAL_HIGH(j-1));
            DELTA_MAGN_HIGH(j)  = 10.*log10((re - REAL_HIGH(j-1))^2 + (im - IMAG_HIGH(j-1))^2);
            if ( (DELTA_MAGN_HIGH(j) > 10*log10(re^2 + im^2) - y) ...
                || (    abs(DELTA_PHASE_HIGH(j-1)-DELTA_PHASE_HIGH(j)) > d1 ...
                    &&  abs(DELTA_PHASE_HIGH(j-1)-DELTA_PHASE_HIGH(j)) < d2) ...
                || (w == inf) ) ...
            && (DELTA_OMEGA_HIGH(j-1) > eps(OMEGA_HIGH(j-1)))
                w_delta = w_delta/alpha;
                w       = w_prev + w_delta;
                
                % The iteration process can get stuck, this
                % piece of code prevents it, by slowly increasing the
                % acceptable range of phase difference values between the
                % subsequent line segments when a "stuck" condition is detected.
                if w_delta < (1e-6)*w
                    j_stuck = j;
                    d1 = (1+dr)*d1;
                    d2 = (1-dr)*d2;
                end
            else
                REAL_HIGH(j)        = re;
                IMAG_HIGH(j)        = im;
                MAG_HIGH(j)         = 10*log10(re^2 + im^2);
                PHASE_HIGH(j)       = atan2(im,re);
                OMEGA_HIGH(j)       = w;
                DELTA_OMEGA_HIGH(j) = w_delta;
                j = j + 1;
                check_end_condition = true;
            
                w_delta = w_delta*alpha;
                w_prev  = w;
                w = w_prev + w_delta;
                
                
                if j_stuck ~= 0 && j > j_stuck
                    % Iteration process managed to continue and
                    % d1 and d2 is reset
                    d1 = dd1;
                    d2 = dd2;
                    j_stuck = 0;
                end
            end
        end
        
        % 3.) Check if enough values have been calculated. If so, exit the loop.
        % If not, continue from the start of the while loop.
        if (check_end_condition == true) && (w_interval_given == false)
            if magnitude_at_w_inf == 0 || magnitude_at_w_inf == inf
                if angle_inside(lowerPhase_at_wInf,upperPhase_at_wInf,atan2(im,re))
                    sum_change_in_mag = sum_change_in_mag + (MAG_HIGH(j-1) - MAG_HIGH(j-2));
                    if sign(constPhase_stopLength_in_dB_at_wInf)*sum_change_in_mag >= abs(constPhase_stopLength_in_dB_at_wInf)
                        % Extra check for the existance of lower w value
                        % which gives a phase value different than the
                        % termination phase.
                        [re_stop_test,im_stop_test] = calc_freq_resp(extra_stopping_test_w_scale*w,a,an,b,bn,L);
                        if (re_stop_test < realmax && re_stop_test > -realmax) && (im_stop_test < realmax && im_stop_test > -realmax) && ~angle_inside(lowerPhase_at_wInf,upperPhase_at_wInf,atan2(im_stop_test,re_stop_test))
                            sum_change_in_mag = 0;
                        else
                            stop_flag = true;
                        end
                    end
                else
                    sum_change_in_mag = 0;
                end
            else
                if (   (abs(REAL_HIGH(j-1) - re_at_wInf) < constantEndValue_re_tolerance) ...
                    && (abs(IMAG_HIGH(j-1) - im_at_wInf) < constantEndValue_im_tolerance) ...
                    && (abs(PHASE_HIGH(j-1) - phase_at_w_inf) < constantEndValue_phase_tolerance_in_rad) )
                    stop_flag = true;
                end
            end
            check_end_condition = false;
        end
    end % while loop
    
    j_high_stop = j-1;
    
    %% Calc freq resp for w = [0,startpoint_w]
    j                   = 1;
    w                   = startpoint_w;
    w_delta             = w_delta_start;
    stop_flag           = false;
    check_end_condition = false;
    sum_change_in_mag   = 0;
    
    while ( ((w_interval_given == false) && (stop_flag ~= true)) || ...
            ((w_interval_given == true)  && (w > w_start))         )
        % 1.) Calculate frequency response value at w
        [re,im] = calc_freq_resp(w,a,an,b,bn,L);
        
        % The first two iterations of the loop is calculated without
        %  i. change in step length for w
        % ii. stopping condition check
        if j <= 2
            REAL_LOW(j)        = re;
            IMAG_LOW(j)        = im;
            MAG_LOW(j)         = 10*log10(re^2 + im^2);
            PHASE_LOW(j)       = atan2(im,re);
            OMEGA_LOW(j)       = w;
            DELTA_OMEGA_LOW(j) = w_delta;
            if j == 2
                DELTA_PHASE_LOW(j) = atan2(im - IMAG_LOW(j-1), re - REAL_LOW(j-1));
                DELTA_MAGN_LOW(j)  = 10.*log10((re - REAL_LOW(j-1))^2 + (im - IMAG_LOW(j-1))^2);
            end
            j = j + 1;
            w_prev = w;
            w = w_prev - w_delta;
        else
            % 2.) Check phase difference of the last line segment and the line
            % drawn between the previous point and the newly calculated
            % point. Also checks the increase in magnitude with
            % regards to the magnitude value at the current step.
            % If the checked values satisfy the given conditions then
            % increase the step size for w and save the frequency response
            % value. If the conditions are not meet, then decrease step
            % size for w and try a new step point by starting from the top
            % of the while loop.
            DELTA_PHASE_LOW(j) = atan2(im - IMAG_LOW(j-1), re - REAL_LOW(j-1));
            DELTA_MAGN_LOW(j)  = 10.*log10((re - REAL_LOW(j-1))^2 + (im - IMAG_LOW(j-1))^2);
            if (DELTA_MAGN_LOW(j) > 10*log10(re^2 + im^2) - y ...
                || (   abs(DELTA_PHASE_LOW(j-1)-DELTA_PHASE_LOW(j)) > d1 ...
                    && abs(DELTA_PHASE_LOW(j-1)-DELTA_PHASE_LOW(j)) < d2) ) ...
            && (DELTA_OMEGA_HIGH(j-1) > eps(OMEGA_HIGH(j-1)))
                w_delta = w_delta/alpha;
                w       = w_prev - w_delta;
                
                % The iteration process can get stuck, this
                % piece of code prevents it, by slowly increasing the
                % acceptable range of phase difference values between the
                % subsequent line segments when a "stuck" condition is detected.
                if w_delta < (1e-6)*w
                    j_stuck = j;
                    d1 = (1+dr)*d1;
                    d2 = (1-dr)*d2;
                end
            else
                REAL_LOW(j)        = re;
                IMAG_LOW(j)        = im;
                MAG_LOW(j)         = 10*log10(re^2 + im^2);
                PHASE_LOW(j)       = atan2(im,re);
                OMEGA_LOW(j)       = w;
                DELTA_OMEGA_LOW(j) = w_delta;
                j = j + 1;
                check_end_condition = true;
            
                w_delta = w_delta*alpha;
                w_prev  = w;
                w = w_prev - w_delta;
                
                % 
                if j_stuck ~= 0 && j > j_stuck
                    % Iteration process managed to continue and
                    % d1 and d2 is reset
                    d1 = dd1;
                    d2 = dd2;
                    j_stuck = 0;
                end
            end
        end
        
        while w < 0
            w_delta = w_delta/alpha;
            w       = w_prev - w_delta;
        end
        
        % 3.) Check if enough values have been calculated. If so, exit the loop.
        % If not, continue from the start of the while loop.
        if (check_end_condition == true) && (w_interval_given == false)
            if magnitude_at_w_zero == 0 || magnitude_at_w_zero == inf
                if angle_inside(lowerPhase_at_wZero,upperPhase_at_wZero,atan2(im,re))
                    sum_change_in_mag = sum_change_in_mag + (MAG_LOW(j-1) - MAG_LOW(j-2));
                    if sign(constPhase_stopLength_in_dB_at_wZero)*sum_change_in_mag >= abs(constPhase_stopLength_in_dB_at_wZero)
                        % Extra check for the existance of lower w value
                        % which gives a phase different than the
                        % termination phase
                        [re_stop_test,im_stop_test] = calc_freq_resp(1/extra_stopping_test_w_scale*w,a,an,b,bn,L);
                        if (re_stop_test < realmax && re_stop_test > -realmax) && (im_stop_test < realmax && im_stop_test > -realmax) && ~angle_inside(lowerPhase_at_wZero,upperPhase_at_wZero,atan2(im_stop_test,re_stop_test))
                            sum_change_in_mag = 0;
                        else
                            stop_flag = true;
                        end
                    end
                else
                    sum_change_in_mag = 0;
                end
            else
                if (    (abs(REAL_LOW(j-1) - re_at_wZero) < constantEndValue_re_tolerance) ...
                     && (abs(IMAG_LOW(j-1) - im_at_wZero) < constantEndValue_im_tolerance) ...
                     && (abs(PHASE_LOW(j-1) - phase_at_w_zero) < constantEndValue_phase_tolerance_in_rad) )
                    stop_flag = true;
                end
            end
            check_end_condition = false;
        end
    end % while loop
    
    j_low_stop = j-1;
    
    %% Add more points to include important approching zero or inf dynamics
    if w_interval_given == false
        min_mag = min([MAG_LOW(j_low_stop:-1:1), MAG_HIGH(1:j_high_stop)]);
        max_mag = max([MAG_LOW(j_low_stop:-1:1), MAG_HIGH(1:j_high_stop)]);
        
        % Calculate the wished lower value for zero approching
        if min_mag == MAG_LOW(j_low_stop) || min_mag == MAG_HIGH(j_high_stop)
            wished_lower_mag_dB = min_mag;
        else
            wished_lower_mag_dB = min_mag - 10;            
        end
        wished_lower_mag = 10^(wished_lower_mag_dB/20);
        
        % Calculate the wished upper value for inf approaching
        if max_mag == MAG_LOW(j_low_stop) || max_mag == MAG_HIGH(j_high_stop)
            wished_upper_mag_dB = max_mag;
        else
            wished_upper_mag_dB = max_mag + 10;            
        end
        wished_upper_mag = 10^(wished_upper_mag_dB/20);
        
        % Add extra point at w = 0
        switch magnitude_at_w_zero
            case 0
                if min_mag ~= MAG_LOW(j_low_stop)
                    j = j_low_stop + 1;
                    REAL_LOW(j)  = wished_lower_mag*cos(phase_at_w_zero);
                    IMAG_LOW(j)  = wished_lower_mag*sin(phase_at_w_zero);
                    MAG_LOW(j)   = wished_lower_mag_dB;
                    PHASE_LOW(j) = phase_at_w_zero;
                    OMEGA_LOW(j) = (abs(a(an_li)/b(bn_li))*wished_lower_mag)^(1/(bn(bn_li) - an(an_li)));
                    j_low_stop = j;
                end
            case inf
                if max_mag ~= MAG_LOW(j_low_stop)
                    j = j_low_stop + 1;
                    REAL_LOW(j)  = wished_upper_mag*cos(phase_at_w_zero);
                    IMAG_LOW(j)  = wished_upper_mag*sin(phase_at_w_zero);
                    MAG_LOW(j)   = wished_upper_mag_dB;
                    PHASE_LOW(j) = phase_at_w_zero;
                    OMEGA_LOW(j) = (abs(a(an_li)/b(bn_li))*wished_upper_mag)^(1/(bn(bn_li) - an(an_li)));
                    j_low_stop = j;
                end
            otherwise
                % Do Nothing, the frequency response is good enough as it is
        end
        
        % Add extra point at w = inf
        switch magnitude_at_w_inf
            case 0
                if min_mag ~= MAG_HIGH(j_high_stop)
                    j = j_high_stop + 1;
                    REAL_HIGH(j)  = wished_lower_mag*cos(phase_at_w_inf);
                    IMAG_HIGH(j)  = wished_lower_mag*sin(phase_at_w_inf);
                    MAG_HIGH(j)   = wished_lower_mag_dB;
                    PHASE_HIGH(j) = phase_at_w_inf;
                    OMEGA_HIGH(j) = (abs(a(an_hi)/b(bn_hi))*wished_lower_mag)^(1/(bn(bn_hi) - an(an_hi)));
                    j_high_stop = j;
                end
            case inf
                if max_mag ~= MAG_HIGH(j_high_stop)
                    j = j_high_stop + 1;
                    REAL_HIGH(j)  = wished_upper_mag*cos(phase_at_w_inf);
                    IMAG_HIGH(j)  = wished_upper_mag*sin(phase_at_w_inf);
                    MAG_HIGH(j)   = wished_upper_mag_dB;
                    PHASE_HIGH(j) = phase_at_w_inf;
                    OMEGA_HIGH(j) = (abs(a(an_hi)/b(bn_hi))*wished_upper_mag)^(1/(bn(bn_hi) - an(an_hi)));
                    j_high_stop = j;
                end
            otherwise
                % Do Nothing, the frequency response is good enough as it is
        end
    end
    
    %% Add w = 0 and w = inf points
    if w_interval_given == false
        % Add point at w = 0
        j = j_low_stop + 1;
        REAL_LOW(j)  = magnitude_at_w_zero*cos(phase_at_w_zero);
        IMAG_LOW(j)  = magnitude_at_w_zero*sin(phase_at_w_zero);
        MAG_LOW(j)   = 20*log10(magnitude_at_w_zero);
        PHASE_LOW(j) = phase_at_w_zero;
        OMEGA_LOW(j) = 0;
        j_low_stop = j;
        
        % Add point at w = inf
        j = j_high_stop + 1;
        REAL_HIGH(j)  = magnitude_at_w_inf*cos(phase_at_w_inf);
        IMAG_HIGH(j)  = magnitude_at_w_inf*sin(phase_at_w_inf);
        MAG_HIGH(j)   = 20*log10(magnitude_at_w_inf);
        PHASE_HIGH(j) = phase_at_w_inf;
        OMEGA_HIGH(j) = inf;
        j_high_stop = j;
    end
    
    %% Return the calculated frequency response
    re_out    =           [ REAL_LOW(j_low_stop:-1:1),  REAL_HIGH(1:j_high_stop)];
    im_out    =           [ IMAG_LOW(j_low_stop:-1:1),  IMAG_HIGH(1:j_high_stop)];
    mag_out   =           [  MAG_LOW(j_low_stop:-1:1),   MAG_HIGH(1:j_high_stop)];
    phase_out = (180/pi).*[PHASE_LOW(j_low_stop:-1:1), PHASE_HIGH(1:j_high_stop)];
    w_out     =           [OMEGA_LOW(j_low_stop:-1:1), OMEGA_HIGH(1:j_high_stop)];
    
end % astep_fotf_freqresp()


%% Helper functions
function [re, im] = calc_freq_resp(w, a, na, b, nb, L)
    freq_resp_num = b*((1i*w).^nb.');
    freq_resp_den = a*((1i*w).^na.');
    freq_resp = freq_resp_num/freq_resp_den;
    if L > 0
        freq_resp = freq_resp.*exp(-L*1i*w);
    end
    
    re = real(freq_resp);
    im = imag(freq_resp);
end

function isInside = angle_inside(lower_angle,upper_angle,test_angle)
    upper_angle = upper_angle - lower_angle;
    if upper_angle < 0
        upper_angle = upper_angle + 360;
    end
    test_angle = test_angle - lower_angle;
    if test_angle < 0
        test_angle = test_angle + 360;
    end
    isInside = (test_angle < upper_angle);
end
