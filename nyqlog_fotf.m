function nyqlog_fotf(G,w1,w2)
    % NYQLOG_FOTF Function for plotting a logarithmic amplitude polar diagram
    % of a fractional-order transfer function (fotf object) or
    % integer-order transfer function (tf object).
    %
    % nyqlog_fotf(G) G is the transfer function to plot. Uses the adaptive stepping
    % function astep_fotf_freqresp(...) to calculate frequency response.
    %
    % nyqlog_fotf(G,w1) A vector of frequencies w1 can be specified if the diagram
    % should be plotted for a set of specific frequencies. Uses the fotf/bode(...) function
    % from the fotf toolbox to calculate frequency response.
    % 
    % nyqlog_fotf(G,w1,w2) A frequency interval can be specified with
    % w1 (lower) and w2 (upper) if only a certain frequency range is needed.
    %
    % AUTHOR: Andreas Hanssen Moltumyr

    %% Input guard
    if isa(G, 'tf')
        G = fotf(G);
    elseif ~isa(G, 'fotf')
        error(['Argument G has unsupported type ', class(G), '.']);
    end

    switch nargin
        case 1
            [re, im, mag_dB, phase_deg, w] = astep_fotf_freqresp(G);
            if abs(mag_dB(1)) == inf
                re        = re(2:end);
                im        = im(2:end);
                mag_dB    = mag_dB(2:end);
                phase_deg = phase_deg(2:end);
                w         = w(2:end);
            end
            if abs(mag_dB(end)) == inf
                re        = re(1:end-1);
                im        = im(1:end-1);
                mag_dB    = mag_dB(1:end-1);
                phase_deg = phase_deg(1:end-1);
                w         = w(1:end-1);
            end
            phase_rad = (pi/180)*phase_deg;
        case 2
            H = bode(G,w1);
            re        = real(squeeze(H.ResponseData))';
            im        = imag(squeeze(H.ResponseData))';
            mag_dB    = 10.*log10(re.^2 + im.^2);
            phase_rad = atan2(im,re);
            w         = w1;
        case 3
            [re, im, mag_dB, phase_deg, w] = astep_fotf_freqresp(G,w1,w2);
            if abs(mag_dB(1)) == inf
                re        = re(2:end);
                im        = im(2:end);
                mag_dB    = mag_dB(2:end);
                phase_deg = phase_deg(2:end);
                w         = w(2:end);
            end
            if abs(mag_dB(end)) == inf
                re        = re(1:end-1);
                im        = im(1:end-1);
                mag_dB    = mag_dB(1:end-1);
                phase_deg = phase_deg(1:end-1);
                w         = w(1:end-1);
            end
            phase_rad = (pi/180)*phase_deg;
        otherwise
            error('Number of arguments not supported');
    end
    
    plot_logarithmic_nyquist(G, re, im, mag_dB, phase_rad, w);
end


%% Helper functions
function plot_logarithmic_nyquist(G,re,im,mag_dB,phase,w)    
    for j = 1:length(mag_dB)
        if mag_dB(j) == inf
            if j == 1 && length(mag_dB) > 1
                mag_dB(j) = mag_dB(j+1);
            elseif j == length(mag_dB) && length(mag_dB) > 1
                mag_dB(j) = mag_dB(j-1);
            else
                mag_dB(j) = (mag_dB(j+1) + mag_dB(j-1))/2;
            end
        end
    end
    
    %% Detect/find resonace peaks in data
    w_resonance_intervals_start = [];
    w_resonance_intervals_end   = [];
    
    numeric_precission_threshold_gain = 10^3;
    angle_threshold_deg = 10;
      
    on_resonance_peak = false;
    for j = 1:length(w)-1
        if on_resonance_peak == false && ...
        w(j+1) - w(j) <= numeric_precission_threshold_gain*eps(w(j))
            w_resonance_intervals_start = [w_resonance_intervals_start j];
            on_resonance_peak = true;
        elseif on_resonance_peak == true && ...
        w(j+1) - w(j) > numeric_precission_threshold_gain*eps(w(j))
            on_resonance_peak = false;
            w_resonance_intervals_end = [w_resonance_intervals_end j];
        end
    end
    
    w_dim_diff_res_int = length(w_resonance_intervals_start) - length(w_resonance_intervals_end);
    if w_dim_diff_res_int > 0
        w_resonance_intervals_end = [w_resonance_intervals_end, length(w)];
    end

    w_resonance_intervals_start = w_resonance_intervals_start(...
        find(w(w_resonance_intervals_start) ~= 1)); %#ok<FNDSB>
    w_resonance_intervals_end   = w_resonance_intervals_end(...
        find(w(w_resonance_intervals_end) ~= 1)); %#ok<FNDSB>
    
    w_res_int_start = [];
    w_res_int_end   = [];
    for j = 1:length(w_resonance_intervals_start)
        phs_start = phase(w_resonance_intervals_start(j));
        phs_end   = phase(w_resonance_intervals_end(j));
        phs_diff = distance_between_angles(phs_start, phs_end);
        if phs_diff > pi - (pi/180)*angle_threshold_deg
            w_res_int_start = [w_res_int_start, w_resonance_intervals_start(j)];
            w_res_int_end = [w_res_int_end, w_resonance_intervals_end(j)];
        end
    end
    
    %% Remove areas around resonance peaks with low numerical precission
    for j = 1:length(w_res_int_start)
        idx_start = w_res_int_start(j);
        idx_end   = w_res_int_end(j);
        
        re     = [    re(1:idx_start),     re(idx_end:end)];
        im     = [    im(1:idx_start),     im(idx_end:end)];
        mag_dB = [mag_dB(1:idx_start), mag_dB(idx_end:end)];
        phase  = [ phase(1:idx_start),  phase(idx_end:end)];
        w      = [     w(1:idx_start),      w(idx_end:end)];
        
        % Change indicies for remaining intervals in w_res_int_(start/end)
        num_removed_elements = idx_end - idx_start - 1;
        if j ~= length(w_res_int_start)
            w_res_int_start(j+1:end) = w_res_int_start(j+1:end) - num_removed_elements;
        end
        w_res_int_end(j:end)   = w_res_int_end(j:end) - num_removed_elements;
    end
    
    %% Find size of background grid
    mag_db_min = min_ignInf(mag_dB);
    mag_db_max = max_ignInf(mag_dB);
    
    mag_lower_grid_line = 10*floor(mag_db_min/10);
    mag_upper_grid_line = 10*ceil(mag_db_max/10);
    if mag_upper_grid_line < 0
        mag_upper_grid_line = 0;
    end
    if mag_lower_grid_line > 0
        mag_lower_grid_line = 0;
    end
    
    % Calculate value at w = 0
    [a,an,b,bn,~] = fotfdata(G);
    [~, an_li] = min(an);
    [~, bn_li] = min(bn);
    [~, an_hi] = max(an);
    [~, bn_hi] = max(bn);
    
    w_zero = 0;
    mag_at_w_zero   = abs(b(bn_li)/a(an_li))*w_zero^(bn(bn_li) - an(an_li));
    phase_at_w_zero = wrapToPi(pi*(sign(b(bn_li)) ~= sign(a(an_li)))...
                                   + pi/2*(bn(bn_li) - an(an_li)));
    % Calculate value at w = inf
    w_inf  = inf;
    mag_at_w_inf   = abs(b(bn_hi)/a(an_hi))*w_inf^(bn(bn_hi) - an(an_hi));
    phase_at_w_inf = wrapToPi(pi*(sign(b(bn_hi)) ~= sign(a(an_hi)))...
                                  + pi/2*(bn(bn_hi) - an(an_hi)));

    approx_number_of_circles_in_grid = 5;
    approx_distance_between_circles = (mag_upper_grid_line - mag_lower_grid_line)/approx_number_of_circles_in_grid;
    
    
    mag_length_between_lines = 10*floor(approx_distance_between_circles/10);
    %mag_length_between_lines = 20;

    grid_lines_at_dB = mag_lower_grid_line:mag_length_between_lines:mag_upper_grid_line;
    if grid_lines_at_dB(1:end) ~= 0
        grid_lines_at_dB = [grid_lines_at_dB, 0];
        grid_lines_at_dB = sort(grid_lines_at_dB);
    end
    
    center_space_dB = 30;
    
    min_circ_mag = min(grid_lines_at_dB);
    nyquist_curve_offset_dB = center_space_dB - min_circ_mag;
    
    plot_mag_dB = mag_dB + nyquist_curve_offset_dB;
    
    re_dB = plot_mag_dB.*cos(phase);
    im_dB = plot_mag_dB.*sin(phase);
    
    %% semi-circles at inf should all be placed at different hights for better visual impression
    base_inf_offset_dB = 20;
    step_inf_offset_dB = 10;
    inf_offset_dB = base_inf_offset_dB;
    
    %% Add wrap around at inf magnitude for resonance peaks
    for j = 1:length(w_res_int_start)
        idx_start = w_res_int_start(j);
        idx_end   = w_res_int_end(j);
        
        %Decide on rotation direction with symbolic computation
        syms w_start w_end;
        N_w = 100;
        w_sym_start = [w_start, w_start + (w_end - w_start)/N_w];
        w_sym_end   = [w_end + (w_start - w_end)/N_w, w_end];
        
        w_sym_start = subs(w_sym_start, [w_start,w_end], [w(idx_start),w(idx_end)]);
        w_sym_end   = subs(w_sym_end,   [w_start,w_end], [w(idx_start),w(idx_end)]);
        [~, ~, ~, G_phase] = calc_freq_resp_symbolic([w_sym_start, w_sym_end], a,an,b,bn,0);
        
        phase_start_delta = G_phase(2) - G_phase(1);
        phase_end_delta   = G_phase(4) - G_phase(3);
        if phase_start_delta >= 0 && phase_end_delta >= 0
            dir = 1; % CCW
            draw_half_circle = true;
        elseif phase_start_delta <= 0 && phase_end_delta <= 0
            dir  = -1; % CW
            draw_half_circle = true;
        else
            warning(['Direction of rotation for resonance peak ', ...
                     'at w = ', num2str(w(idx_start)), ...
                     ' rad/s can not be decided!', ' Drawing of half '...
                     'circle at inf have been suppressed for this peak.']);
            draw_half_circle = false;
        end

        if draw_half_circle == true
            from_angle = phase(idx_start);
            to_angle   = phase(idx_end);
            angle_step_deg = 1;
            
            if dir == 1 && to_angle > from_angle
                angle_step = (pi/180)*angle_step_deg;
                angles_semi_circle = [from_angle:angle_step:to_angle, to_angle];
            elseif dir == 1 && to_angle < from_angle
                angle_step =  -(pi/180)*angle_step_deg;
                angles_semi_circle = [from_angle:angle_step:to_angle, to_angle];
            elseif dir == -1 && to_angle > from_angle
                angle_step = -(pi/180)*angle_step_deg;
                angles_semi_circle = [from_angle:angle_step:(to_angle-2*pi), (to_angle-2*pi)];
            elseif dir == -1 && to_angle < from_angle
                angle_step = -(pi/180)*angle_step_deg;
                angles_semi_circle = [from_angle:angle_step:to_angle, to_angle];
            else
                error('Can not decide on rotation');
            end
            
            mag_semi_circle = max_ignInf(plot_mag_dB) + inf_offset_dB;
            inf_offset_dB = inf_offset_dB + step_inf_offset_dB;
            re_half_circle = mag_semi_circle.*cos(angles_semi_circle);
            im_half_circle = mag_semi_circle.*sin(angles_semi_circle);
            
            re_dB  = [ re_dB(1:idx_start), re_half_circle,  re_dB(idx_end:end)];
            im_dB  = [ im_dB(1:idx_start), im_half_circle,  im_dB(idx_end:end)];
        
            w_half_circle = ((w(idx_end)+w(idx_start))/2)*ones(size(angles_semi_circle));
            uncalculateable_values = NaN(size(angles_semi_circle));
        
            re     = [    re(1:idx_start), uncalculateable_values,     re(idx_end:end)];
            im     = [    im(1:idx_start), uncalculateable_values,     im(idx_end:end)];
            mag_dB = [mag_dB(1:idx_start), uncalculateable_values, mag_dB(idx_end:end)];
            phase  = [ phase(1:idx_start), uncalculateable_values,  phase(idx_end:end)];
            w      = [     w(1:idx_start),          w_half_circle,      w(idx_end:end)];
            
            inserted_vector_length = length(angles_semi_circle);
        else
            mag_peak_val = max_ignInf(plot_mag_dB) + inf_offset_dB;
            inf_offset_dB = inf_offset_dB + step_inf_offset_dB;
            re_peak_fix_val_start = mag_peak_val.*cos(phase(idx_start));
            im_peak_fix_val_start = mag_peak_val.*sin(phase(idx_start));
            re_peak_fix_val_end   = mag_peak_val.*cos(phase(idx_end));
            im_peak_fix_val_end   = mag_peak_val.*sin(phase(idx_end));
        
            re_peak_subst = [re_peak_fix_val_start, inf, re_peak_fix_val_end];
            im_peak_subst = [im_peak_fix_val_start, inf, im_peak_fix_val_end];
            
            re_dB  = [ re_dB(1:idx_start-1), re_peak_subst,  re_dB(idx_end+1:end)];
            im_dB  = [ im_dB(1:idx_start-1), im_peak_subst,  im_dB(idx_end+1:end)];
        
            uncalculateable_values = NaN(size(re_peak_subst));
            w_peak_subst = w(idx_start).*ones(size(re_peak_subst));
        
            re     = [    re(1:idx_start), uncalculateable_values,     re(idx_end:end)];
            im     = [    im(1:idx_start), uncalculateable_values,     im(idx_end:end)];
            mag_dB = [mag_dB(1:idx_start), uncalculateable_values, mag_dB(idx_end:end)];
            phase  = [ phase(1:idx_start), uncalculateable_values,  phase(idx_end:end)];
            w      = [     w(1:idx_start),           w_peak_subst,      w(idx_end:end)];
            
            inserted_vector_length = length(re_peak_subst);
        end
        
        % Change indicies for remaining intervals in w_res_int_(start/end)
        if j ~= length(w_res_int_start)
            w_res_int_start(j+1:end) = w_res_int_start(j+1:end) + inserted_vector_length;
        end
        w_res_int_end(j:end) = w_res_int_end(j:end) + inserted_vector_length;
    end
    
    %% Add point in zero or quarter circle at inf
    if mag_at_w_zero == 0
        re_dB  = [0,     re_dB];
        im_dB  = [0,     im_dB];
        
        w      = [0,         w];
        re     = [0,        re];
        im     = [0,        im];
        phase  = [NaN,   phase];
        mag_dB = [-inf, mag_dB];
        
    elseif mag_at_w_zero == inf
        if sign(b(bn_li)/a(an_li)) == 1
            from_angle = 0;
        else
            from_angle = pi;
        end
        to_angle = phase_at_w_zero;
        angle_step_deg = 1;
        if from_angle <= to_angle
            angle_step = (pi/180)*angle_step_deg;
        else
            angle_step = -(pi/180)*angle_step_deg;
        end
        angles_quarter_circle = from_angle:angle_step:to_angle;
        mag_quarter_circle = max_ignInf(plot_mag_dB) + inf_offset_dB;
        inf_offset_dB = inf_offset_dB + step_inf_offset_dB;
        re_quarter_circle = mag_quarter_circle.*cos(angles_quarter_circle);
        im_quarter_circle = mag_quarter_circle.*sin(angles_quarter_circle);
        re_dB = [re_quarter_circle, re_dB];
        im_dB = [im_quarter_circle, im_dB];
        
        re     = [inf.*cos(angles_quarter_circle),       re];
        im     = [inf.*sin(angles_quarter_circle),       im];
        mag_dB = [inf(size(angles_quarter_circle)),  mag_dB];
        phase  = [NaN(size(angles_quarter_circle)),   phase];
        w      = [zeros(size(angles_quarter_circle)),      w];
    end
    
    if mag_at_w_inf == 0
        re_dB  = [re_dB,     0];
        im_dB  = [im_dB,     0];
        
        w      = [w,       inf];
        re     = [re,        0];
        im     = [im,        0];
        phase  = [phase,   NaN];
        mag_dB = [mag_dB, -inf];
        
    elseif mag_at_w_inf == inf
        if sign(b(bn_hi)/a(an_hi)) == 1
            to_angle = 0;
        else
            to_angle = pi;
        end
        from_angle = phase_at_w_inf;
        angle_step_deg = 1;
        if from_angle <= to_angle
            angle_step = (pi/180)*angle_step_deg;
        else
            angle_step = -(pi/180)*angle_step_deg;
        end
        angles_quarter_circle = from_angle:angle_step:to_angle;
        mag_quarter_circle = max_ignInf(plot_mag_dB) + inf_offset_dB;
        inf_offset_dB = inf_offset_dB + step_inf_offset_dB;
        re_quarter_circle = mag_quarter_circle.*cos(angles_quarter_circle);
        im_quarter_circle = mag_quarter_circle.*sin(angles_quarter_circle);
        re_dB  = [re_dB, re_quarter_circle];
        im_dB  = [im_dB, im_quarter_circle];
        
        re     = [re,      inf.*cos(angles_quarter_circle)];
        im     = [im,      inf.*sin(angles_quarter_circle)];
        mag_dB = [mag_dB, inf(size(angles_quarter_circle))];
        phase  = [phase,  NaN(size(angles_quarter_circle))];
        w      = [w,      inf(size(angles_quarter_circle))];
    end
    
    
    %% Plot Nyquist curve
    nyquist_curve_pw_handle = plot(re_dB,  im_dB,...
                        'Color',     [0, 0.447, 0.741],...
                        'LineStyle', '-',...
                        'LineWidth', 1.5);
    hold on;

    nyquist_curve_nw_handle = plot(re_dB, -im_dB,...
                                   'Color',     [0, 0, 0],...
                                   'LineStyle', ':',...
                                   'LineWidth', 1.5);
    hold off; axis off; axis equal;
    
    set(nyquist_curve_pw_handle, 'UserData', 'positive_nyquist_curve');
    set(nyquist_curve_nw_handle, 'UserData', 'negative_nyquist_curve');
    
    plot_grid(grid_lines_at_dB, center_space_dB);
    
    fig_handle = gcf;
    set(fig_handle, 'Color', [0.98, 0.98, 0.98]);
    
    %% Add Data Cursor calculation function (Enables readout of values on the curve)
    dcm_obj = datacursormode(fig_handle);
    set(dcm_obj, 'UpdateFcn', {@fotf_nyqlog_dcm_update_fcn,...
                               w, mag_dB, phase, re, im});
    set(dcm_obj, 'DisplayStyle', 'datatip');
    set(dcm_obj, 'Enable', 'on');
    
    title({'Logarithmic Nyquist Diagram'; ' '});
end

function plot_grid(circ_num_and_mag, center_space)
    min_circ_mag = min(circ_num_and_mag);
    print_circ_num_and_mag = circ_num_and_mag - min_circ_mag + center_space;
    
    % Options
    circle_resolution                  = 1024;
    num_straight_lines                 = 8;
    start_offset_degree_straight_lines = 0;
    magnitude_text_angle               = -22.5;
    magnitude_text_distance_from_line  = 1;
    phase_text_distance_from_circle    = 4;
    text_FontSize                      = 8;
    
    hold on;
    % Draw circles at given magnitudes
    omega = 0:2*pi/circle_resolution:2*pi;
    for j = 1:length(print_circ_num_and_mag)
        x = print_circ_num_and_mag(j).*cos(omega);
        y = print_circ_num_and_mag(j).*sin(omega);
        circle_handle = plot(x,y, 'r-.');
        set(circle_handle, 'PickableParts', 'none');
        if (circ_num_and_mag(j) == 0)
            set(circle_handle, 'LineStyle', '-');
            set(circle_handle, 'LineWidth', 1.0);
        end
    end
    
    % Draw straight lines streaching outwards from the origin
    line_angles = start_offset_degree_straight_lines ...
                + ((0:num_straight_lines-1)./num_straight_lines).*360;
    mag_lower  = min(print_circ_num_and_mag);
    mag_higher = max(print_circ_num_and_mag);
    for j = 1:num_straight_lines
        x = [mag_lower, mag_higher].*cos((pi/180)*line_angles(j));
        y = [mag_lower, mag_higher].*sin((pi/180)*line_angles(j));
        line_handle = plot(x,y, 'r-.');
        set(line_handle, 'PickableParts', 'none');
    end
    
    % Draw Magnitude labels for grid
    for j = 1:length(circ_num_and_mag)
        x = (print_circ_num_and_mag(j)+magnitude_text_distance_from_line)...
            .*cos((pi/180)*magnitude_text_angle);
        y = (print_circ_num_and_mag(j)+magnitude_text_distance_from_line)...
            .*sin((pi/180)*magnitude_text_angle);
        txt = [num2str(circ_num_and_mag(j)), ' dB'];
        text(x, y, txt, 'FontSize', text_FontSize);
    end
    
    % Draw Phase labels for grid
    degree_symbol = char(176);
    for j = 1:num_straight_lines
        angle = (pi/180)*line_angles(j);
        x = (mag_higher+phase_text_distance_from_circle).*cos(angle);
        y = (mag_higher+phase_text_distance_from_circle).*sin(angle);
        txt = [num2str(line_angles(j)), degree_symbol];
        text_handle = text(x, y, txt, 'FontSize', text_FontSize);
        if 0.05 < cos(angle)
            set(text_handle, 'HorizontalAlignment', 'left');
        elseif cos(angle) < -0.05
            set(text_handle, 'HorizontalAlignment', 'right');
        else
            set(text_handle, 'HorizontalAlignment', 'center');
        end
    end
    
    % Draw -1 point (mag=0dB and phase = -180)
    x = (0+ center_space - min_circ_mag).*cos(pi);
    y = (0+ center_space - min_circ_mag).*sin(pi);
    minus_one_point_handle = plot(x,y, 'Color',     [0, 0, 0], ...
              'LineStyle', 'none', ...
              'Marker',    'o', ...
              'LineWidth', 1.5);
    set(minus_one_point_handle, 'PickableParts', 'none');
    hold off;
end

function max_val = max_ignInf(val)
    max_val = -inf;
    for j = 1:length(val)
        if val(j) > max_val && val(j)~= inf
            max_val = val(j);
        end
    end
end

function min_val = min_ignInf(val)
    min_val = inf;
    for j = 1:length(val)
        if val(j) < min_val && val(j)~= -inf
            min_val = val(j);
        end
    end
end

function diff = distance_between_angles(angle1, angle2)
    psi = rem(abs(angle1 - angle2),2*pi);
    if psi > pi
        diff = 2*pi - psi;
    else
        diff = psi;
    end
end

function [re, im, mag, phase] = calc_freq_resp_symbolic(w, a, na, b, nb, L)
    syms s;
    freq_resp_num = b*(s.^nb.');
    freq_resp_den = a*(s.^na.');
    freq_resp = freq_resp_num/freq_resp_den;
    if L > 0
        freq_resp = freq_resp.*exp(-L*s);
    end
    
    freq_resp = subs(freq_resp, s, 1i.*w);
    
    phase = double(angle(freq_resp));
    mag   = double(abs(freq_resp));
    re    = double(real(freq_resp));
    im    = double(imag(freq_resp));
end
