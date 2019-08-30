%% FILE: assess_stability.m
%% AUTHOR: Andreas Hanssen Moltumyr

function [isStable, GMs, GMs_dB, neg_sgn_GMs_dB, PMs, data] = assess_stability(G,Np,w,re,im)
% assess_stability(G,Np,w) - Gives info about the stability of the
% transfer function G. Np is the number of poles in right half plane (rhp)
% of G. w is an array of frequencies used when calculating the frequency
% response of G during the stability assessment.
% assess_stability(G,Np,{w1,w2}) - Instead of giving an array of
% frequencies a lower (w1) and upper (w2) value of interesting frequencies
% can be given and the assess_stability will adaptively choose a set of
% frequencies to calculate the frequency respose.
% assess_stability(G,Np,w,re,im) - Another option is to give a set of
% frequencies and the associated frequency response in re and im.

switch nargin
    case 2
        [REAL,IMAG,~,~,w] = astep_fotf_freqresp(G);
    case 3
        %% Calculate frequency response when arguments re and im are not given
        if isa(w,'cell') && length(w) == 2
            w1 = w{1};
            w2 = w{2};
            if w1 > w2
                error('When passing w = {w1,w2}, w1 must be greater than w2');
            end
            [REAL,IMAG,~,~,w] = astep_fotf_freqresp(G, w1, w2);
        else
            H = bode(G,w);
            [REAL,IMAG,w] = nyquist(H);
            REAL = squeeze(REAL);
            IMAG = squeeze(IMAG);
        end
    case 5
        REAL = re;
        IMAG = im;
    otherwise
        error('Number of arguments unsupported')
end

%% 1.) Find real axis crossings
[re_x_val, w_x_val, x_dir] = find_zero_crossings(REAL,IMAG,w);

%% 2.) Add data for w = 0 and w = inf
[a,an,b,bn,~] = fotfdata(G);

[an_lowest,  an_lowest_index]  = min(an);
[bn_lowest,  bn_lowest_index]  = min(bn);
[an_highest, an_highest_index] = max(an);
[bn_highest, bn_highest_index] = max(bn);

% Check if system is improper, biproper or strictly proper
% Calculate value at w = 0
w_x_val_w_zero  = 0;
if an_lowest > bn_lowest
    re_x_val_w_zero = sign(b(bn_lowest_index))*sign(a(an_lowest_index))*inf; % Can be positive or negative, has to check this
elseif an_lowest < bn_lowest
    re_x_val_w_zero = 0;
else
    re_x_val_w_zero = b(bn_lowest_index)/a(an_lowest_index);
end

% Calculate value at w = inf
w_x_val_w_inf  = inf;
if an_highest > bn_highest
    % Strictly proper
    re_x_val_w_inf = 0;
elseif an_highest < bn_highest
    % Improper
    re_x_val_w_inf = sign(b(bn_highest_index))*sign(a(an_highest_index))*inf;
else
    % Biproper
    re_x_val_w_inf = b(bn_highest_index)/a(an_highest_index);
end

if isempty(x_dir)
    dPsi_w_zero = calc_phase_change(G,0);
    dPsi_w_inf  = calc_phase_change(G,inf);
    dPsi_w_ninf  = calc_phase_change(G,-inf);
else
    x_dir_w_zero = -x_dir(1);
    x_dir_w_inf  = -x_dir(end);
end

re_x_val = [re_x_val_w_zero, re_x_val, re_x_val_w_inf];
w_x_val  = [ w_x_val_w_zero,  w_x_val,  w_x_val_w_inf];
x_dir    = [   x_dir_w_zero,    x_dir,    x_dir_w_inf];

%% 3.) Calculate encirclements
encircs = zeros(1,length(w_x_val)+1);

% Create lookup table to map frequency order with real axis crossing order
lutbl = zeros(4,length(re_x_val));
lutbl(1,:) = 1:length(re_x_val);
lutbl(2,:) = re_x_val;
lutbl(3,:) = w_x_val;
lutbl(:,:) = sortrows(lutbl(:,:)',2)';
lutbl(4,:) = 1:length(re_x_val);
lutbl_2(:,:) = sortrows(lutbl(:,:)',1)';

for j = 1:length(w_x_val)-1
    % Find which values in encircs that should add rotations
    k1 = lutbl_2(4,j);
    k2 = lutbl_2(4,j+1);
    if k1 > k2
        k_values = k2:k1;
    else
        k_values = k1:k2;
    end
    k_values = k_values(1:end-1);
    
    % Check real axis crossing direction and if the value from one crossing
    % is lesser than or greater than the next. Add encirclements according
    % to these checks.
    if (x_dir(j) == -1 && re_x_val(j) > re_x_val(j+1)) || ...
       (x_dir(j) ==  1 && re_x_val(j) < re_x_val(j+1))
        encircs(k_values+1) = encircs(k_values+1) - 1;
    elseif (x_dir(j) ==  1 && re_x_val(j) > re_x_val(j+1)) || ...
           (x_dir(j) == -1 && re_x_val(j) < re_x_val(j+1))
        encircs(k_values+1) = encircs(k_values+1) + 1;
    end
end

%% 4.) Find stable regions
region_data_encircs           = encircs;
region_data_lower_real_value  = [-inf, lutbl(2,1:end)];
region_data_higher_real_value = [lutbl(2,1:end),  inf];

if region_data_lower_real_value(1)  == -inf && ...
   region_data_higher_real_value(1) == -inf
    region_data_encircs           = region_data_encircs(2:end);
    region_data_lower_real_value  = region_data_lower_real_value(2:end);
    region_data_higher_real_value = region_data_higher_real_value(2:end);
    
end
if region_data_lower_real_value(end)  == inf && ...
   region_data_higher_real_value(end) == inf
    region_data_encircs           = region_data_encircs(1:end-1);
    region_data_lower_real_value  = region_data_lower_real_value(1:end-1);
    region_data_higher_real_value = region_data_higher_real_value(1:end-1);
end

stable_regions = [];
for j = 1:length(region_data_encircs)
    if region_data_encircs(j) == Np
        stable_regions = [stable_regions,...
                          [region_data_lower_real_value(j);...
                           region_data_higher_real_value(j)]];
    end
end

% Return information of the regions and their encirclements as value data
data = [region_data_encircs; region_data_lower_real_value; region_data_higher_real_value];

%% 5.) Check stability based on given Np value
make_stable_gain_pairs = [];
isStable = false;

for j = 1:size(stable_regions,2)
    k1 = stable_regions(1,j);
    k2 = stable_regions(2,j);
    e1 = 1/abs(k1);
    e2 = 1/abs(k2);
    if sign(k1) == -1 && sign(k2) == -1
        g_lower = e1;
        g_upper = e2;
        make_stable_gain_pairs = [make_stable_gain_pairs, [g_lower; g_upper]];
    elseif sign(k1) == -1 && sign(k2) == 1
        g_lower = e1;
        g_upper = inf;
        make_stable_gain_pairs = [make_stable_gain_pairs, [g_lower; g_upper]];
        
        g_lower = -e2;
        g_upper = 0;
        make_stable_gain_pairs = [make_stable_gain_pairs, [g_lower; g_upper]];
    elseif sign(k1) == 1 && sign(k2) == 1
        g_lower = -e2;
        g_upper = -e1;
        make_stable_gain_pairs = [make_stable_gain_pairs, [g_lower; g_upper]];
    else
        if k1 < inf && k2 < inf
            error(['k1 > k2, is not possible, ', 'k1 = ', num2str(k1), ', k2 = ', num2str(k2)]);
        else
            
        end
    end
    
    if k1 < -1 && -1 < k2
        % Sets the isStable output value to true if some of the stable
        % regions contain the -1 point.
        isStable = true;
    end
end

%% 6.) Calculate gain margins
GMs_dB         = [];
neg_sgn_GMs_dB = [];
for j = 1:size(make_stable_gain_pairs,2)
    if sign(make_stable_gain_pairs(1,j)) == -1 || sign(make_stable_gain_pairs(2,j)) == -1
        neg_sgn_GMs_dB = [neg_sgn_GMs_dB, 20*log10(abs(make_stable_gain_pairs(:,j)))];
    else
        GMs_dB         = [GMs_dB, 20*log10(abs(make_stable_gain_pairs(:,j)))];
    end
end

GMs = make_stable_gain_pairs;

%% 7.) Calculate phase margins
[phase0, w0] = find_unity_gain_crossings(REAL,IMAG,w);

phase_margins = sort(180-abs(phase0),'ascend');

if isStable == true
    PMs = phase_margins;
else
    PMs = [];
end

end % assess_stability()

%% Helper functions
function [phase0, w0] = find_unity_gain_crossings(re, im, w)
    unity_gain_crossing_pairs = [];
    for j = 1:length(re)-1
        z1 = abs(re(j)  +1i*im(j));
        z2 = abs(re(j+1)+1i*im(j+1));
        if (z1 < 1 && z2 > 1) || (z1 > 1 && z2 < 1)
            new_unity_gain_crossing_pair = [j; j+1];
            unity_gain_crossing_pairs = [unity_gain_crossing_pairs, ...
                                         new_unity_gain_crossing_pair];
        end
    end

    unity_gain_crossing_data = zeros(2,size(unity_gain_crossing_pairs,2));
    for j = 1:size(unity_gain_crossing_data,2)
        % Find approximate unity gain crossing through linear interpolation
        x1 = re(unity_gain_crossing_pairs(1,j));
        y1 = im(unity_gain_crossing_pairs(1,j));
        x2 = re(unity_gain_crossing_pairs(2,j));
        y2 = im(unity_gain_crossing_pairs(2,j));
        
        x  = (x1+x2)/2;
        y  = (y1+y2)/2;
        
        % Find frequency value at approximate unity crossing point
        w1 = w(unity_gain_crossing_pairs(1,j));
        w2 = w(unity_gain_crossing_pairs(2,j));
        
        w_x = (w1+w2)/2;
        unity_gain_crossing_data(2,j) = w_x;
        
        phase = angle(x+1i*y);
        unity_gain_crossing_data(1,j) = phase;
    end
    
    phase0 = (180/pi)*unity_gain_crossing_data(1,:);
    w0     = unity_gain_crossing_data(2,:);
end


function [re0, w0, crossing_direction] = find_zero_crossings(re,im,w)
    % Find approximate crossings of real axis
    sign_shift_pairs = [];
    for i = 1:length(im)-1
        if (im(i) > 0 && im(i+1) < 0)
            % Detected falling sign shift (from + to -)
            new_sign_shift_pair = [i; i+1; -1];
            sign_shift_pairs = [sign_shift_pairs, new_sign_shift_pair];
        elseif (im(i) < 0 && im(i+1) > 0)
            % Detected rising sign shift (from - to +)
            new_sign_shift_pair = [i; i+1; 1];
            sign_shift_pairs = [sign_shift_pairs, new_sign_shift_pair];
        end
    end

    % Linear interpolation of real axis crossings
    zero_crossing_data = zeros(3,size(sign_shift_pairs,2));
    for i = 1:size(zero_crossing_data,2)
        % Find approximate real axis crossing
        x1 = re(sign_shift_pairs(1,i));
        y1 = im(sign_shift_pairs(1,i));
        x2 = re(sign_shift_pairs(2,i));
        y2 = im(sign_shift_pairs(2,i));
        
        x  = x1 - y1*((x2-x1)/(y2-y1));
        zero_crossing_data(1,i) = x;
        
        % Find frequency value at crossing point
        w1 = w(sign_shift_pairs(1,i));
        w2 = w(sign_shift_pairs(2,i));
        
        w_x = w1 - y1*((w2-w1)/(y2-y1));
        zero_crossing_data(2,i) = w_x;
        
        % Add info about crossing direction
        zero_crossing_data(3,i) = sign_shift_pairs(3,i);
    end
    
    re0 = zero_crossing_data(1,:);
    w0  = zero_crossing_data(2,:);
    crossing_direction = zero_crossing_data(3,:);
end


function dPsi = calc_phase_change(G,w_val)
    [a,an,b,bn,~] = fotfdata(G);
    w = w_val;
    
    % Calculate real and imag part of nominator and denominator, of fotf G
    B_mag = b.*(w.^bn);
    B_theta = pi/2*bn;
    Br = B_mag*cos(B_theta)';
    Bi = B_mag*sin(B_theta)';
    
    A_mag = a.*(w.^an);
    A_theta = pi/2*an;
    Ar = A_mag*cos(A_theta)';
    Ai = A_mag*sin(A_theta)';
    
    % Calculate derivative with respect of w for Br, Bi, Ar and Ai
    dB_mag = b.*bn.*(w.^(bn-1));
    dBr = dB_mag*cos(B_theta)';
    dBi = dB_mag*sin(B_theta)';
    
    dA_mag = a.*an.*(w.^(an-1));   
    dAr = dA_mag*cos(A_theta)';
    dAi = dA_mag*sin(A_theta)';
    
    % Find derivative of phase with respect to w
    F1  = Ar*Br + Ai*Bi;
    dF1 = dAr*Br + Ar*dBr + dAi*Bi + Ai*dBi;
    F2  = Ar*Bi - Ai*Br;
    dF2 = dAr*Bi + Ar*dBi - dAi*Br - Ai*dBr;
    
    dPsi = (dF2*F1 - F2*dF1)/(F1^2 + F2^2);
end
