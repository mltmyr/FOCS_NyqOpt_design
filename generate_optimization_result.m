%% FILE: generate_optimization_result.m
%% AUTHOR: Andreas Hanssen Moltumyr

%% Configuration/Options
show_plots = true;
save_results_and_plots = true;

%% Process configs/opts
if show_plots
    plot_visibility = 'on';
else
    plot_visibility = 'off';
end

%% Get open-loop (L), closed-loop/complementary sensitivity fnc (T) and output sensitivity fnc (So)
[L, T, So, misc] = construct_system_fotfs(G, regulator_type, x_optimal);

%% Check regulator open loop stability so to adjust needed encirclement number Np
    % N_p: Needed encirclements in nyquist diagram for closed-loop stability
    switch regulator_type
        case 'FO-PPF_2'
            zeta_d = x_optimal(3);
            beta_d = x_optimal(5);
            if zeta_d > -cos(pi/2*beta_d)
                N_p = num_unstable_poles_open_loop;
            else
                N_p = num_unstable_poles_open_loop + 2;
            end
        case 'FO-PPF_3'
            zeta_d = x_optimal(4);
            beta_d = x_optimal(6);
            if zeta_d > -cos(pi/2*beta_d)
                N_p = num_unstable_poles_open_loop;
            else
                N_p = num_unstable_poles_open_loop + 2;
            end
        otherwise
            N_p = num_unstable_poles_open_loop;
    end

%% Assess stability (Nyquist Criterion Approach)
[L_isStable, L_GMs, L_GMs_dB, L_neg_sgn_GMs_dB, L_PMs, L_data] = assess_stability(L,N_p);
disp(['Stable: ', num2str(L_isStable)]);
disp('Gain Margins:'); disp(L_GMs);
disp('Gain Margins(0deg) [dB]:'); disp(L_GMs_dB);
disp('Gain Margins(180deg) [dB]:'); disp(L_neg_sgn_GMs_dB);
if length(L_PMs) >= 1
disp(['Phase Margin: ', num2str(L_PMs(1))]);
end

%% Plot Logarithmic Nyquist Diagram of open-loop transfer function (L)
nyq_fig_handle  = figure('visible', plot_visibility);
nyqlog_fotf(L);

%% Plot Bode Diagram of plant, open- and closed-loop transfer functions
bode_fig_handle = figure('visible', plot_visibility);
w = logspace(-1, 5, 10000);
bode(G, w); hold on;
bode(L, w);
bode(T, w); hold off;
legend('G(s)', 'L(s)', 'T(s)'); grid on;

%% Simulink Simulation of System with Regulator given in regulator_type
simStepAmplitude         = 1;
simStepSimulationTime    = 1;

simFractionalFilterOrder      = 21;
%simFractionalOrderFrequencies = [1e-10 1e6];
simFractionalOrderFrequencies = [1e-5 1e5];

switch regulator_type
    case {'PPF', 'FO-PPF_1', 'FO-PPF_2', 'FO-PPF_3'}
        C_t = misc(1);
        C_d = misc(2);
        sim('step_response_simulation_foppf.slx');
    case {'PID', 'FO-PID'}
        C = misc(1);
        sim('step_response_simulation_fopid.slx');
    otherwise
            error('Choose a supported regulator or implement support for new regulator');
end

t_sim  = sim_data.Time;
y_sim  = sim_data.Data(:,1);
y0_sim = sim_data.Data(:,2);
u_sim  = sim_data.Data(:,3);

color_orange = [0.8500, 0.3250, 0.0980];
color_blue   = [0, 0.4470, 0.7410];

sim_fig1_handle = figure('visible', plot_visibility);
plot(t_sim, y0_sim, 'color', color_orange); hold on;
plot(t_sim, y_sim, 'color', color_blue); hold off;
grid on; xlabel('Time [sec]'); ylabel('y(t) [V]');

%% Save Figures and Data
if save_results_and_plots == true

% Store results in folder 'Results', create directory if it does not exist
results_path = fullfile(pwd, 'Results');
if ~exist(results_path, 'file')
       mkdir(results_path);
end

% Store results for different regulators in different folders. Create
% folder if it does not exist
reg_path = fullfile(results_path, regulator_type);
if ~exist(reg_path, 'file')
       mkdir(reg_path);
end

% A metadata.mat file is put into all the different regulator folders to
% add a counter. The counter is used to differentiate the different runs so
% results from multiple runs can be saved.
reg_metadata_path = fullfile(reg_path, 'metadata.mat');
if exist(reg_metadata_path, 'file')
    load(reg_metadata_path, 'num_results_saved');
else
    num_results_saved = 0;
end
num_results_saved = num_results_saved + 1;
save(reg_metadata_path, 'num_results_saved');

% Save the different plots in fig, eps and png format
new_results_folder_path = fullfile(reg_path, ['Run', num2str(num_results_saved)]);
mkdir(new_results_folder_path);

saveas(nyq_fig_handle,  fullfile(new_results_folder_path, 'nyqlog_plot'), 'fig');
saveas(bode_fig_handle, fullfile(new_results_folder_path, 'bode_plot'),   'fig');
saveas(sim_fig1_handle, fullfile(new_results_folder_path, 'step_plot'),   'fig');

saveas(nyq_fig_handle,  fullfile(new_results_folder_path, 'nyqlog_plot'), 'epsc');
saveas(bode_fig_handle, fullfile(new_results_folder_path, 'bode_plot'),   'epsc');
saveas(sim_fig1_handle, fullfile(new_results_folder_path, 'step_plot'),   'epsc');

saveas(nyq_fig_handle,  fullfile(new_results_folder_path, 'nyqlog_plot'), 'png');
saveas(bode_fig_handle, fullfile(new_results_folder_path, 'bode_plot'),   'png');
saveas(sim_fig1_handle, fullfile(new_results_folder_path, 'step_plot'),   'png');

% List of variables to try to save to new_results_data_path
vars_to_save = {'G', 'L', 'So', 'misc', ...
    'regulator_type', 'x_optimal', 'lower_bound_x', 'upper_bound_x', ...
    'num_unstable_poles_open_loop', 'L_isStable', 'L_GMs', ...
    'L_GMs_dB', 'L_neg_sgn_GMs_dB', 'L_PMs', 'L_data', ...
    'simStepAmplitude', 'simStepSimulationTime', ...
    'simFractionalFilterOrder', 'simFractionalOrderFrequencies', ...
    'Population', 'optimization_elapsed_time'};

% Guarding against the attempted storage of variables that does not exist
for k = 1:length(vars_to_save)
    if exist(vars_to_save{k}, 'var')
        optimization_and_result_data.(vars_to_save{k}) = eval(vars_to_save{k});
    end
end

% Save data to data.mat file
new_results_data_path = fullfile(new_results_folder_path, 'data');
save(new_results_data_path, '-struct', 'optimization_and_result_data');

end
