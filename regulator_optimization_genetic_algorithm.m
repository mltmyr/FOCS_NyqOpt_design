%% FILE: regulator_optimization_genetic_algorithm.m
%% PURPOSE: Run optimization based tuning of controller parameters.
%% AUTHOR: Andreas Hanssen Moltumyr

%% Load Plant Model
ident_tfs = load('plant_models/afm_identified_tfs.mat');
G = ident_tfs.tf9;
G = fotf(G);

%% Choice of controller/regulator to control plant, and optimize
%regulator_type = 'PPF';
%regulator_type = 'FO-PPF_1';
%regulator_type = 'FO-PPF_2';
%regulator_type = 'FO-PPF_3';
regulator_type = 'FO-PID';
%regulator_type = 'PID';

%% Regulator Options
switch regulator_type
    case 'PPF'
        num_unstable_poles_open_loop = 0;
        %               [   k_t,    k_d,  zeta_d,  omega_d]
        lower_bound_x = [     1,      1,       0,   1*10^3];
        upper_bound_x = [1*10^3, 1*10^3,       5,   3*10^4];
    case 'FO-PPF_1'
        num_unstable_poles_open_loop = 0;
        %               [   k_t, alpha_t,    k_d,  zeta_d,  omega_d]
        lower_bound_x = [     1,    0.01,      1,       0,   5*10^2];
        upper_bound_x = [1*10^4,     2.0, 1*10^3,      10,   2*10^4];
    case 'FO-PPF_2'
        num_unstable_poles_open_loop = 0;
        %               [   k_t,     k_d, zeta_d, omega_d,   beta_d]
        lower_bound_x = [     1, 1*10^-3,    0.1,  5*10^2,     0.01];
        upper_bound_x = [1*10^5,  1*10^3,      6,  1*10^6,     1.99];
    case 'FO-PPF_3'
        num_unstable_poles_open_loop = 0;
        %               [    k_t, alpha_t,     k_d, zeta_d, omega_d, beta_d]
        lower_bound_x = [1*10^-2,    0.9,  1*10^-3,    0.1,  1*10^0,   0.01];
        upper_bound_x = [ 1*10^6,    1.99,  1*10^3,     10,  1*10^6,   1.99];
    case 'FO-PID'
        num_unstable_poles_open_loop = 0;
        %               [    k_p,     k_i,     k_d, mu_i,  mu_d,    tau_f]
        lower_bound_x = [1*10^-9,  1*10^-2, 1*10^-6, 0.01,  0.01,   1*10^2];
        upper_bound_x = [ 1*10^2,  1*10^7,  2*10^6, 1.99,  1.99,   5*10^5];
    case 'PID'
        num_unstable_poles_open_loop = 0;
        %               [    k_p,     k_i,     k_d,    tau_f]
        lower_bound_x = [1*10^-4, 1*10^-4, 1*10^-4,  1*10^-4];
        upper_bound_x = [ 1*10^4,  1*10^5,  1*10^4,   1*10^2];
    otherwise
        error('Choose a supported regulator or implement support for new regulator');
end

%% Set Objectivity and Constraint functions
objective_function  = @(x) regulator_optimization_opt(x,G,regulator_type);
constraint_function = @(x) regulator_optimization_con(x,G,regulator_type,num_unstable_poles_open_loop);

%% Genetic Algorithm (GA) optimization Options
population_size       = 100;

% Reproduction Options: Controls number of elites, crossover and mutation
%                       individuals in next generation.
num_elite_individuals = ceil(0.05*population_size);
crossover_fraction    = 0.6;

% Stopping Conditions
max_num_generations   = 10;

max_stall_generations = 10;
function_tolerance    = 1e-20;

fitness_limit         = -inf;

time_limit            = inf;

stall_time_limit      = inf;

% Use previous end population as start population
use_prev_pop_as_start_pop = true;

% Use a previously saved random seed to enable repeatability
use_previous_rng_state = false;

%%
if use_prev_pop_as_start_pop == true
    initial_population = POPULATION;
else
    initial_population = [];
end
initial_population_range = [lower_bound_x; upper_bound_x];

ga_options = optimoptions(@ga,...
                'Display', 'diagnose',...
                'PlotFcn', {@gaplotbestf},...
                'PlotInterval', 1,...
                'PopulationSize',           population_size,...
                'InitialPopulationMatrix',  initial_population,...
                'InitialPopulationRange',   initial_population_range,...
                ...
                'FitnessScalingFcn',        @fitscalingrank,...
                'SelectionFcn',             @selectionstochunif,...
                'CrossoverFcn',             @crossoverscattered,...
                'MutationFcn',              @mutationadaptfeasible, ...
                'EliteCount',               num_elite_individuals,...
                'CrossoverFraction',        crossover_fraction,...
                'NonlinearConstraintAlgorithm', 'auglag',...
                ...
                'MaxGenerations',           max_num_generations,...
                'MaxStallGenerations',      max_stall_generations,...
                'FunctionTolerance',        function_tolerance,...
                'FitnessLimit',             fitness_limit,...
                'MaxTime',                  time_limit,...
                'MaxStallTime',             stall_time_limit...
                );

%% Run Optimization with Genetic Algorithm
if use_previous_rng_state == true
    rng(rng_saved_state);
end
rng_saved_state = rng;

length_x = length(lower_bound_x);

start_optimization = tic;
[x_optimal,FVAL,EXITFLAG,OUTPUT,POPULATION,SCORES] = ga(...
                                 objective_function, length_x,...
                                 [],[],[],[],...
                                 lower_bound_x,upper_bound_x,...
                                 constraint_function,...
                                 ga_options);
optimization_elapsed_time = toc(start_optimization);
