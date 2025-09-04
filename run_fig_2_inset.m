%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: run_fig_2_inset.m
% Purpose: Generate Figure 2-inset from the paper 
%          "Leakage at Zero Temperature from Changes in Chemical Potential
%           in Majorana Qubits"
%
% Author: Marcus C. Goffage
% Affiliation: University of New South Wales
%
% Paper: "Decoherence in Majorana Qubits by 1/f Noise"
% Paper Authors: M. C. Goffage^1, A. Alase^2, M. C. Cassidy^1, 
%                S. N. Coppersmith^1
% Affiliations: ^1 University of New South Wales
%               ^2 University of Sydney 
%               
%
% -------------------------------------------------------------------------
% ABOUT THIS SCRIPT
% -------------------------------------------------------------------------
% Description:
%   - Generates Figure 2 Inset. 
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves Figure_2b in the ./results
%   - Saves entire worskpace in the ./results directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running run_fig_2_inset.m');

%% Environment Setup
clear all; 
%close all;
addpath('QPP_Library_P1_submit')


%% Simulation Parameters
% Ramp properties
ramp_height = 0.03;             % ramp_height = mu_fin - mu_in 
N_vec = 2:2:100;               % Kitaev chain lengths 
ramp_rate_constant = 2e-2;      % Constant ramp rate

% Tetron parameters
mu_init   = 0;       
mu_offset = 0.0; % DC offset in chemical potential between the chains    
w         = 0.5;     
delta     = 0.5;     
BC        = "OBC";  

% Simulation time
t_init           = 0;
t_final_no_ramp  = 100; 

% Noise / timestep parameters
num_trials             = 1;       % Keep fixed at 1
delta_mu_time_step_max = 1e-3;    % Max Δμ increase per timestep
delta_t_max            = 1e-2;    % Maximum allowed timestep
min_num_steps          = 5;       % Minimum steps per simulation

% Covariance matrix method parameters
use_parity_loss_or_subspace_leakage = 'parity_loss'; 
simulation_type = 'ramp_up_stop'; 
basis           = 'instant_hamil_basis'; 

%% Preallocate Results
gap_1_min_vec = zeros(length(N_vec), 1);
ramp_rate_vec = zeros(length(N_vec), 1);

BdG_leakage_mat = zeros(1, length(N_vec)); 
L_BdG_mat       = zeros(length(N_vec), 1);
L_BdG_comp_cell = cell(length(N_vec), 1);


%% Compute Leakage using Covariance Matrix Method
delta_t_cov_sim_vec = zeros(length(N_vec), 1); 
L_p_mat    = zeros(length(N_vec), 1); % Parity Leakage (Leakage out of even 
                                      % MZM parity sector
L_gs_mat   = zeros(length(N_vec), 1); % Ground state leakage (leakage out of 
                                      % the ground state manifold)
L_gs_mat_init = zeros(length(N_vec), 1); 
L_p0_mat   = zeros(length(N_vec), 1); 
L_c0_mat   = zeros(length(N_vec), 1);

for ii = 1:length(N_vec)

    fprintf('%i ', ii);
    if ~mod(ii, 10), fprintf('\n'); end

    N_curr         = N_vec(ii); 
    ramp_rate_curr = ramp_rate_constant;

    % Ramp time and timestep
    ramp_time = (ramp_height - mu_init) ./ ramp_rate_curr;
    t_final   = ramp_time + t_init;   

    delta_t_curr = min([ ...
        delta_mu_time_step_max ./ ramp_rate_curr, ...
        delta_t_max, ...
        t_final ./ min_num_steps ]);
    delta_t_curr = t_final ./ ceil(t_final ./ delta_t_curr);  
    delta_t_cov_sim_vec(ii) = delta_t_curr; 
    delta_t = delta_t_curr;

    % Time vector
    t_vec = t_init:delta_t:t_final;    

    % Define noise trajectory
    noise_top = (ramp_height - mu_init) .* t_vec ./ ramp_time;
    noise_traj_curr = [noise_top(:), noise_top(:)];

    % Run covariance matrix simulation
    init_condition = 'initial_state_no_noise';
    [omega_0_from_fit, T2_from_fit, omega_0_guess, X_exp_mean, ...
     P_y1y2y3y4, P_y1y2y3y4_v2, overlap_with_init_state_mean, ...
     overlap_fixed_basis_states_mean, t_vec, top_power_mat, bot_power_mat, ...
     mzm_parity_decay, comparison_determinants_mat, gap_mean, ...
     e_vals_mat, ~, overlap_inst_basis_states_mean] = ...
        charge_noise_engine_inst_basis_setup(mu_init, mu_offset, ...
        w, delta, N_curr, BC, t_init, t_final, delta_t, ...
        'custom', noise_traj_curr, num_trials, ...
        'time_ev_diagonalise', 'do_not_save_output', ...
        'do_not_plot_results', init_condition, basis); 

    % Extract leakage metrics
    final_parity      = P_y1y2y3y4(end); 
    L_p_mat(ii)       = 0.5 * (1 - final_parity); 
    L_gs_mat_init(ii) = 1 - overlap_fixed_basis_states_mean(end); 
    L_gs_mat(ii)      = 1 - overlap_inst_basis_states_mean(end);

    % Sudden quench limits
    P_0        = get_parity_quench_limit(mu_init, ramp_height, mu_offset, w, delta, N_curr, BC);
    L_p0_mat(ii) = 0.5 * (1 - P_0); 

    L_c0_mat(ii) = 1 - get_overlap_sq_quench_limit(mu_init, ramp_height, mu_offset, w, delta, N_curr, BC);
end

toc


%% Plot Results (Adiabatic)
L_odd  = L_p_mat;
L_even = L_gs_mat - L_odd;

% Plot settings
marker_size_c = 4;
line_width    = 2;
fig_border_LW = 2.3;
fig_ax_font_size = 22; 

col_1 = [0.9290 0.6940 0.1250];
col_2 = [0.8500 0.3250 0.0980];

legend_str   = {};
plot_handles = [];

currFig = figure();
plot_handles(end+1) = plot(N_vec, L_even, 'o', 'Color', col_1, ...
    'MarkerFaceColor', col_1, 'LineWidth', line_width, ...
    'MarkerSize', marker_size_c);
legend_str{end+1} = sprintf('L_{even}, \\mu_{fin} = %.3g', ramp_height);

hold on;
plot_handles(end+1) = plot(N_vec, L_odd, 'o', 'Color', col_2, ...
    'MarkerFaceColor', col_2, 'LineWidth', line_width, ...
    'MarkerSize', marker_size_c);
legend_str{end+1} = sprintf('L_{odd}, \\mu_{fin} = %.3g', ramp_height);

% Axis formatting
ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;

legend(plot_handles, legend_str, 'Location', 'NorthWest');
ylim('tight');
set(gca,'box','on')
currFig.Position = [300 251 672 466];

xlabel('Kitaev Chain Length, N');
ylabel('Leakage');

ylim([0, 1.2*max([L_even(:); L_odd(:)])]);


%% Save Figure
saveas(currFig, 'Results/fig_2_inset.png');    
saveas(currFig, 'Results/fig_2_inset.fig'); 


%% Save Workspace
clear currFig plot_handles ax;
save('Results/fig_2_inset_data.mat');
