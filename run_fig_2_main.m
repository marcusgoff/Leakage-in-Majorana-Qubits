%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: run_fig_2_main.m
% Purpose: Generate Figure 2 from the paper 
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
%----------------------------------------%
%% Clear
clear all; 
%close all;

%% Get QPP Library
addpath('QPP_Library_P1_submit');

%% Input Parameters
debug_plots_on = 0;
plot_straight_line_fits = 1;

simulation_type = 'ramp_up_stop'; % 'ramp_up_stop', 'ramp_up_hold_ramp_down'
basis = 'instant_hamil_basis';    % 'instant_hamil_basis', 'init_hamil_basis'

%--------------------------%
% Simulation Properties
%--------------------------%
initial_state = "+";

ramp_rate_vec = logspace(-4.5, 1, 500); 

ramp_height_vec = [0.03, 0.1, 0.01]; % Run 0.01 as well, as this data is accessed
                                     % in run_fig_S3.m in Supplementary
                                     % Plots.

% TETRON PARAMETERS
mu_init  = 0;   % Initial chemical potential
mu_offset = 0;  % Offset
w      = 1/2;  % Hopping
delta  = 1/2;  % Pairing
N      = 40;   % Number of sites in each Kitaev chain
BC     = "OBC"; % Boundary Conditions

% NOISE PARAMETERS 
t_init           = 0;
t_final_no_ramp  = 100; % Not used for ramp-up stop

% Time step controls
delta_mu_time_step_max = 1e-3; % Maximum Δμ per step
delta_t_max            = 1;    % Max Δt
min_num_steps          = 50;   % Minimum number of time steps

% MINIMUM RECOMMENDED VALUES:
% delta_mu_time_step_max = 1e-3; delta_t_max = 1; min_num_steps = 5;

num_trials = 1; % Keep fixed as 11

% Effective Gap Parameters
adiabatic_regime_low_bound = 3.8; % log10(g_1^2/v) lower bound
n_chosen = 1; % Chosen μ_final index for reference g_1

%% Ramped Noise Sweep
run_sim = 1;

if run_sim
    final_parity_mat              = zeros(length(ramp_rate_vec), length(ramp_height_vec));
    leakage_comp_mat              = zeros(length(ramp_rate_vec), length(ramp_height_vec));
    gap_1_min_vec                 = zeros(1, length(ramp_height_vec)); % Bulk gap
    gap_2_min_vec                 = zeros(1, length(ramp_height_vec)); % Excited gap
    delta_ramp_height_vec         = ramp_height_vec - mu_init;
    exp_N_site_basis_vec          = zeros(1, length(ramp_height_vec));
    ramp_rate_mat                 = zeros(length(ramp_rate_vec), length(ramp_height_vec));
    delta_t_mat                   = zeros(length(ramp_rate_vec), length(ramp_height_vec));
    num_points_mat                = zeros(length(ramp_rate_vec), length(ramp_height_vec));
    qp_population_final_mat       = zeros(length(ramp_rate_vec), length(ramp_height_vec));
    excited_qp_population_final_mat = zeros(length(ramp_rate_vec), length(ramp_height_vec));

    for ramp_height_ind = 1:length(ramp_height_vec)
        fprintf('ramp_height_ind = %i\n', ramp_height_vec(ramp_height_ind));
        
        % Compute <N> in site basis
        mu_curr = delta_ramp_height_vec(ramp_height_ind) + mu_init + mu_offset * [+1, -1];
        H_tetron_final = get_tetron_BdG_Hamiltonian(mu_curr, w, delta, N, BC, 'chain_1_chain_2');

        [e_vecs_final, e_vals_final, majorana_zero_modes_ref] = ...
            diagonalise_uncoupled_tetron_via_Kitaev_chains( ...
                H_tetron_final(1:2*N, 1:2*N), ...
                H_tetron_final((2*N+1):end, (2*N+1):end), ...
                'dirac_zero_modes');

        [~, corr_site_final, ~] = get_tetron_init_cov_mat_general(e_vecs_final, initial_state);
        temp_vec = diag(corr_site_final);
        exp_N_site_basis_vec(ramp_height_ind) = sum(temp_vec([1:N, (2*N+1):3*N]));

        ramp_rate_mat(:, ramp_height_ind) = ramp_rate_vec;

        for ramp_rate_ind = 1:length(ramp_rate_vec)
            fprintf('%i', ramp_rate_ind);
            if mod(ramp_rate_ind, 10) == 0, fprintf('\n'); end

            ramp_time_curr = delta_ramp_height_vec(ramp_height_ind) ./ ramp_rate_vec(ramp_rate_ind);
            delta_t_curr   = min([ ...
                delta_mu_time_step_max ./ ramp_rate_vec(ramp_rate_ind), ...
                delta_t_max, ...
                ramp_time_curr ./ min_num_steps]);

            delta_t_curr = ramp_time_curr ./ ceil(ramp_time_curr ./ delta_t_curr);
            num_points_mat(ramp_rate_ind, ramp_height_ind) = ceil(ramp_time_curr ./ delta_t_curr);

            % Ramp generation
            switch simulation_type
                case 'ramp_up_hold_ramp_down'
                    t_final_w_ramp = t_final_no_ramp + 2 * ramp_time_curr;
                    num_steps      = ceil(t_final_w_ramp ./ delta_t_curr);
                    delta_t_curr   = t_final_w_ramp ./ num_steps;
                    t_vec_curr     = t_init:delta_t_curr:t_final_w_ramp;
                    ramp           = get_ramp_envelope(t_init, t_final_w_ramp, delta_t_curr, ramp_time_curr);
                    noise_top      = delta_ramp_height_vec(ramp_height_ind) .* ramp;
                    noise_bot      = noise_top;

                case 'ramp_up_stop'
                    t_final_w_ramp = ramp_time_curr;
                    num_steps      = ceil(t_final_w_ramp ./ delta_t_curr);
                    delta_t_curr   = t_final_w_ramp ./ num_steps;
                    t_vec_curr     = t_init:delta_t_curr:t_final_w_ramp;
                    ramp           = t_vec_curr ./ ramp_time_curr;
                    noise_top      = delta_ramp_height_vec(ramp_height_ind) .* ramp;
                    noise_bot      = noise_top;

                    if norm(t_vec_curr(end) - ramp_time_curr) > 1e-10
                        warning('delta_t not commensurate with ramp time');
                    end

                otherwise
                    error('Invalid simulation_type');
            end

            delta_t_mat(ramp_rate_ind, ramp_height_ind) = delta_t_curr;
            noise_traj_curr = [noise_top(:), noise_bot(:)];

            if debug_plots_on
                figure();
                plot(t_vec_curr, noise_top);
                xlabel('t'); ylabel('\Delta \mu'); grid on;
            end

            % Run simulation
            t_final = t_final_w_ramp;
            init_condition = 'initial_state_no_noise';
            [omega_0_from_fit, T2_from_fit, omega_0_guess, X_exp_mean, ...
                P_y1y2y3y4, P_y1y2y3y4_v2, overlap_with_init_state_mean, ...
                overlap_fixed_basis_states_mean, t_vec, top_power_mat, ...
                bot_power_mat, mzm_parity_decay, comparison_determinants_mat, ...
                gap_mean, ~, quasi_occ_t, overlap_inst_basis_states_mean] = ...
                charge_noise_engine_inst_basis_setup( ...
                    mu_init, mu_offset, w, delta, N, BC, ...
                    t_init, t_final, delta_t_curr, 'custom', noise_traj_curr, ...
                    num_trials, 'time_ev_diagonalise', ...
                    'do_not_save_output', 'do_not_plot_results', ...
                    init_condition, basis, initial_state);

            final_parity_curr = P_y1y2y3y4(end);

            final_parity_mat(ramp_rate_ind, ramp_height_ind)  = final_parity_curr;
            leakage_comp_mat(ramp_rate_ind, ramp_height_ind) = 1 - overlap_inst_basis_states_mean(end);
            qp_population_final_mat(ramp_rate_ind, ramp_height_ind) = sum(quasi_occ_t(:, end));

            quasi_occ_final_excited_only = quasi_occ_t([1:N-1, (N+1):end-1], end);
            excited_qp_population_final_mat(ramp_rate_ind, ramp_height_ind) = sum(quasi_occ_final_excited_only);
        end
    end

    L_p_mat = 0.5 * (1 - final_parity_mat);
    L_c_mat = leakage_comp_mat;
end

%% Plot Parameters
fac = 0.8;

line_width     = 4; % 4 for legends, 6 for plots
fig_border_LW  = 3;
fig_ax_font_size = 40;

marker_size_s = 35 * fac;
marker_size_p = 40 * fac;
marker_LW     = 4;

% Colors
plot_colors_dark_all = [0,      0.4470, 0.7410; ...
                        0.8500, 0.3250, 0.0980; ...
                        0.1333, 0.5451, 0.1333];

plot_colors_light_all = [0.3010, 0.7450, 0.9330; ...
                         0.9290, 0.6940, 0.1250; ...
                         0.4670, 0.8670, 0.4670];

order = [2, 3, 1];
plot_colors_dark  = plot_colors_dark_all(order, :);
plot_colors_light = plot_colors_light_all(order, :);

currFig = figure();
legend_str   = {};
plot_handles = [];

% Selected curves
selected_curves = [1, 2];

L_gc_mat = leakage_comp_mat - L_p_mat;

for ramp_height_ind = selected_curves
    plot_handles(end+1) = loglog( ...
        ramp_rate_mat(:, ramp_height_ind), L_gc_mat(:, ramp_height_ind), '-', ...
        'Color', plot_colors_light(ramp_height_ind, :), 'LineWidth', line_width); hold on;

    legend_str{end+1} = sprintf('L_{even}, \\mu_{fin} = %.3g', ramp_height_vec(ramp_height_ind));
end

for ramp_height_ind = selected_curves
    plot_handles(end+1) = loglog( ...
        ramp_rate_mat(:, ramp_height_ind), L_p_mat(:, ramp_height_ind), '-.', ...
        'Color', plot_colors_dark(ramp_height_ind, :), 'LineWidth', line_width); hold on;

    legend_str{end+1} = sprintf('L_{odd}, \\mu_{fin} = %.3g', ramp_height_vec(ramp_height_ind));
end

xlabel('Ramp Rate, v');
ylabel('Leakage');
axis tight;

ax = gca;
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
ax.XAxis.MinorTick = 'off';
ax.YAxis.MinorTick = 'off';

currFig.Position = [156, 164, 1109, 641];

yticks([1e-8, 1e-5, 1e-2]);
xticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1]);

xlim([0, 12.8]);
ylim([0, 0.1]);

legend_order = [1, 3, 2, 4];
legend(plot_handles(legend_order), legend_str(legend_order), 'Location', 'SouthEast');

%% Save Figure
saveas(currFig, 'Results/fig_2_main.png');    
saveas(currFig, 'Results/fig_2_main.fig'); 


%% Save Workspace
clear currFig plot_handles ax;
save('Results/fig_2_main_data.mat');


% - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Helper Function
function ramp = get_ramp_envelope(t_init, t_final, delta_t, ramp_time)
    t_vec = t_init:delta_t:t_final;
    ramp  = zeros(length(t_vec), 1);

    t_ramp_down = t_vec(end) - ramp_time;

    ind_ramp_up   = (t_vec < ramp_time);
    ind_ramp_down = (t_vec > t_ramp_down);

    ramp(ind_ramp_up)   = t_vec(ind_ramp_up) ./ ramp_time;
    ramp(ind_ramp_down) = 1 - (t_vec(ind_ramp_down) - t_ramp_down) ./ ramp_time;
    ramp(~ind_ramp_up & ~ind_ramp_down) = 1;
end
