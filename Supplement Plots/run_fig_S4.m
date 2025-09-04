%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: run_fig_S4.m
% Purpose:
%   Generate Figure S4 from the supplementary material of the paper:
%   "Leakage at Zero Temperature from Changes in Chemical Potential
%    in Majorana Qubits"
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
%   - Generates Figure S4 
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves figures ./results
%   - Saves entire worskpace in the ./results directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running run_fig_S4.m');

%% Environment Setup
clear all; 
%close all;
addpath('../QPP_Library_P1_submit')


%% Global Parameters
mu_init  = 0;
mu_offset = 0;
w       = 1/2;
delta   = 1/2;
N       = 40;
BC      = "OBC";

ramp_rate_vec   = logspace(-4, -3, 500); 

%% Top Row mu = 0.03
ramp_height_vec_a = [0.03];

[L_p_mat_a, L_c_mat_a] = run_tetron_ramp_local(ramp_rate_vec, ramp_height_vec_a,...
    mu_init, mu_offset, w, delta, N, BC);

[fig_odd, fig_even, fig_loglog, odd_fit_params, even_fit_params] = ...
         plot_S4_panels(L_c_mat_a, L_p_mat_a, ramp_rate_vec.', ramp_height_vec_a, N);

saveas(fig_loglog, 'Results/fig_S4_a.png');    
saveas(fig_loglog, 'Results/fig_S4_a.fig'); 
clear fig_loglog plot_handles ax;

saveas(fig_odd, 'Results/fig_S4_b.png');    
saveas(fig_odd, 'Results/fig_S4_b.fig'); 
clear fig_odd plot_handles ax;

saveas(fig_even, 'Results/fig_S4_c.png');    
saveas(fig_even, 'Results/fig_S4_c.fig');
clear fig_even plot_handles ax;

%% Bottom Row mu = 0.01
ramp_height_vec_b = [0.01];

[L_p_mat_b, L_c_mat_b] = run_tetron_ramp_local(ramp_rate_vec, ramp_height_vec_b,...
    mu_init, mu_offset, w, delta, N, BC);

[fig_odd, fig_even, fig_loglog, odd_fit_params, even_fit_params] = ...
         plot_S4_panels(L_c_mat_b, L_p_mat_b, ramp_rate_vec.', ramp_height_vec_b, N);

saveas(fig_loglog, 'Results/fig_S4_d.png');    
saveas(fig_loglog, 'Results/fig_S4_d.fig'); 
clear fig_loglog plot_handles ax;

saveas(fig_odd, 'Results/fig_S4_e.png');    
saveas(fig_odd, 'Results/fig_S4_e.fig'); 
clear fig_odd plot_handles ax;

saveas(fig_even, 'Results/fig_S4_f.png');    
saveas(fig_even, 'Results/fig_S4_f.fig');
clear fig_even plot_handles ax;


%% Save Worskpace
save('Results/fig_S4_data.mat');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Local Functions
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [fig_odd, fig_even, fig_loglog, odd_fit_params, even_fit_params] = ...
         plot_S4_panels(L_c_mat, L_p_mat, ramp_rate_mat, ramp_height_vec, N)
% Wraps original script into a function
% Returns handles to the three generated figures and the fit parameters

%% Original Script Starts

up_bound = -2; % because I've flipped x_1 -> -x_1.
low_bound = -4; % log10(4e-4); % log10(4e-4);

fit_dom_low = log10(4e-4);

%% Control Plot Params

plot_colors_dark_all = [0 0.4470 0.7410; ...
                    0.8500 0.3250 0.0980;...
                    0.1333 0.5451 0.1333];

plot_colors_light_all = [0.3010 0.7450 0.9330; ...
                    0.9290 0.6940 0.1250; ...
                    0.467 0.867 0.467];

col_grey = 0.5*[1 1 1];

order = [1 2 3];
plot_colors_dark = plot_colors_dark_all(order, :);
plot_colors_light = plot_colors_light_all(order, :);

legend_str = {};
plot_handles = [];

line_width = 6; %6
fig_border_LW = 3;
fig_ax_font_size = 35; 

show_fit = 1; 

%% Function Fit - L_odd

L_even_mat = L_c_mat - L_p_mat;
L_odd_mat = L_p_mat;

legend_str = {};
plot_handles = [];

selected_curves = [1];
for ramp_height_ind = selected_curves  

    v = ramp_rate_mat(:, ramp_height_ind);
    v_plot_ind =  log10(v) <= up_bound &...
                  log10(v) >= low_bound;
    v_temp = v(v_plot_ind);
    L_temp = L_odd_mat(v_plot_ind, ramp_height_ind);

    v_fit_ind = log10(v) <= up_bound &...
                  log10(v) >= fit_dom_low;
    x = 1./v(v_fit_ind);
    
    [y_mid, y_diff, top_envelope, bot_envelope, min_mat, max_mat,...
        fit_top, fit_bot, mean_peak_diff, fit_mid, fit_diff, max_ind_2] = ...
        fit_top_and_bottom(x, L_odd_mat(v_fit_ind, ramp_height_ind));

    k1 = 0.5*(10^fit_top(2) +  10^fit_bot(2));
    k2_v2 = 10.^(fit_mid(2));

    k2 = 0.5*(10^fit_top(2) -  10^fit_bot(2));

    fig_odd = figure();
    plot_handles(end+1) = plot(v_temp, L_temp,...
         'Color', plot_colors_dark(2, :), 'LineWidth', line_width); hold on;

    m1 = -fit_mid(1);
    m2 = -fit_diff(1);
    
    k1 = 10.^(fit_mid(2));
    k2 = 10.^(fit_diff(2));
    k3 = (2*pi)./mean_peak_diff; 

    fit_vec = k1.*v_temp.^m1 - k2.*v_temp.^(m2).*cos(k3./v_temp);

    if show_fit
        plot_handles(end+1) = plot(v_temp, fit_vec,...
                    '-.', 'Color', plot_colors_dark(1, :), 'LineWidth', 0.7*line_width); hold on;
        legend(sprintf('$L_{\\rm odd}, \\mu_{\\rm fin} = %.3g$', ramp_height_vec), ...
              sprintf('Fit: $%.3g v^{%.3g} - %.3g v^{%.3g} {\\rm cos} (%.2g/v)$', ...
                k1, m1, k2, m2, abs(k3)), 'Interpreter', 'latex', 'Location', 'northwest');       
    end
    xlabel('Ramp Rate, v'); ylabel('L_{odd}');  
     
    ax = gca;
    ax.YAxis.Exponent = 0;
    ax.FontSize = fig_ax_font_size;
    ax.LineWidth = fig_border_LW;
    fig_odd.Position = [254 325 989 420];
    ytickformat('%.0e');

    ax.XAxis.Exponent = -4;
    ylim([0, 6e-6]);
end

odd_fit_params = [k1 k2 k3 m1 m2];

%% Function Fit - L_even

L_even_mat = L_c_mat - L_p_mat;
L_odd_mat = L_p_mat;

legend_str = {};
plot_handles = [];

selected_curves = [1];
for ramp_height_ind = selected_curves  

    v = ramp_rate_mat(:, ramp_height_ind);
    v_plot_ind =  log10(v) <= up_bound &...
                  log10(v) >= low_bound;
    v_temp = v(v_plot_ind);
    L_temp = L_even_mat(v_plot_ind, ramp_height_ind);

    v_fit_ind = log10(v) <= up_bound &...
                  log10(v) >= fit_dom_low;
    x = 1./v(v_fit_ind);
    
    [y_mid, y_diff, top_envelope, bot_envelope, min_mat, max_mat,...
        fit_top, fit_bot, mean_peak_diff, fit_mid, fit_diff, max_ind_2] = ...
        fit_top_and_bottom(x, L_even_mat(v_fit_ind, ramp_height_ind));

    k1 = 0.5*(10^fit_top(2) +  10^fit_bot(2));
    k2_v2 = 10.^(fit_mid(2));

    k2 = 0.5*(10^fit_top(2) -  10^fit_bot(2));

    fig_even = figure();
    plot_handles(end+1) = plot(v_temp, L_temp,...
         'Color', plot_colors_light(2, :), 'LineWidth', line_width); hold on;

    m1 = -fit_mid(1);
    m2 = -fit_diff(1);
    
    k1 = 10.^(fit_mid(2));
    k2 = 10.^(fit_diff(2));
    k3 = (2*pi)./mean_peak_diff; 

    fit_vec = k1.*v_temp.^m1 - k2.*v_temp.^(m2).*cos(k3./v_temp);

    if show_fit 
        plot_handles(end+1) = plot(v_temp, fit_vec,...
                   '-.', 'Color', plot_colors_dark(1, :),'LineWidth', 0.7*line_width); hold on;
    end
    xlabel('Ramp Rate, v'); ylabel('L_{even}');  
     
    ax = gca;
    ax.YAxis.Exponent = 0;
    ax.FontSize = fig_ax_font_size;
    ax.LineWidth = fig_border_LW;
    fig_even.Position = [254 325 989 420];
    ytickformat('%.0e');
    ax.XAxis.Exponent = -4;
    ylim([0, 6e-6]);

    hold on;
    plot(v_temp, (N).*v_temp.^2./8, '-.', 'Color', col_grey, 'LineWidth', 0.7*line_width);

    if show_fit
        legend(sprintf('$L_{\\rm even}, \\mu_{\\rm fin} = %.3g$', ramp_height_vec), ...
              sprintf('Fit: $%.3g v^{%.3g} - %.3g v^{%.3g} {\\rm cos} (%.2g/v)$', k1, m1, k2, m2, abs(k3)),...
              sprintf('$\\tilde{L}_{\\rm even, near-ad}, \\mu_{\\rm fin} = %.3g$', ramp_height_vec),...
              'Interpreter', 'latex', 'Location', 'northwest');
    end
end

even_fit_params = [k1 k2 k3 m1 m2];

%% LogLog Plot

fig_loglog = figure();
plot_handles(end+1) = loglog(v(v_plot_ind), L_even_mat(v_plot_ind, ramp_height_ind),...
     'Color', plot_colors_light(2, :), 'LineWidth', line_width); hold on;
plot_handles(end+1) = loglog(v(v_plot_ind), L_odd_mat(v_plot_ind, ramp_height_ind), '-.',...
     'Color', plot_colors_dark(2, :), 'LineWidth', line_width); hold on;
plot_handles(end+1) = loglog(v(v_plot_ind), (N).*v(v_plot_ind).^2./8, '-.', 'Color', col_grey,...
    'LineWidth', 0.7*line_width);

ax = gca;
ax.YAxis.Exponent = 0;
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
axis tight;

xlabel('Ramp Rate, v');
ylabel('Leakage');
fig_loglog.Position = [251   308   785   437];

%% Local Function

function [y_mid, y_diff, top_envelope, bot_envelope, min_mat, max_mat,...
    fit_top, fit_bot, mean_peak_diff, fit_mid, fit_diff, max_ind_2] =  fit_top_and_bottom(x, y)

    min_ind = islocalmin(y);
    max_ind = islocalmax(y);
    
    y_max = y(max_ind);
    x_max = x(max_ind);
    max_mat = [x_max(:), y_max(:)];
    
    y_min = y(min_ind);
    x_min = x(min_ind); 
    min_mat = [x_min(:), y_min(:)];
    
    fit_top = polyfit(log10(x_max), log10(y_max), 1);
    fit_bot = polyfit(log10(x_min), log10(y_min), 1);
    
    top_envelope = (10^fit_top(2).*x.^fit_top(1));
    bot_envelope = (10^fit_bot(2).*x.^fit_bot(1));
    
    y_mid = (top_envelope + bot_envelope)./2;
    y_diff = (top_envelope - bot_envelope)./2;

    fit_mid = polyfit(log10(x), log10(y_mid), 1);
    fit_diff = polyfit(log10(x), log10(y_diff), 1);

    max_ind_2 = islocalmax(y - y_mid);
    mean_peak_diff = mean(diff(x(max_ind_2))); 
end

end




function [L_p_mat, L_c_mat] = run_tetron_ramp_local(ramp_rate_vec, ramp_height_vec, mu_init, mu_offset, w, delta, N, BC)
%% Hardcoded Parameters
debug_plots_on = 0;
plot_straight_line_fits = 1;

simulation_type = 'ramp_up_stop'; % 'ramp_up_stop', 'ramp_up_hold_ramp_down'
basis = 'instant_hamil_basis';    % 'instant_hamil_basis', 'init_hamil_basis'

initial_state = "+";

% NOISE PARAMETERS 
t_init           = 0;
t_final_no_ramp  = 100; % Not used for ramp-up stop

% Time step controls
delta_mu_time_step_max = 1e-3; % Maximum Δμ per step
delta_t_max            = 1;    % Max Δt
min_num_steps          = 50;   % Minimum number of time steps

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
        end
    end

    L_p_mat = 0.5 * (1 - final_parity_mat);
    L_c_mat = leakage_comp_mat;
end

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
end

