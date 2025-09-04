%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: run_fig_S2.m
% Purpose: Generate Figure S2 from the supplementary material of the paper
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

%% Setup
clear all;
%close all;

%% Get QPP Library
addpath('../QPP_Library_P1_submit');


%% Global Parameters
% Tetron Parameters
mu_init = 0;
mu_offset = 0;
w = 1/2;
delta = 1/2;
BC = "OBC";

N_vec = 2:1:100;

%% Subplot A
ramp_height = 0.01;
[~, L_even_a, L_odd_a, L_odd_approx_a, L_even_approx_a, fig_a] = ...
    figure_S3_subplot_func(mu_init, mu_offset, w, delta, BC, N_vec, ramp_height);

saveas(fig_a, 'Results/fig_S2_a.png');    
saveas(fig_a, 'Results/fig_S2_a.fig'); 
clear fig_a plot_handles ax;
%% Subplot B
ramp_height = 0.03;
[~, L_even_b, L_odd_b, L_odd_approx_b, L_even_approx_b, fig_b] = ...
    figure_S3_subplot_func(mu_init, mu_offset, w, delta, BC, N_vec, ramp_height);

saveas(fig_b, 'Results/fig_S2_b.png');    
saveas(fig_b, 'Results/fig_S2_b.fig'); 
clear fig_b plot_handles ax;
%% Subplot C
ramp_height = 0.1;
[~, L_even_c, L_odd_c, L_odd_approx_c, L_even_approx_c, fig_c] = ...
    figure_S3_subplot_func(mu_init, mu_offset, w, delta, BC, N_vec, ramp_height);

saveas(fig_c, 'Results/fig_S2_c.png');    
saveas(fig_c, 'Results/fig_S2_c.fig'); 
clear fig_c plot_handles ax;
%% Subplot D
ramp_height = 0.5;
[~, L_even_d, L_odd_d, L_odd_approx_d, L_even_approx_d, fig_d] = ...
    figure_S3_subplot_func(mu_init, mu_offset, w, delta, BC, N_vec, ramp_height);

saveas(fig_d, 'Results/fig_S2_d.png');    
saveas(fig_d, 'Results/fig_S2_d.fig'); 
clear fig_d plot_handles ax;


%% Save Worskpace
save('Results/fig_S2_data.mat');





%% Local Functions

function [N_vec, L_even_mat, L_odd_mat, L_odd_approx_mat, L_even_approx, currFig] = ...
    figure_S3_subplot_func(mu_init, mu_offset, w, delta, BC, N_vec, ramp_height)
%FIGURE_S3_SUBPLOT_FUNC 
%   Generates leakage plots for the tetron model.
%
%   Inputs:
%       mu_init     - Initial chemical potential
%       mu_offset   - Offset chemical potential
%       w           - Coupling parameter
%       delta       - Superconducting gap parameter
%       BC          - Boundary condition ("OBC" or "PBC")
%       N_vec       - Vector of chain lengths
%       ramp_height - Ramp height
%
%   Outputs:
%       N_vec             - Vector of chain lengths
%       L_even_mat        - Even-parity leakage (numerical)
%       L_odd_mat         - Odd-parity leakage (numerical)
%       L_odd_approx_mat  - Odd-parity leakage (approximate)
%       L_even_approx     - Even-parity leakage (approximate)
%       currFig           - Handle to the generated figure

    %% Hard-coded simulation parameters
    num_trials             = 1; 
    delta_mu_time_step_max = 1e-3;
    delta_t_max            = 1e-2; 
    min_num_steps          = 5; 

    use_parity_loss_or_subspace_leakage = 'parity_loss'; 
    simulation_type = 'ramp_up_stop'; %#ok<NASGU>
    basis = 'instant_hamil_basis'; %#ok<NASGU>

    %% Calculate leakages using covariance matrices
    L_gs_mat        = zeros(length(N_vec), 1); 
    L_odd_mat       = zeros(length(N_vec), 1); 
    P_0_approx_mat  = zeros(length(N_vec), 1); 
    L_odd_approx_mat= zeros(length(N_vec), 1); 

    for ii = 1:length(N_vec)
        fprintf('%i ', ii);
        if mod(ii,10) == 0, fprintf('\n'); end

        % Sudden limit parameters
        [P_0, P_0_approx] = get_parity_quench_limit_with_approx( ...
            mu_init, ramp_height, mu_offset, w, delta, N_vec(ii), BC);

        L_odd_mat(ii) = 0.5 * (1 - P_0); 

        L_gs_mat(ii) = 1 - get_overlap_sq_quench_limit( ...
            mu_init, ramp_height, mu_offset, w, delta, N_vec(ii), BC);

        P_0_approx_mat(ii)  = P_0_approx;
        L_odd_approx_mat(ii)= 0.5 * (1 - P_0_approx); 
    end
    fprintf('\n');
    L_even_mat = L_gs_mat - L_odd_mat;

    %% Approximation curves
    N_vec_cont    = linspace(min(N_vec), max(N_vec)); 
    L_even_approx = ((N_vec_cont - 2) .* ramp_height.^2 / 8);

    %% Plot style parameters
    font_size_scale   = 0.8; 

    marker_size_c     = 30;
    marker_size_cross = 6;
    marker_LW         = 4;

    line_width        = 1.5;
    fig_border_LW     = 2.3;
    fig_ax_font_size  = 45; 
    legend_font_size  = 20;

    %% Colors
    col_1     = [0.9290 0.6940 0.1250];
    col_2     = [0.8500 0.3250 0.0980];
    col_3     = [0 0 0];
    col_analy = 0.3 * ones(3,1);

    %% Plotting
    legend_str   = {};
    plot_handles = [];

    currFig = figure();

    % Analytical L_even approx
    plot_handles(end+1) = plot(N_vec_cont, L_even_approx, '-', ...
        'Color', col_analy, 'MarkerFaceColor', col_analy, ...
        'LineWidth', 8*line_width, 'MarkerSize', marker_size_cross);
    legend_str{end+1} = sprintf('$\\tilde{L}_{\\rm even}, \\mu_{\\rm fin} = %.3g, v = \\infty$', ramp_height);
    hold on;

    % Numerical L_even
    plot_handles(end+1) = plot(N_vec, L_even_mat, '.', ...
        'Color', col_1, 'MarkerFaceColor', col_1, ...
        'LineWidth', line_width, 'MarkerSize', marker_size_c);
    legend_str{end+1} = sprintf('$L_{\\rm even}, \\mu_{\\rm fin} = %.3g, v = \\infty$', ramp_height);

    % Numerical L_odd
    plot_handles(end+1) = plot(N_vec, L_odd_mat, '.', ...
        'Color', col_2, 'MarkerFaceColor', col_2, ...
        'LineWidth', line_width, 'MarkerSize', marker_size_c);
    legend_str{end+1} = sprintf('$L_{\\rm odd}, \\mu_{\\rm fin} = %.3g, v = \\infty$', ramp_height);

    % Analytical L_odd approx
    plot_handles(end+1) = plot(N_vec, L_odd_approx_mat, '*', ...
        'Color', col_3, 'MarkerFaceColor', col_3, ...
        'LineWidth', line_width, 'MarkerSize', marker_size_cross);
    legend_str{end+1} = sprintf('$\\tilde{L}_{\\rm odd}, \\mu_{\\rm fin} = %.3g, v = \\infty$', ramp_height);

    %% Axis formatting
    ax = gca;
    ax.YAxis.Exponent = 0;
    ax.FontSize       = fig_ax_font_size;
    ax.LineWidth      = fig_border_LW;
    set(gca,'box','on')

    currFig.Position = [125 172 1070 645];

    ylabel('Leakage');
    xlabel('Kitaev Chain Length, N');

    %% Legend
    legend_order = [2,3,1,4];
    h = legend(plot_handles(legend_order), legend_str(legend_order), 'Location', 'east');
    set(h,'interpreter','Latex','FontSize',legend_font_size);

end
