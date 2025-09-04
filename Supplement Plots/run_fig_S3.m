%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: run_fig_S3.m
%
% Author: Marcus C. Goffage
% Affiliation: University of New South Wales
%
% Purpose:
%   Generate Figure S3 from the supplementary material of the paper:
%   "Leakage at Zero Temperature from Changes in Chemical Potential
%    in Majorana Qubits"
%
% Author: Marcus C. Goffage
% Affiliation: University of New South Wales
%
% Paper Reference:
%   "Decoherence in Majorana Qubits by 1/f Noise"
%   Authors: M. C. Goffage^1, A. Alase^2, M. C. Cassidy^1, S. N. Coppersmith^1
%   Affiliations: ^1 University of New South Wales, ^2 University of Sydney
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load QPP Library
addpath('../QPP_Library_P1_submit');

%% Load Precomputed Data
try
    load('../Results/fig_2_main_data.mat');
catch ME
    error(['You must execute run_fig_2_main.m prior to run_fig_S3.m. ' ...
           'Could not load results file: %s\nOriginal error: %s'], ...
           '../Results/fig_2_main_data.mat', ME.message);
end

% Upper bound of adiabatic regime (due to flipped x_1 → -x_1)
adiabatic_regime_up_bound = -2; 

%% Plotting Parameters
% Define color schemes (dark and light variants)
plot_colors_dark_all = [0.0000 0.4470 0.7410;  ... % blue
                        0.8500 0.3250 0.0980;  ... % orange
                        0.1333 0.5451 0.1333];     % green
plot_colors_light_all = [0.3010 0.7450 0.9330;  ... % light blue
                         0.9290 0.6940 0.1250;  ... % yellow
                         0.4670 0.8670 0.4670];    % light green

order            = [2 3 1]; % re-ordering for consistency
plot_colors_dark = plot_colors_dark_all(order, :);
plot_colors_light= plot_colors_light_all(order, :);

% Figure style parameters
line_width      = 6;   % use 6 for plots, 4 for legend
fig_border_LW   = 3;
fig_ax_font_size= 35;

%% Compute Limiting Values
L_p0_vec = zeros(1, length(ramp_height_vec));
L_c0_vec = zeros(1, length(ramp_height_vec));

for ramp_height_ind = 1:length(ramp_height_vec)
    % Odd-parity leakage limit
    P_0 = get_parity_quench_limit(mu_init, ...
        ramp_height_vec(ramp_height_ind), mu_offset, w, delta, N, BC); 
    L_p0_vec(ramp_height_ind) = 0.5 * (1 - P_0);

    % Total leakage limit
    L_c0_vec(ramp_height_ind) = 1 - get_overlap_sq_quench_limit(mu_init, ...
        ramp_height_vec(ramp_height_ind), mu_offset, w, delta, N, BC);
end

% Even-parity leakage limit
L_gc0_vec = L_c0_vec - L_p0_vec;

%% Plot: Approaching Limiting Values (L_p0 and L_c0)
currFig    = figure();
legend_str = {};

for ramp_height_ind = 1:length(ramp_height_vec)
    % Extract limiting values
    L_p0  = L_p0_vec(ramp_height_ind);
    L_c0  = L_c0_vec(ramp_height_ind);
    L_gc0 = L_gc0_vec(ramp_height_ind);

    % Extract current values
    L_p_curr = L_p_mat(:, ramp_height_ind);
    L_c_curr = leakage_comp_mat(:, ramp_height_ind);
    L_gc_curr= L_c_curr - L_p_curr;

    % Compute scaled differences
    L_p_diff_scaled  = (L_p0  - L_p_curr) ./ L_p0;
    L_c_diff_scaled  = (L_c0  - L_c_curr) ./ L_c0;
    L_gc_diff_scaled = (L_gc0 - L_gc_curr) ./ L_gc0;

    % Ramp rate values
    v_curr = ramp_rate_mat(:, ramp_height_ind); 

    % Plot odd-parity leakage
    loglog(v_curr, L_p_diff_scaled, '-.', ...
        'Color', plot_colors_dark(ramp_height_ind, :), ...
        'LineWidth', line_width); 
    hold on;
    legend_str{end+1} = sprintf('L_{odd}, \\mu_{fin} = %.3g', ...
        ramp_height_vec(ramp_height_ind));

    % Plot even-parity leakage
    loglog(v_curr, L_gc_diff_scaled, '-', ...
        'Color', plot_colors_light(ramp_height_ind, :), ...
        'LineWidth', line_width);
    hold on;
    legend_str{end+1} = sprintf('L_{even}, \\mu_{fin} = %.3g', ...
        ramp_height_vec(ramp_height_ind));   
end

%% Final Figure Formatting
legend(legend_str, 'Location', 'southwest');
grid off;
axis tight;

ax = gca;
ax.FontSize       = fig_ax_font_size;
ax.LineWidth      = fig_border_LW;
ax.XAxis.MinorTick= 'off';  
ax.YAxis.MinorTick= 'off'; 
ylim([2e-5, 1e0]);
yticks([1e-4, 1e-2, 1e0]);

% Figure size
currFig.Position = [161 151 1230 639];

% Axis labels
xlabel('Ramp Rate, v'); 
ylabel('(L_{x}(\infty) - L_x)/L_{x}(\infty)');

%% Save

saveas(currFig, 'Results/fig_S3.png');    
saveas(currFig, 'Results/fig_S3.fig'); 
clear currFig plot_handles ax;


%% Save Worskpace
save('Results/fig_S3_data.mat');
