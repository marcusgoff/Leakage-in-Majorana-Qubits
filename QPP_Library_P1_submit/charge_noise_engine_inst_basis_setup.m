%--------------------------------------%
% Charge Noise Engine - General
%--------------------------------------%
% Author: Marcus C. Goffage
% Affiliation: University of New South Wales
%
% General version of the charge noise engine for different inbuilt noise
% types OR with your own imported noise trajectory.
%
% INPUTS:
%   noise_type        : 'gaussian', 'two_state_exp', or 'custom'
%   noise_parameters  : Array or matrix, depending on noise_type:
%                       Gaussian: [t_noise_factor, fluctuation_amplitude, 
%                                 characteristic_noise_corr_freq]
%                       two_state_exp: [sigma, tau, fluctuation_amplitude]
%                       custom: matrix of size "num_time_steps x 2*num_trials",
%                               columns contain independent noise trajectories
%                               (first num_trials = top chain, next num_trials = bottom)
%   time_evolution_method : 'time_ev_diagonalise', 'time_ev_exp', 
%                           'time_ev_1st_order', 'time_ev_Magnus' (TO BE IMPLEMENTED)
%   save_output       : 'save_output' or 'do_not_save_output'
%   plot_results      : 'plot_all_results', 'plot_X_exp_only', or 'do_not_plot_results'
%   noise_trajectories: optional (varargin{1}) 2*num_trials x num_time_steps_noise
%                       matrix of noise realizations.
%   init_condition    : 'initial_state_no_noise', 'initial_state_noise'
%                       (not included in this version)
%   basis_choice      : 'init_hamil_basis', 'instant_hamil_basis'
%   delta_t_vec       : scalar or vector of delta_t values
%
% OUTPUTS:
%   e_vals_mat                   : eigenvalues for first trajectory
%   omega_0_from_fit, T2_from_fit: fit results from X_exp
%   omega_0_guess                : estimate of initial frequency
%   X_exp_mean, Z_exp_mean       : mean expectation values
%   P_y1y2y3y4_v1_mean, P_y1y2y3y4_v2_mean : parity expectation values
%   overlap_with_init_state_mean, overlap_fixed_basis_states_mean
%   overlap_inst_basis_states_mean
%   t_vec                        : time vector
%   top_power_mat, bot_power_mat : total noise power per trial
%   mzm_parity_decay             : decay time of total MZM parity
%   comparison_determinants_cell : cell containing comparison determinants
%   gap_mean                     : mean energy gap
%   quasi_occ_t                  : quasi-particle occupation evolution

function [omega_0_from_fit, T2_from_fit, omega_0_guess, X_exp_mean,...
          P_y1y2y3y4_v1_mean, P_y1y2y3y4_v2_mean, overlap_with_init_state_mean,...
          overlap_fixed_basis_states_mean, t_vec, top_power_mat, bot_power_mat,...
          mzm_parity_decay, comparison_determinants_cell, gap_mean, e_vals_mat,...
          quasi_occ_t, overlap_inst_basis_states_mean, Z_exp_mean] = ...
    charge_noise_engine_inst_basis_setup(mu_mean, mu_offset, w, delta, N, BC, ...
        t_init, t_final, delta_t, noise_type, noise_parameters, num_trials, ...
        time_evolution_method, save_output, plot_results, init_condition, ...
        basis_choice, varargin)

    %--------------------------------------%
    % Check valid inputs
    %--------------------------------------%
    check_valid_inputs(save_output, plot_results);

    %--------------------------------------%
    % Optional input processing (varargin)
    %--------------------------------------%
    if nargin == 17
        init_state = "+";
    elseif nargin == 18
        init_state = varargin{1};
    elseif nargin > 19
        error('Too many input parameters');
    end

    if length(delta_t) > 1 && ~strcmp(noise_type, 'custom')
        error('When delta_t is a vector, noise_type must be custom');
    end

    %--------------------------------------%
    % Process Input Parameters
    %--------------------------------------%
    mu = mu_mean + mu_offset * [+1, -1];

    %--------------------------------------%
    % Get Initial Hamiltonian (if no noise)
    %--------------------------------------%
    if strcmp(init_condition, 'initial_state_no_noise')
        H_tetron_init = get_tetron_BdG_Hamiltonian(mu, w, delta, N, BC, 'chain_1_chain_2');
        [e_vecs_init, e_vals_init, majorana_zero_modes_ref] = ...
            diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
            H_tetron_init((2*N+1):end,(2*N+1):end), 'dirac_zero_modes');
        
        [~, ~, corr_qp_init] = get_tetron_init_cov_mat_general(e_vecs_init, init_state);
    end

    %--------------------------------------%
    % Time vector and initialization
    %--------------------------------------%
    if isscalar(delta_t)
        t_vec = t_init:delta_t:t_final;
    else
        t_vec = t_init + [t_init, cumsum(delta_t(:)).'];
        if abs((t_vec(end) - t_final)/t_final) > 10*eps
            warning('delta_t vector does not match t_init and t_final');
            t_final = t_vec(end);
        end
    end

    % Initialize matrices to store results
    Z_exp_mat = zeros(num_trials, length(t_vec));
    X_exp_mat = zeros(num_trials, length(t_vec));
    Y_exp_mat = zeros(num_trials, length(t_vec));
    P_y1y2y3y4_v1_mat = zeros(num_trials, length(t_vec));
    P_y1y2y3y4_v2_mat = zeros(num_trials, length(t_vec));
    mu_mat_top = zeros(num_trials, length(t_vec));
    mu_mat_bot = zeros(num_trials, length(t_vec));
    top_power_mat = zeros(num_trials,1);
    bot_power_mat = zeros(num_trials,1);
    comparison_determinants_mat_top = zeros(num_trials, length(t_vec));
    comparison_determinants_mat_bot = zeros(num_trials, length(t_vec));
    overlap_with_init_state_mat = zeros(num_trials, length(t_vec));
    overlap_fixed_basis_states_mat = zeros(num_trials, length(t_vec));
    overlap_inst_basis_states_mat = zeros(num_trials, length(t_vec));

    tStart = tic;

    if num_trials > 1
        fprintf('Trial Number: \n');
    end

    %--------------------------------------%
    % Loop over trials
    %--------------------------------------%
    for trial_no = 1:num_trials
        if trial_no > 1
            fprintf('%i ', trial_no);
            if mod(trial_no, 10) == 0
                fprintf('\n');
            end
        end

        %--------------------------------------%
        % Noise trajectories
        %--------------------------------------%
        if ~strcmp(noise_type, 'custom')
            error(['In-house noise generation not included in this code version.\n', ...
                   'Please proceed with a custom noise trajectory.']);
        else
            noise_top = noise_parameters(:, trial_no).';
            noise_bot = noise_parameters(:, trial_no + num_trials).';
            top_power = sqrt(sum(abs(noise_top).^2) * (t_vec(2) - t_vec(1)));
            bot_power = sqrt(sum(abs(noise_bot).^2) * (t_vec(2) - t_vec(1)));
            if ~isscalar(delta_t)
                warning('Variable delta_t gives incorrect powers');
            end
            fluctuation_amplitude = sqrt(top_power);
        end

        top_power_mat(trial_no) = top_power;
        bot_power_mat(trial_no) = bot_power;

        % Construct chemical potential vector for top and bottom chains
        mu_vec = zeros(2, length(t_vec));
        mu_vec(1,:) = mu_mean + mu_offset + noise_top;
        mu_vec(2,:) = mu_mean - mu_offset + noise_bot;
        mu_mat_top(trial_no, :) = mu_vec(1,:);
        mu_mat_bot(trial_no, :) = mu_vec(2,:);

        %--------------------------------------%
        % Initial Hamiltonian for noisy start (not implemented)
        %--------------------------------------%
        if strcmp(init_condition, 'initial_state_noise')
            error('initial_state_noise functionality not included in this version');
        end

        gap_mean = sort(abs(e_vals_init));
        gap_mean = gap_mean(5); % omit zero modes

        %--------------------------------------%
        % Time evolution calculation
        %--------------------------------------%
        tetron_params = {w, delta, N, BC};

        if nargout >= 17
            [Z_exp, X_exp, Y_exp, P_y1y2y3y4_v1, P_y1y2y3y4_v2, ...
             overlap_with_init_state, overlap_fixed_basis_states, ~,~,~, ...
             comparison_determinants_mat, e_vals_mat_curr, ~, quasi_occ_t, ...
             overlap_inst_basis_states] = ...
                calc_noisy_tetron_ev_inst_basis_setup(corr_qp_init, e_vecs_init, ...
                e_vals_init, tetron_params, mu_vec, delta_t, t_init, t_final, ...
                time_evolution_method, basis_choice, majorana_zero_modes_ref, mu);
        else
            [Z_exp, X_exp, Y_exp, P_y1y2y3y4_v1, P_y1y2y3y4_v2, ...
             overlap_with_init_state, overlap_fixed_basis_states, ~,~,~, ...
             comparison_determinants_mat, e_vals_mat_curr, ~, quasi_occ_t] = ...
                calc_noisy_tetron_ev_inst_basis_setup(corr_qp_init, e_vecs_init, ...
                e_vals_init, tetron_params, mu_vec, delta_t, t_init, t_final, ...
                time_evolution_method, basis_choice, majorana_zero_modes_ref, mu);
            overlap_inst_basis_states = NaN * ones(size(overlap_with_init_state));
        end

        if trial_no == 1
            e_vals_mat = e_vals_mat_curr; % only save first trial
        end

        comparison_determinants_mat_top(trial_no, :) = comparison_determinants_mat(1,:);
        comparison_determinants_mat_bot(trial_no, :) = comparison_determinants_mat(2,:);

        Z_exp_mat(trial_no,:) = Z_exp;
        X_exp_mat(trial_no,:) = X_exp;
        Y_exp_mat(trial_no,:) = Y_exp;
        P_y1y2y3y4_v1_mat(trial_no,:) = P_y1y2y3y4_v1;
        P_y1y2y3y4_v2_mat(trial_no,:) = P_y1y2y3y4_v2;
        overlap_with_init_state_mat(trial_no,:) = overlap_with_init_state;
        overlap_fixed_basis_states_mat(trial_no,:) = overlap_fixed_basis_states;
        overlap_inst_basis_states_mat(trial_no,:) = overlap_inst_basis_states;
    end

    %--------------------------------------%
    % Collect results
    %--------------------------------------%
    comparison_determinants_cell = {comparison_determinants_mat_top, comparison_determinants_mat_bot};
    fprintf('\n');
    run_time = toc(tStart);

    if strcmp(plot_results, 'plot_all_results') || strcmp(save_output, 'save_output')
        error('plot_all_results and save_output not supported in this code version');
    end

    %--------------------------------------%
    % Fit X_exp data
    %--------------------------------------%
    X_exp_mean = mean(X_exp_mat,1);
    Z_exp_mean = mean(Z_exp_mat,1);
    P_y1y2y3y4_v1_mean = mean(P_y1y2y3y4_v1_mat,1);
    P_y1y2y3y4_v2_mean = mean(P_y1y2y3y4_v2_mat,1);
    overlap_with_init_state_mean = mean(overlap_with_init_state_mat,1);
    overlap_fixed_basis_states_mean = mean(overlap_fixed_basis_states_mat,1);
    overlap_inst_basis_states_mean = mean(overlap_inst_basis_states_mat,1);

    if norm(imag(X_exp_mean)) > 1e-10
        warning('non-negligible imaginary component in X_exp_mean');
    end
    X_exp_mean = real(X_exp_mean);

    cosine_damped = 'exp(-a*x)*cos(b*x)';
    start_point = [fluctuation_amplitude*1e-2, 2*min(abs(e_vals_init))];
    weights = ones(size(t_vec));

    try
        X_exp_fit = fit(t_vec(:), X_exp_mean(:), cosine_damped, 'Start', start_point, 'Weights', weights);
        best_fit = exp(-X_exp_fit.a .* t_vec) .* cos(X_exp_fit.b .* t_vec);
        omega_0_from_fit = X_exp_fit.b;
        T2_from_fit = 1 ./ X_exp_fit.a;
    catch
        % Note this fit is not used in the uploaded scripts for our paper. 
        % So it is not an issue if this catch is reached. I use this in
        % some other (not-yet-uploaded) code. 
        T2_from_fit = NaN;
        omega_0_from_fit = NaN;
        best_fit = zeros(size(t_vec));
    end

    if strcmp(plot_results, 'plot_X_exp_only') || strcmp(plot_results, 'plot_all_results')
        warning('plot_all_results and save_output not supported in this code version');
    end

    %--------------------------------------%
    % MZM parity decay
    %--------------------------------------%
    linear_symb = '-a*x';
    start_point = [-1/1000];
    weights = ones(size(t_vec));
    log_mzm_parity = log(P_y1y2y3y4_v1_mean);
    parity_leakage_fit = fit(t_vec(:), log_mzm_parity(:), linear_symb, 'Start', start_point, 'Weights', weights);
    mzm_parity_decay = 1 ./ parity_leakage_fit.a;

    %--------------------------------------%
    % Frequency estimate
    %--------------------------------------%
    a_temp = sort(e_vals_init, 'ascend', 'ComparisonMethod','abs');
    omega_0_guess = a_temp(1) * 2;

    %--------------------------------------%
    % Save output
    %--------------------------------------%
    if strcmp(save_output, 'save_output')
        warning('save_output option not available in this version');
    end

end


%% Helper Functions
%--------------------------------------%
% Check valid inputs
%--------------------------------------%
function check_valid_inputs(save_output, plot_results)
    if ~strcmp(save_output, 'save_output') && ~strcmp(save_output, 'do_not_save_output')
        error('Invalid Input for save_output');
    end

    if ~strcmp(plot_results, 'plot_all_results') && ~strcmp(plot_results, 'plot_X_exp_only') && ...
       ~strcmp(plot_results, 'do_not_plot_results')
        error('Invalid Input for plot_results');
    end
end

