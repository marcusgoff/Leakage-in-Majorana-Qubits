%-------------------------------------------------%
% Tetron Time Evolution - QP Correlation with Noise
% using various time evolution methods.
%-------------------------------------------------%
% Author: Marcus Goffage
%
% Adapted from calculate_noisy_tetron_evolution_dirac_qp_corr_mat to
% include the input time_evolution_method.
%
% Notes on inputs:
%   e_vecs - Tetron eigenvectors, ordered by chain-1 chain-2 first,
%            and by descending eigenvalues second.
%            This is achieved by using
%            diagonalise_uncoupled_tetron_via_Kitaev_chains.
%            Make sure you use this ordering, or this
%            function will return an error.
%
%   tetron_params = {mu_offset, w, delta, N, BC}
%
%   mu_vec - [2 x num_time_steps] values of mu (noisy) at each time step.
%            First row = top Kitaev chain,
%            Second row = bottom Kitaev chain.
%            Length as t_init:delta_t:t_final.
%
%   time_evolution_method - Choose method for obtaining the unitary U_t
%       'time_ev_diagonalise' - Diagonalise H_t to construct U_t
%       'time_ev_exp'         - Use MATLAB function expm to obtain U_t
%       'time_ev_1st_order'   - 1st order Taylor expansion for U_t
%       'time_ev_Magnus'      - Magnus expansion (Not Implemented)
%
%   basis_choice - Define Pauli operators based on:
%       'init_hamil_basis'    - Hamiltonian at t = 0 (Not Implemented)
%       'instant_hamil_basis' - Instantaneous Hamiltonian H(t)
%
%   delta_t - Either a scalar or vector of delta_t values.
%
% Outputs:
%   P_y1y2y3y4_v1 - Using localised Majoranas
%   P_y1y2y3y4_v2 - Using delocalised Majoranas (not implemented)
%   e_vals_mat    - (i,j) element is j-th instantaneous eigenvalue
%                   at time step t_i. Calculated only if
%                   basis_choice = 'instant_hamil_basis'.
%                   Otherwise NaN matrix.
%   quasi_occ_final - N x 2 matrix of quasiparticle occupations at t_final.
%                     Col 1 = KC1, Col 2 = KC2.
%                     Rows ordered by decreasing eigenvalue
%                     [d_n, ..., d_1].
%   quasi_occ_t     - 2N x length(t_vec) matrix of quasiparticles vs time.
%                     quasi_occ_t(:,end) = quasi_occ_final(:).
%                     Large matrix—don’t output unless needed.
%
%   Note: "fixed" basis states are computational basis states defined using
%   eigenvectors of H at t=0, as opposed to instantaneous basis states.
%
%   Overlaps outputted are abs(overlap).^2.

function [Z_exp, X_exp, Y_exp, P_y1y2y3y4_v1, P_y1y2y3y4_v2, overlap_with_init_state, ...
          overlap_fixed_basis_states, t_vec, first_site_occ, mu_vec, ...
          comparison_determinants_mat, e_vals_mat, quasi_occ_final, ...
          quasi_occ_t, overlap_inst_basis_states] = ...
          calc_noisy_tetron_ev_inst_basis_setup( ...
              corr_mat_init, e_vecs_init, e_vals_init, tetron_params, mu_vec, ...
              delta_t, t_init, t_final, time_evolution_method, basis_choice, ...
              majorana_zero_modes_ref, mu)

    %-------------------------------------------------%
    % Stub outputs for to-be-implemented variables
    %-------------------------------------------------%
    P_y1y2y3y4_v2 = NaN;


    %-------------------------------------------------%
    % Process Inputs
    %-------------------------------------------------%    
    if length(tetron_params) ~= 4
        error('tetron_params input must be a length 4 vector');
    end
    w     = tetron_params{1};
    delta = tetron_params{2};
    N     = tetron_params{3};
    BC    = tetron_params{4};

    N = length(e_vals_init) ./ 4;

    %-------------------------------------------------%
    % Time evolution setup
    %-------------------------------------------------%
    if isscalar(delta_t)
        t_vec = t_init:delta_t:t_final;
    else
        t_vec = t_init + [t_init, cumsum(delta_t(:)).'];
        if t_vec(end) ~= t_final
            warning('delta_t vector does not match t_init and t_final');
            t_final = t_vec(end);
        end
    end

    Z_exp                     = zeros(size(t_vec));
    X_exp                     = zeros(size(t_vec));
    Y_exp                     = zeros(size(t_vec));
    P_y1y2y3y4_v1             = zeros(size(t_vec));
    overlap_with_init_state   = zeros(size(t_vec));
    overlap_fixed_basis_states = zeros(size(t_vec));
    first_site_occ            = zeros(size(t_vec));
    quasi_occ_final           = zeros(N, 2);
    quasi_occ_t               = zeros(2*N, length(t_vec));

    omega_total = get_corr_cov_conversion_matrix(N);
    [~, ~, state_0_fixed_corr_qp] = get_tetron_init_cov_mat_general(e_vecs_init, "0");
    [~, ~, state_1_fixed_corr_qp] = get_tetron_init_cov_mat_general(e_vecs_init, "1");

    % Computation of state_0_fixed_corr_qp does not depend on fixed basis
    % eigenvectors as it uses a quasiparticle basis.
    state_0_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_0_fixed_corr_qp, omega_total);
    state_1_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_1_fixed_corr_qp, omega_total);

    corr_mat_curr = corr_mat_init; % In QP basis
    cov_mat_qp_init = -1i .* omega_total * (2*corr_mat_curr - eye(4*N)) * omega_total';
    cov_mat_qp_curr = cov_mat_qp_init;

    Z_exp(1) = cov_mat_qp_curr(1, 1+N);
    X_exp(1) = cov_mat_qp_curr(1, 1+2*N);
    Y_exp(1) = -cov_mat_qp_curr(1, 1+3*N);

    % Leakage
    mzm_indices_cov_qp = [1, N+1, 1+2*N, 1+3*N];
    P_y1y2y3y4_v1(1) = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp));

    % Overlaps
    overlap_with_init_state(1) = calc_overlap_cov_mats(cov_mat_qp_init, cov_mat_qp_curr, 2*N);
    overlap_fixed_basis_states(1) = calc_overlap_cov_mats(cov_mat_qp_curr, state_0_fixed_cov_qp, 2*N) ...
                                  + calc_overlap_cov_mats(cov_mat_qp_curr, state_1_fixed_cov_qp, 2*N);

    % Instantaneous and initial bases same at first time step
    overlap_inst_basis_states = overlap_fixed_basis_states;

    %-------------------------------------------------%
    % Delocalised MZM code
    %-------------------------------------------------%
    mu_init   = mu_vec(:,1);
    H_tetron  = get_tetron_BdG_Hamiltonian(mu, w, delta, N, BC, 'chain_1_chain_2');
    [e_vecs_deloc_init, e_vals_deloc_init, ~, comparison_determinants_init] = ...
        diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim( ...
            H_tetron(1:2*N,1:2*N), ...
            H_tetron((2*N+1):end, (2*N+1):end), ...
            'dirac_zero_modes', majorana_zero_modes_ref);

    comparison_determinants_mat      = zeros(2, length(t_vec));
    comparison_determinants_mat(:,1) = comparison_determinants_init(:);

    % Store eigenvalues
    switch basis_choice
        case 'instant_hamil_basis'
            e_vals_mat = zeros(length(t_vec), 4*N);
            e_vals_mat(1,:) = e_vals_deloc_init(:).';
        case 'init_hamil_basis'
            error('init_hamil_basis not implemented');
    end

    % Initial quasiparticle occupation (instantaneous = initial basis)
    quasi_occ_t(1:N, 1)     = diag(corr_mat_curr(1:N, 1:N));
    quasi_occ_t((1:N)+N, 1) = diag(corr_mat_curr((1:N)+2*N, (1:N)+2*N));
    %-------------------------------------------------%

    if nargout >= 5
        corr_mat_site_curr = conj(e_vecs_init) * corr_mat_curr * e_vecs_init.';
        first_site_occ(1) = corr_mat_site_curr(1,1);
    end

    %-------------------------------------------------%
    % Main loop over time steps
    %-------------------------------------------------%
    for ii = 2:length(t_vec)
        % Get current time evolution unitary
        mu_curr  = mu_vec(:,ii);
        H_tetron = get_tetron_BdG_Hamiltonian(mu_curr, w, delta, N, BC, 'chain_1_chain_2');
        H_evolve = e_vecs_init' * H_tetron * e_vecs_init; % In initial QP basis

        if isscalar(delta_t)
            delta_t_curr = delta_t;
        else
            delta_t_curr = delta_t(ii-1);
        end

        switch time_evolution_method
            case 'time_ev_exp'
                U_delta_t = expm(1i * H_evolve .* delta_t_curr);

            case 'time_ev_diagonalise' % Preferred method
                H_e = (H_evolve + H_evolve') / 2;
                [V_QP, D_evolve] = eig(H_e);
                D_evolve = real(D_evolve);
                U_delta_t = V_QP * diag(exp(1i * diag(D_evolve) .* delta_t_curr)) * V_QP';

            case 'time_ev_1st_order'
                H_e = (H_evolve + H_evolve') / 2;
                U_delta_t = eye(4*N) + 1i * H_e * delta_t_curr;

            otherwise
                error('Invalid value for time_evolution_method');
        end

        % Time evolution
        corr_mat_next = U_delta_t * corr_mat_curr * U_delta_t';
        corr_mat_curr = corr_mat_next; % In initial QP basis

        %-------------------------------------------------%
        % Calculate quantities
        %-------------------------------------------------%
        if strcmp(basis_choice, 'init_hamil_basis')
            error('init_hamil_basis not implemented');

        elseif strcmp(basis_choice, 'instant_hamil_basis')
            % Get covariance matrix in instantaneous Hamiltonian eigenbasis
            % (time evolution still done in initial basis)

            % Convert correlation matrix in initial QP basis to site basis
            corr_mat_curr_site = conj(e_vecs_init) * corr_mat_curr * e_vecs_init.';
            if ~strcmp(time_evolution_method, 'time_ev_diagonalise')
                [V_QP, D_evolve] = eig(H_e); % In initial Hamiltonian eigenbasis
            end

            % Convert from site basis to instantaneous QP basis
            [e_vecs_deloc_curr, e_vals_deloc_curr, ~, comparison_determinants_curr] = ...
                diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim( ...
                    H_tetron(1:2*N,1:2*N), ...
                    H_tetron((2*N+1):end, (2*N+1):end), ...
                    'dirac_zero_modes', majorana_zero_modes_ref);

            comparison_determinants_mat(:,ii) = comparison_determinants_curr(:);

            corr_mat_curr_deloc = e_vecs_deloc_curr.' * corr_mat_curr_site * conj(e_vecs_deloc_curr);
            cov_mat_qp_curr = -1i .* omega_total * (2*corr_mat_curr_deloc - eye(4*N)) * omega_total';

            % Extract instantaneous QP basis populations
            quasi_occ_t(1:N, ii)     = diag(corr_mat_curr_deloc(1:N, 1:N));
            quasi_occ_t((1:N)+N, ii) = diag(corr_mat_curr_deloc((1:N)+2*N, (1:N)+2*N));

        else
            error('Invalid Input for basis_choice');
        end

        %-------------------------------------------------%
        % Observables
        %-------------------------------------------------%
        Z_exp(ii) = cov_mat_qp_curr(1, 1+N);
        X_exp(ii) = cov_mat_qp_curr(1, 1+2*N);
        Y_exp(ii) = -cov_mat_qp_curr(1, 1+3*N);
        P_y1y2y3y4_v1(ii) = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp));

        % Overlaps
        overlap_with_init_state(ii) = calc_overlap_cov_mats(cov_mat_qp_init, cov_mat_qp_curr, 2*N);

        cov_mat_qp_init_basis_curr = -1i .* omega_total * (2*corr_mat_curr - eye(4*N)) * omega_total';
        overlap_fixed_basis_states(ii) = calc_overlap_cov_mats(cov_mat_qp_init_basis_curr, state_0_fixed_cov_qp, 2*N) ...
                                       + calc_overlap_cov_mats(cov_mat_qp_init_basis_curr, state_1_fixed_cov_qp, 2*N);

        % Instantaneous computational basis state overlap
        state_0_inst_cov_qp = state_0_fixed_cov_qp;
        state_1_inst_cov_qp = state_1_fixed_cov_qp;
        overlap_inst_basis_states(ii) = calc_overlap_cov_mats(cov_mat_qp_curr, state_0_inst_cov_qp, 2*N) ...
                                      + calc_overlap_cov_mats(cov_mat_qp_curr, state_1_inst_cov_qp, 2*N);

        if strcmp(basis_choice, 'instant_hamil_basis')
            e_vals_mat(ii,:) = e_vals_deloc_curr(:).';
        end

        if nargout >= 5
            corr_mat_site_curr = conj(e_vecs_init) * corr_mat_curr * e_vecs_init.';
            first_site_occ(ii) = corr_mat_site_curr(1,1);
        end
    end

    %-------------------------------------------------%
    % Post-processing
    %-------------------------------------------------%
    Z_exp                   = remove_negligible_imag_parts(Z_exp);
    X_exp                   = remove_negligible_imag_parts(X_exp);
    Y_exp                   = remove_negligible_imag_parts(Y_exp);
    P_y1y2y3y4_v1           = remove_negligible_imag_parts(P_y1y2y3y4_v1);
    overlap_with_init_state = remove_negligible_imag_parts(overlap_with_init_state);
    overlap_fixed_basis_states = remove_negligible_imag_parts(overlap_fixed_basis_states);
    overlap_inst_basis_states  = remove_negligible_imag_parts(overlap_inst_basis_states);

    % Extract final QP populations
    quasi_occ_t          = remove_negligible_imag_parts(quasi_occ_t);
    quasi_occ_final(:,1) = quasi_occ_t(1:N, end);
    quasi_occ_final(:,2) = quasi_occ_t((1+N):end, end);
end

%-------------------------------------------------%
% Helper functions
%-------------------------------------------------%
function omega_total = get_corr_cov_conversion_matrix(N)
    omega_s = 1/sqrt(2) * [eye(N), eye(N); 1i*eye(N), -1i*eye(N)];
    lam = [fliplr(eye(N)), zeros(N); zeros(N), eye(N)];
    omega_total = kron(eye(2), omega_s*lam);
end

function cov_mat_qp = convert_corr_qp_to_cov_qp(corr_mat_qp, omega_total)
    N = size(corr_mat_qp,1) ./ 4;
    % QP covariance matrix
    cov_mat_qp = -1i .* omega_total * (2*corr_mat_qp - eye(4*N)) * omega_total';
end

function y = remove_negligible_imag_parts(x)
    if max(abs(imag(x))) < 1e-12
        y = real(x);
    else
        y = x;
    end
end

% Calculate 4x4 Pfaffian
% Input: M must be a 4x4 covariance matrix in Majorana basis
% (quasiparticle or site basis).
function pf = calc_pfaffian_4x4(M)
    pf = M(1,2)*M(3,4) - M(1,3)*M(2,4) + M(2,3)*M(1,4);
end

% Calculate Gaussian state overlap from covariance matrices
function overlap = calc_overlap_cov_mats(A, B, N)
    overlap = 2^(-N) * sqrt(det(A + B));
end
