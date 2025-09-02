%--------------------------------%
% Get Parity Quench Limit
%--------------------------------%
% This function calculates the parity of the Kitaev-tetron immediately
% following a quench, using the expression:
%
% <P> = <psi_init | (-i2)^2 gamma_1 gamma_2 gamma_3 gamma_4 | psi_init>
% where |psi_init> is the state immediately before the quench and
%       gamma_i is the i-th MZM immediately after the quench.

function [P, gap_init, gap_final] = get_parity_quench_limit(mu_mean_init, mu_mean_final, ...
    mu_offset, w, delta, N, BC)

    %------------------------------%
    % Steps:
    % - Get initial basis eigenvectors (U_1)
    % - Express |+> state as covariance matrix in initial basis (M_1)
    % - Get final basis eigenvectors (U_2)
    % - Rotate M_1 to basis of new eigenvectors (M_2)
    % - Restrict M_2 to gamma_1 to gamma_4 subspace
    % - Calculate parity using calc_pfaffian_4x4(M_2_restr)
    %-----------------------------%
    
    mu_init = mu_mean_init(:) + mu_offset * [+1, -1].';
    mu_final = mu_mean_final(:) + mu_offset * [+1, -1].';

    % Initial tetron Hamiltonian
    H_tetron_init = get_tetron_BdG_Hamiltonian(mu_init, w, delta, N, BC, 'chain_1_chain_2');
    [e_vecs_init, e_vals_init, majorana_zero_modes_ref_init] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
        H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');   
    gap_init = sort(abs(e_vals_init)); 
    gap_init = gap_init(5);

    [~, corr_mat_site_init_basis, corr_mat_qp_init_basis] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "+"); 

    % Final tetron Hamiltonian
    H_tetron_final = get_tetron_BdG_Hamiltonian(mu_final, w, delta, N, BC, 'chain_1_chain_2');
 
    [e_vecs_final, e_vals_final, majorana_zero_modes_ref_final, comparison_determinants] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim(H_tetron_final(1:2*N,1:2*N), ...
        H_tetron_final((2*N+1):end, (2*N+1):end), 'dirac_zero_modes', ...
        majorana_zero_modes_ref_init);
    gap_final = sort(abs(e_vals_final)); 
    gap_final = gap_final(5);

    % Rotate initial site-basis covariance matrix to final QP basis
    omega_total = get_corr_cov_conversion_matrix(N);
    corr_mat_qp_final_basis = e_vecs_final.' * corr_mat_site_init_basis * conj(e_vecs_final);
    cov_mat_qp_curr = -1i * omega_total * (2 * corr_mat_qp_final_basis - eye(4*N)) * omega_total'; 

    % Restrict to MZM 4x4 subspace and take the Pfaffian
    mzm_indices_cov_qp = [1, N+1, 1+2*N, 1+3*N];
    P = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp));
    P = remove_negligible_imag_parts(P);
    
end


%% Calculate 4x4 Pfaffian 
% INPUT: M - 4x4 covariance matrix in Majorana basis (quasiparticle or site basis)
function pf = calc_pfaffian_4x4(M)
    pf = M(1,2)*M(3,4) - M(1,3)*M(2,4) + M(2,3)*M(1,4); 
end


%% Get omega_total conversion matrix
function omega_total = get_corr_cov_conversion_matrix(N)
    omega_s = 1/sqrt(2) * [eye(N), eye(N); 1i*eye(N), -1i*eye(N)]; 
    lam = [fliplr(eye(N)), zeros(N); zeros(N), eye(N)]; 
    omega_total = kron(eye(2), omega_s * lam);
end


%% Remove negligible imaginary parts
function y = remove_negligible_imag_parts(x)
    if max(abs(imag(x))) < 1e-12
        y = real(x);
    else
        y = x;
    end  
end
