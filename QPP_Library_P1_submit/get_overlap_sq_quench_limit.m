%--------------------------------%
% Get Overlap Squared Quench Limit
%--------------------------------%
% This function calculates the overlap-mod-squared of the Kitaev-tetron immediately
% following a quench, using the following expression:
%
% overlap_sq = |<0_f|psi_init>|.^2 + |<1_f|psi_init>|.^2.
%
% Covariance matrices are used to represent all the states and to calculate overlaps.
% The quantity L_c0 is then computed as L_c0 = 1 - overlap_sq.
%
% All calculations are performed in the instantaneous basis.

function [overlap_sq, gap_init, gap_final] = get_overlap_sq_quench_limit(mu_mean_init, mu_mean_final, ...
    mu_offset, w, delta, N, BC)
 
    %------------------------------%
    % From parity version of this code:
    % - Get initial basis eigenvectors (U_1)
    % - Express |+> state as a covariance matrix in initial basis (M_1)
    % - Get final basis eigenvectors (U_2)
    % - Covariance matrices for |0_f> and |1_f> in final basis (generic forms)
    % - Rotate M_1 to the new basis using U_2 (M_2)
    % - Calculate |<0_f|psi_init>|.^2 and |<1_f|psi_init>|.^2
    % - Return overlap_sq = |<0_f|psi_init>|.^2 + |<1_f|psi_init>|.^2
    %-----------------------------%
    
    overlap_sq = 0; % STUB

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

    % Covariance matrices for |0_f> and |1_f>
    omega_total = get_corr_cov_conversion_matrix(N);
    [~, ~, state_0_final_corr_qp] = get_tetron_init_cov_mat_general(e_vecs_init, "0"); 
    [~, ~, state_1_final_corr_qp] = get_tetron_init_cov_mat_general(e_vecs_init, "1"); 
    state_0_final_cov_qp = convert_corr_qp_to_cov_qp(state_0_final_corr_qp, omega_total);
    state_1_final_cov_qp = convert_corr_qp_to_cov_qp(state_1_final_corr_qp, omega_total);

    % Rotate initial site-basis covariance matrix to final QP basis
    corr_mat_qp_final_basis = e_vecs_final.' * corr_mat_site_init_basis * conj(e_vecs_final);
    cov_mat_qp_curr = -1i * omega_total * (2 * corr_mat_qp_final_basis - eye(4*N)) * omega_total'; 

    % Calculate the modulus-squared overlap
    overlap_sq = calc_overlap_cov_mats(state_0_final_cov_qp, cov_mat_qp_curr, 2*N) + ...
                 calc_overlap_cov_mats(state_1_final_cov_qp, cov_mat_qp_curr, 2*N);
   
    overlap_sq = remove_negligible_imag_parts(overlap_sq);

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


%% Calculate Gaussian state overlap from covariance matrices
function overlap = calc_overlap_cov_mats(A, B, N) 
    overlap = 2^(-N) * sqrt(det(A + B)); 
end


%% Convert correlation QP matrix to covariance QP matrix
function cov_mat_qp = convert_corr_qp_to_cov_qp(corr_mat_qp, omega_total)
    N = size(corr_mat_qp, 1) / 4;
    cov_mat_qp = -1i * omega_total * (2 * corr_mat_qp - eye(4*N)) * omega_total'; 
end
