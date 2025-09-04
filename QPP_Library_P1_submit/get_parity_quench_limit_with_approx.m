% -------------------------------------------- %
% L_p (L_odd) Sudden Limit Approximation
% -------------------------------------------- %


function [P, P_approx, gap_init, gap_final] = get_parity_quench_limit_with_approx(mu_mean_init, mu_mean_final, ...
    mu_offset, w, delta, N, BC)

    %------------------------------%
    % Get initial basis eigenvectors. Store in U_1.

    % Express |+> state as a covariance matrix in the basis of the initial
    % basis eigenvectors. Store in M_1.

    % Get final basis eigenvectors. store in U_2.

    % Rotate M_1 to basis of new eigenvectors using U_2. Store in M_2.

    % Restrict M_2 to gamma_1 to gamma_4 subspace in the new basis. 

    % Use calc_pfaffian_4x4(M_2_restr) to get the parity of the initial state in
    % the final basis. Store result in P. 
    %-----------------------------%
    mu_init = mu_mean_init(:) + mu_offset*[+1, -1].';
    mu_final = mu_mean_final(:) + mu_offset*[+1, -1].';


    H_tetron_init = get_tetron_BdG_Hamiltonian(mu_init, w, delta, N, BC, 'chain_1_chain_2');
    [e_vecs_init, e_vals_init, mzm_mat_init] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
        H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');   
    gap_init = sort(abs(e_vals_init)); 
    gap_init = gap_init(5);


    [~, corr_mat_site_init_basis, corr_mat_qp_init_basis] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "+"); 

    H_tetron_final = get_tetron_BdG_Hamiltonian(mu_final, w, delta, N, BC, 'chain_1_chain_2');

    [e_vecs_final, e_vals_final, mzm_mat_final] = ...
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_final(1:2*N,1:2*N), ...
            H_tetron_final((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');    

    gap_final = sort(abs(e_vals_final)); 
    gap_final = gap_final(5);

    % Now rotate corr_mat_site_init_basis to cov_mat_qp_final_basis. Same
    % state different basis.

    omega_total = get_corr_cov_conversion_matrix(N);
    corr_mat_qp_final_basis = e_vecs_final.'*corr_mat_site_init_basis*conj(e_vecs_final);
    cov_mat_qp_curr = -1i.*omega_total*(2*corr_mat_qp_final_basis - eye(4*N))*omega_total'; 


    % Restrict to the MZM 4x4 matrix and take the Pfaffian
    mzm_indices_cov_qp = [1,N+1, 1+2*N, 1+3*N];
    P = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp));
    P = remove_negligible_imag_parts(P);


    % Get Approximate Parity Leakage 
    P_approx_a = diag(mzm_mat_final'*mzm_mat_init);
    P_approx_a = prod(remove_negligible_imag_parts(P_approx_a));

    P_approx = P_approx_a;
end



%% Calculate 4x4 Pfaffian 
% INPUT: M - must be a 4x4 covariance matrix in majorana basis
% (quasiparticle or site basis is good). 

function pf = calc_pfaffian_4x4(M)
        pf = M(1,2)*M(3,4) - M(1,3)*M(2,4) + M(2,3)*M(1,4); 
end

%% Get omega_total 

function omega_total = get_corr_cov_conversion_matrix(N)
    omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
    lam = [fliplr(eye(N)), zeros(N); zeros(N), eye(N)]; 
    omega_total = kron(eye(2),omega_s*lam);
end

%% Remove negligible imag parts

function y = remove_negligible_imag_parts(x)

    if max(abs(imag(x))) < 1e-12 % & ...
           % norm(imag(x))/norm(real(x)) < 1e-8
        y = real(x);
    else
        y = x;
    end  
    
end