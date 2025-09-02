%-------------------------------------------%
% Get Tetron Initialised Covariance Matrix
%-------------------------------------------%
% Author: Marcus Goffage
%
% Computes the covariance matrices for the computational basis states of the
% tetron:
%      |0> = |even>|even>
%      |1> = |odd>|odd>
%
% INPUTS:
%   e_vecs - eigenvectors of the tetron qubit in the chain_1_chain_2
%            basis. For each chain the eigenvalues must be ordered in descending
%            order. 
%   state  - allowed inputs: "+", "-", "0", "1"
%      
% OUTPUTS:
%   cov_site  - covariance matrix in the site Majorana basis
%   corr_site - correlation matrix in the site Dirac fermion basis
%   corr_qp   - correlation matrix in the Dirac fermion quasiparticle basis
%
% Comments:
%   corr_qp is in the quasiparticle basis, and so it does not depend on
%   the e_vecs.

function [cov_site, corr_site, corr_qp] = get_tetron_init_cov_mat_general(e_vecs, state)
    
    check_inputs(state); % only checks state at present   
   
    N = size(e_vecs,1) / 4;
    
    cov_0 = zeros(4*N);
    cov_1 = zeros(4*N);
     
    % Create correlation matrices in QP basis (Dirac fermions)
    A = zeros(N); A(N,N) = 1; 
    B = eye(N); B(1,1) = 0; 
   
    corr_qp_KC_even = [zeros(N), zeros(N); zeros(N), eye(N)];
    corr_qp_KC_odd  = [A, zeros(N); zeros(N), B];
   
    corr_qp_0 = [corr_qp_KC_even, zeros(2*N); zeros(2*N), corr_qp_KC_even];
    corr_qp_1 = [corr_qp_KC_odd,  zeros(2*N); zeros(2*N), corr_qp_KC_odd];   
   
    % Create off-diagonal block matrices for correlations in +/- state
    C = zeros(2*N);
    C(N+1, N) = 1i; 
    C(N, N+1) = 1i;  
    gamma_plus_off_diag = [zeros(2*N), C; -C, zeros(2*N)];
    % gamma_plus_off_diag is the off-diagonal terms (x2) for the "+" state
    % The "-" state simply has the negative of the off-diagonal terms.
      
    switch state
        case "0"
            corr_qp = corr_qp_0;            
        case "1"
            corr_qp = corr_qp_1;
        case "+"
            corr_qp = 0.5*(corr_qp_0 + corr_qp_1 + gamma_plus_off_diag);
        case "-"
            corr_qp = 0.5*(corr_qp_0 + corr_qp_1 - gamma_plus_off_diag);
        case "+i"
            error('case not yet coded');        
        case "-i"
            error('case not yet coded');
    end      
   
    % Convert correlation matrix to site basis (Dirac fermions)
    corr_site = conj(e_vecs) * corr_qp * e_vecs.';
    
    % Convert to covariance matrices in site basis (Majorana fermions)
    omega_s = 1/sqrt(2) * [eye(N), eye(N); 1i*eye(N), -1i*eye(N)]; 
    cov_site = -1i * kron(eye(2), omega_s) * (2*corr_site - eye(4*N)) * kron(eye(2), omega_s'); 
   
end


%% Input check function
function check_inputs(state)
    valid_states = {'0', '1', '+', '-', '+i', '-i'};
    if ~ismember(state, valid_states)
        error('Invalid input for parameter ''state''');
    end
end
