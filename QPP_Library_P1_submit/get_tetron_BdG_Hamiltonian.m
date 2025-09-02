%-------------------------------%
% Get Tetron BdG Hamiltonian
%-------------------------------%
% Author: Marcus
% Date: 28.03.24
%
% This function returns the BdG Hamiltonian for the tetron qubit, modelled
% as two uncoupled Kitaev chains.
%
% Inputs:
%   mu         - chemical potentials; to specify different values for top/bottom
%                chain use a 2x1 vector
%   w          - hopping amplitude
%   delta      - pairing amplitude
%   N          - number of sites in one Kitaev chain (total sites = 2N)
%   BC         - boundary conditions, valid inputs: 
%                "OBC" (open b.c.) or "PBC" (periodic b.c.)
%   basis_str  - 'chain_1_chain_2' or 'c_c_dag'
%
% Output:
%   H          - 4N x 4N BdG Hamiltonian of the tetron

function H = get_tetron_BdG_Hamiltonian(mu, w, delta, N, BC, basis_str)

    if length(mu) == 2
        mu_1 = mu(1);
        mu_2 = mu(2); 
    elseif length(mu) == 1
        mu_1 = mu;
        mu_2 = mu;
    else
        error('mu must be a 1x1 or 2x1 vector');
    end

    [A_1, B_1] = get_A_B_matrices(mu_1, w, delta, N, BC);
    [A_2, B_2] = get_A_B_matrices(mu_2, w, delta, N, BC);
    
    if strcmp(basis_str, 'c_c_dag')
        A = [A_1, zeros(N); zeros(N), A_2];
        B = [B_1, zeros(N); zeros(N), B_2]; 
        H = [A, -conj(B); B, -conj(A)];
        
    elseif strcmp(basis_str, 'chain_1_chain_2')
        H_KC_1 = [A_1, -conj(B_1); B_1, -conj(A_1)];
        H_KC_2 = [A_2, -conj(B_2); B_2, -conj(A_2)];
        H = [H_KC_1, zeros(2*N); zeros(2*N), H_KC_2];
        
    else 
        error('Invalid input for: basis_str');
    end
        
end


%% Get A and B matrices for a single Kitaev chain
function [A, B] = get_A_B_matrices(mu, w, delta, N, BC)

    % Construct tridiagonal hopping and pairing matrices
    A = -mu * diag(ones(N,1)) - w * diag(ones(N-1,1),1) - w * diag(ones(N-1,1),-1);
    B = delta * diag(ones(N-1,1),1) - delta * diag(ones(N-1,1),-1);
    
    % Apply periodic boundary conditions if required
    if BC == "PBC"
        A(1,end) = -w; 
        A(end,1) = -w;
        B(1,end) = -delta; 
        B(end,1) = delta;
    end
    
end
