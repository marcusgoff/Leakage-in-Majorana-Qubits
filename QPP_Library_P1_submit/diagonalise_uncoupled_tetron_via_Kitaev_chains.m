%-------------------------------------------------%
% Diagonalise Uncoupled Tetron via Kitaev Chains
%-------------------------------------------------%
%
% Given two Kitaev chain BdG Hamiltonians, return the eigenvectors and
% eigenvalues of the tetron corresponding to these two uncoupled Kitaev
% chains.
% This code also runs process_majorana_zero_modes_kitaev_chain to ensure
% the eigenvectors of the near-zero degenerate manifold (if it exists)
% correspond to the localised Majoranas on the ends of the chain.
%
% This code works in the "chain_1_chain_2" basis only. 
% A future addition would be to add a "basis_str" input and write code for 
% both basis options.
%
% INPUTS:
%   H_1              - 2N x 2N BdG Hamiltonian of the 1st Kitaev Chain
%                      i.e. H_tetron(1:N,1:N)
%   H_2              - 2N x 2N BdG Hamiltonian of the 2nd Kitaev Chain
%                      i.e. H_tetron((N+1):end, (N+1):end)
%   zero_mode_type   - either 'majorana_zero_modes' or 'dirac_zero_modes'
%
% OUTPUTS:
%   e_vecs            - 2N eigenvectors of full tetron in chain_1_chain_2 basis
%   e_vals            - 2N eigenvalues of full tetron in chain_1_chain_2 basis
%   majorana_zero_modes - of full tetron in chain_1_chain_2 basis
%

function [e_vecs, e_vals, majorana_zero_modes] = ...
    diagonalise_uncoupled_tetron_via_Kitaev_chains(H_1, H_2, zero_mode_type)

    % --- Validate inputs --- %
    if mod(size(H_1,2),2) || mod(size(H_1,1),2) || size(H_1,1) ~= size(H_1,2)
        error('H_1 must be a 2N x 2N matrix');
    end

    if mod(size(H_2,2),2) || mod(size(H_2,1),2) || size(H_2,1) ~= size(H_2,2)
        error('H_2 must be a 2N x 2N matrix');
    end

    if ~strcmp(zero_mode_type, 'majorana_zero_modes') && ...
       ~strcmp(zero_mode_type, 'dirac_zero_modes')
        error('Invalid input for zero_mode_type');
    end

    if strcmp(zero_mode_type, 'majorana_zero_modes')
        error('I have not yet coded up majorana_zero_modes input option');
    end

    % --- Initialize --- %
    N = length(H_1)/2;
    e_vecs = zeros(4*N);
    e_vals = zeros(4*N,1);
    majorana_zero_modes = zeros(4*N, 4);

    % --- Diagonalize each Kitaev chain --- %
    [e_vecs_1, e_vals_1] = eig(H_1, 'vector');
    [e_vecs_2, e_vals_2] = eig(H_2, 'vector');

    % --- Process Majorana zero modes --- %
    [e_vecs_1, e_vals_1, mzm_1] = process_majorana_zero_modes_kitaev_chain_wrapper(...
        e_vecs_1, e_vals_1, 'dirac_zero_modes');
    [e_vecs_2, e_vals_2, mzm_2] = process_majorana_zero_modes_kitaev_chain_wrapper(...
        e_vecs_2, e_vals_2, 'dirac_zero_modes');

    % --- Construct eigenvectors of the two uncoupled Kitaev chains --- %
    e_vecs(1:(2*N), 1:(2*N)) = e_vecs_1;
    e_vecs((2*N+1):end, (2*N+1):end) = e_vecs_2;

    % --- Combine eigenvalues --- %
    e_vals = [e_vals_1(:); e_vals_2(:)];

    % --- Assemble Majorana zero modes --- %
    majorana_zero_modes(1:(2*N), 1:2) = mzm_1;
    majorana_zero_modes((2*N+1):end, 3:4) = mzm_2;

end
