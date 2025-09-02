%-------------------------------------------%
% Process Majorana Zero Modes Kitaev Chain Wrapper
%-------------------------------------------%
% Author: Marcus Goffage
%
% This is a basic wrapper function for
% process_majorana_zero_modes_kitaev_chain. Consider the matrix V, obtained 
% from [V, D] = eig(H_BdG) - where H_BdG is a BdG Hamiltonian for a Kitaev 
% chain with at most two Majorana modes. 
% The purpose of this wrapper is to replace the two 
% eigenvectors in the near-degenerate zero energy manifold with either the
% two Majoranas that are localised on the left and right hand sides of the
% chain or the quasiparticle creation and annihilation operators which
% correspond to superpositions of these two modes. 
%
% Inputs:
%         e_vecs - eigenvectors of a BdG Hamiltonian
%         e_vals - eigenvalues of a BdG Hamiltonian
%         zero_mode_type - either 'majorana_zero_modes' or 'dirac_zero_modes'
%
% Outputs:
%         e_vecs_out - processed eigenvectors where any near-zero
%                      energy eigenvectors have been rotated to Majorana
%                      zero modes with maximised weight on each side of the
%                      chain, or to the corresponding Dirac fermion modes.
%         e_vals_out - eigenvalues adjusted to correspond to processed
%                      e_vecs_out (descending order)
%         majorana_zero_modes - the processed Majorana zero modes
%         dirac_zero_modes - the corresponding Dirac fermion modes
%
% Note:
%         If the MZM splitting is non-zero then the outputs are not technically
%         all eigenvectors when zero_mode_type -> 'majorana_zero_modes'.

function [e_vecs_out, e_vals_out, majorana_zero_modes, dirac_zero_modes] = ...
    process_majorana_zero_modes_kitaev_chain_wrapper(e_vecs, e_vals, zero_mode_type)

    % Input validation
    if mod(size(e_vecs,2),2) || mod(size(e_vecs,1),2) || ...
         size(e_vecs,1) ~= size(e_vecs,2)
        error('e_vecs must be a 2N x 2N matrix');
    end
    if mod(size(e_vals,2),2) && mod(size(e_vals,1),2)
        error('e_vals must be a length 2N vector');
    end
    if ~strcmp(zero_mode_type, 'majorana_zero_modes') && ...
            ~strcmp(zero_mode_type, 'dirac_zero_modes') 
        error('Invalid input for zero_mode_type');
    end

    N = length(e_vals)./2;

    % Sort eigenvectors and eigenvalues
    [e_vecs, e_vals] = sort_eigenvectors_and_eigenvalues(e_vecs, e_vals, 'descend');

    e_vecs_out = e_vecs;
    e_vals_out = e_vals;

    % Extract near-zero energy manifold
    zero_modes = e_vecs(:, [N, N+1]);
    zero_modes_energies = e_vals([N, N+1]);

    % Process to Majorana zero modes
    [majorana_zero_modes, deloc_dirac_zero_mode] = ...
        process_majorana_zero_modes_kitaev_chain(zero_modes, zero_modes_energies);

    % Replace eigenvectors depending on zero_mode_type
    if strcmp(zero_mode_type, 'majorana_zero_modes')
        e_vecs_out(:, [N, N+1]) = majorana_zero_modes;
    end
    if strcmp(zero_mode_type, 'dirac_zero_modes')
        e_vecs_out(:, [N, N+1]) = deloc_dirac_zero_mode;
    end   

    dirac_zero_modes = deloc_dirac_zero_mode;

end
