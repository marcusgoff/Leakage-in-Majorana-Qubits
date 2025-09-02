%-------------------------------------------%
% Process Majorana Zero Modes Kitaev Chain Wrapper - No Optimisation
%-------------------------------------------%
% Author: Marcus Goffage
%
% This function runs a simpler and faster version of
% process_majorana_zero_modes_kitaev_chain.m where the optimisation for the
% maximally localised Majoranas has been removed. 
% This means the Majoranas returned here are not guaranteed to be maximally
% localised; however, the MZMs will still be MZMs (not a mix of Majorana and Dirac
% fermions). 
%
% A check is included to ensure the parity of < 2i gamma_1 gamma_2>
% is consistent (there are no reflections, only rotations in the MZM subspace). 
%
% Intended for use at each time step in time evolution for the
% instantaneous Pauli matrices basis.
%
% Inputs:
%         e_vecs - eigenvectors of a BdG Hamiltonian
%         e_vals - eigenvalues of a BdG Hamiltonian
%         zero_mode_type - 'majorana_zero_modes' or 'dirac_zero_modes'
%         majorana_zero_modes_ref - 4N x 2 matrix of reference Majorana zero modes
%
% Outputs:
%         e_vecs_out - processed eigenvectors
%         e_vals_out - corresponding eigenvalues (descending order)
%         majorana_zero_modes - processed Majorana zero modes
%         comp_det - determinant of comparison matrix for parity check

function [e_vecs_out, e_vals_out, majorana_zero_modes, comp_det] = ...
    process_majorana_zero_modes_kitaev_chain_wrapper_no_optim(e_vecs, e_vals, zero_mode_type, majorana_zero_modes_ref)

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

    % Process MZMs without localisation optimisation
    [majorana_zero_modes, deloc_dirac_zero_mode, comp_det] = ...
        process_majorana_zero_modes_kitaev_chain_no_optim(zero_modes, ...
        zero_modes_energies, majorana_zero_modes_ref);

    % Replace eigenvectors depending on zero_mode_type
    if strcmp(zero_mode_type, 'majorana_zero_modes')
        e_vecs_out(:, [N, N+1]) = majorana_zero_modes;
    end
    if strcmp(zero_mode_type, 'dirac_zero_modes')
        e_vecs_out(:, [N, N+1]) = deloc_dirac_zero_mode;
    end

end

%%------------------------------------------------------------------------%
%% Process Majorana Zero Modes Kitaev Chain - No Optimisation
%%------------------------------------------------------------------------%
function [majorana_zero_modes, deloc_dirac_zero_mode, comp_det] = ...
    process_majorana_zero_modes_kitaev_chain_no_optim(zero_modes, zero_modes_energies, majorana_zero_modes_ref)

    %% Input checks
    if size(zero_modes,2) ~= 2 || mod(size(zero_modes,1),2)
        error('zero_modes must be a 2N x 2 matrix');
    end
    if size(zero_modes_energies,2) ~= 1 && size(zero_modes_energies,1) ~= 1 || ...
            length(zero_modes_energies) ~= 2
        error('zero_modes_energy must be a 1 x 2 matrix');
    end

    [e_vecs, e_vals] = sort_eigenvectors_and_eigenvalues(zero_modes, zero_modes_energies, 'descend');
    e_vec_pos = zero_modes(:,1); % positive energy
    e_vec_neg = zero_modes(:,2); % negative energy

    %% Define delocalised Majorana modes
    norm_zero_thresh = 1e-10;
    N = length(zero_modes)/2;

    maj_deloc_1 = e_vec_neg + particle_hole_trans(e_vec_neg);
    norm_sq_maj_deloc = maj_deloc_1'*maj_deloc_1; 
    if norm_sq_maj_deloc < norm_zero_thresh
       maj_deloc_1  = 1i*(e_vec_neg - particle_hole_trans(e_vec_neg)); 
       norm_sq_maj_deloc = maj_deloc_1'*maj_deloc_1; 
       if norm_sq_maj_deloc < norm_zero_thresh
          error('Need to change algorithm for creating delocalised Majorana mode'); 
       end
    end
    maj_deloc_1 = maj_deloc_1./sqrt(norm_sq_maj_deloc);

    % Orthogonal Majorana in near-zero subspace
    proj_near_zero = e_vec_neg*e_vec_neg' + e_vec_pos*e_vec_pos';
    if norm(proj_near_zero^2 - proj_near_zero)./(4*N^2) > 1e-10
        error('Error in proj_near_zero calculation');
    end
    proj_maj_deloc_2 = proj_near_zero - maj_deloc_1*maj_deloc_1';
    [V,D] = eig(proj_maj_deloc_2, 'vector');
    [~, ind] = min(abs(D - 1));
    if abs(D(ind)-1) > 1e-5
        error('Cannot find maj_deloc_2 using projector method');
    end
    maj_deloc_2_guess = V(:,ind);
    maj_deloc_2 = maj_deloc_2_guess + particle_hole_trans(maj_deloc_2_guess);
    norm_sq_maj_2_deloc = maj_deloc_2'*maj_deloc_2; 
    if norm_sq_maj_2_deloc < norm_zero_thresh
       maj_deloc_2  = 1i*(maj_deloc_2_guess - particle_hole_trans(maj_deloc_2_guess)); 
       norm_sq_maj_2_deloc = maj_deloc_2'*maj_deloc_2; 
       if norm_sq_maj_2_deloc < norm_zero_thresh
          error('Need to change algorithm for creating delocalised Majorana mode'); 
       end
    end
    maj_deloc_2 = maj_deloc_2./sqrt(norm_sq_maj_2_deloc);

    %% Ensure consistency with reference MZMs
    gamma_1_ref = majorana_zero_modes_ref(:,1);
    gamma_2_ref = majorana_zero_modes_ref(:,2);

    comp_matrix = [maj_deloc_1'*gamma_1_ref, maj_deloc_2'*gamma_1_ref; ...
                   maj_deloc_1'*gamma_2_ref, maj_deloc_2'*gamma_2_ref];
    comp_det = det(comp_matrix);

    gamma_1 = maj_deloc_1; 
    if sign(real(comp_det)) == 1
        gamma_2 = maj_deloc_2;
    else
        gamma_2 = -maj_deloc_2;
        % Reconstruct comparison matrix for safety check
        comp_matrix_new = [gamma_1'*gamma_1_ref, gamma_2'*gamma_1_ref; ...
                           gamma_1'*gamma_2_ref, gamma_2'*gamma_2_ref];
        if sign(real(det(comp_matrix_new))) == -1
            warning('Comparison matrix method incorrect.');
        end
    end

    if abs(imag(comp_det)./real(comp_det)) > 1e-5
        warning('Comparison matrix method incorrect');
    end

    majorana_zero_modes = [gamma_1, gamma_2];       
    deloc_dirac_zero_mode = zeros(size(majorana_zero_modes));
    deloc_dirac_zero_mode(:,1) = 1/sqrt(2)*(gamma_1 - 1i*gamma_2);
    deloc_dirac_zero_mode(:,2) = 1/sqrt(2)*(gamma_1 + 1i*gamma_2);

end

%%------------------------------------------------------------------------%
%% Local Function: Particle Hole Transformation
%%------------------------------------------------------------------------%
function vec_PHT = particle_hole_trans(vec)
    N = length(vec)/2;
    X = [0 1; 1 0]; 
    vec_PHT = kron(X, eye(N))*conj(vec);
end
