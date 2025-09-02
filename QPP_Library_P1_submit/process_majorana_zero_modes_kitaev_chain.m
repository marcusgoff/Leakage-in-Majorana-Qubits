%-------------------------------------------%
% Process Majorana Zero Modes Kitaev Chain
%-------------------------------------------%
% Inputs:
%   zero_modes         - 2N x 2 matrix where columns are the eigenvectors of
%                        a BdG Hamiltonian with near-zero energy. 
%   zero_modes_energies- 1 x 2N vector with corresponding energies of zero_modes.
%
% Outputs:
%   majorana_zero_modes      - 2N x 2 matrix where the columns are the
%                              rotated eigenvectors of zero_modes in the 
%                              Majorana zero mode basis. 
%   deloc_dirac_zero_mode    - 2N x 2 matrix containing |d_0> and |d_0^dag>
%                              where |d_0> is the delocalised fermionic mode
%                              corresponding to the majorana_zero_modes. 
%   test_neg_overlap_failed_flag - Optional flag indicating potential
%                              negative overlap issue (see notes below).
%
% Notes:
%   This function finds linear combinations of the near-zero energy eigenvectors
%   to produce Majorana zero modes maximally localised at either end of the chain.
%   The flag test_neg_overlap_failed_flag is raised if |d_0> may have a positive 
%   energy and |d_0^dag> may have a negative energy. This often occurs near 
%   the fixed point of the Kitaev chain, and the outputs can still be used.

function [majorana_zero_modes, deloc_dirac_zero_mode, test_neg_overlap_failed_flag] = ...
    process_majorana_zero_modes_kitaev_chain(zero_modes, zero_modes_energies)

    %% Process and Check Inputs
    if size(zero_modes,2) ~= 2 || mod(size(zero_modes,1),2)
        error('zero_modes must be a 2N x 2 matrix');
    end
    if (size(zero_modes_energies,2) ~= 1 && size(zero_modes_energies,1) ~= 1) || ...
            length(zero_modes_energies) ~= 2
        error('zero_modes_energies must be a 1 x 2 matrix');
    end   

    % Sort eigenvectors and eigenvalues
    [e_vecs, e_vals] = sort_eigenvectors_and_eigenvalues(...
        zero_modes, zero_modes_energies, 'descend');
    e_vec_pos = zero_modes(:,1); % positive energy
    e_vec_neg = zero_modes(:,2); % negative energy
    
    %% Locally Defined Variables
    norm_zero_thresh = 1e-10;
    N = length(zero_modes) / 2; 
    X = [0 1; 1 0];  
    
    %% Define a delocalised Majorana mode from e_vec_neg
    maj_deloc_1 = e_vec_neg + particle_hole_trans(e_vec_neg); % delocalised Majorana mode
    
    norm_sq_maj_deloc = maj_deloc_1' * maj_deloc_1; 
    if norm_sq_maj_deloc < norm_zero_thresh
        maj_deloc_1 = 1i * (e_vec_neg - particle_hole_trans(e_vec_neg)); 
        norm_sq_maj_deloc = maj_deloc_1' * maj_deloc_1; 
        if norm_sq_maj_deloc < norm_zero_thresh
            error('Need to change algorithm for creating delocalised Majorana mode'); 
        end
    end
    maj_deloc_1 = maj_deloc_1 / sqrt(norm_sq_maj_deloc);
    
    %% Define maj_deloc_2 as orthogonal eigenstate in near-zero energy subspace
    proj_near_zero = e_vec_neg * e_vec_neg' + e_vec_pos * e_vec_pos';
    if norm(proj_near_zero^2 - proj_near_zero) / (4*N^2) > 1e-10
        error('Error in proj_near_zero calculation');
    end
    proj_maj_deloc_2 = proj_near_zero - maj_deloc_1 * maj_deloc_1';
    [V,D] = eig(proj_maj_deloc_2, 'vector');
    [~, ind] = min(abs(D - 1));
    if abs(D(ind) - 1) > 1e-5
        error('Cannot find maj_deloc_2 using projector method');
    end
    maj_deloc_2_guess = V(:,ind);
    
    % Ensure maj_deloc_2 is Hermitian
    maj_deloc_2 = maj_deloc_2_guess + particle_hole_trans(maj_deloc_2_guess);
    norm_sq_maj_2_deloc = maj_deloc_2' * maj_deloc_2; 
    if norm_sq_maj_2_deloc < norm_zero_thresh
        maj_deloc_2 = 1i * (maj_deloc_2_guess - particle_hole_trans(maj_deloc_2_guess)); 
        norm_sq_maj_2_deloc = maj_deloc_2' * maj_deloc_2; 
        if norm_sq_maj_2_deloc < norm_zero_thresh
            error('Need to change algorithm for creating delocalised Majorana mode'); 
        end
    end
    maj_deloc_2 = maj_deloc_2 / sqrt(norm_sq_maj_2_deloc);
    
    %% Find Majorana operators with weight on either side of chain
    position_op = diag([1:N, 1:N]);
    theta_checks = linspace(-pi, pi, 1e5);
    theta_min_val = 0;
    min_val = Inf;
    exp_vals = zeros(size(theta_checks)); % DEBUG
    
    for ii = 1:length(theta_checks)
        gamma_check = cos(theta_checks(ii)) * maj_deloc_1 + sin(theta_checks(ii)) * maj_deloc_2;
        exp_position_op = real(gamma_check' * position_op * gamma_check);
        exp_vals(ii) = exp_position_op; % DEBUG
        if exp_position_op < min_val
            min_val = exp_position_op;
            theta_min_val = theta_checks(ii);
        end
    end
    
    gamma_1_temp = cos(theta_min_val) * maj_deloc_1 + sin(theta_min_val) * maj_deloc_2;
    gamma_2_temp = cos(theta_min_val + pi/2) * maj_deloc_1 + sin(theta_min_val + pi/2) * maj_deloc_2; 
    
    dirac_zero_mode_temp = 1/sqrt(2) * [(gamma_1_temp + 1i*gamma_2_temp), ...
                                        gamma_1_temp - 1i*gamma_2_temp]; 
    
    test_neg_overlap = abs(e_vec_neg' * dirac_zero_mode_temp);

    if test_neg_overlap(1) >= test_neg_overlap(2)
        gamma_2 = gamma_2_temp;
    else
        gamma_2 = -gamma_2_temp;
    end    
    if max(test_neg_overlap) < 0.98
        test_neg_overlap_failed_flag = 1;
    else
        test_neg_overlap_failed_flag = 0;
    end
    
    gamma_1 = gamma_1_temp; 
  
    majorana_zero_modes = [gamma_1, gamma_2];       
    deloc_dirac_zero_mode = zeros(size(majorana_zero_modes));
    deloc_dirac_zero_mode(:,1) = 1/sqrt(2) * (gamma_1 - 1i * gamma_2);
    deloc_dirac_zero_mode(:,2) = 1/sqrt(2) * (gamma_1 + 1i * gamma_2);
      
end


%% ------------------- %%
%    Local Functions
%-----------------------%

%% Particle Hole Transformation
function vec_PHT = particle_hole_trans(vec)
    N = length(vec) / 2;
    X = [0 1; 1 0]; 
    vec_PHT = kron(X, eye(N)) * conj(vec);
end
