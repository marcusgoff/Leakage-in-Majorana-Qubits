%-----------------------------------%
% Sort Eigenvectors and Eigenvalues
%-----------------------------------%
% Author: Marcus Goffage
%
% Utility function to sort eigenvectors and eigenvalues according to the
% eigenvalues.
%
% Inputs:
%   V, D   - obtained from [V,D] = eig(A,'vector'), where A is a matrix
%   order  - sorting method: 'ascend', 'descend', or 'ascend_mag_pos_neg'.
%            'ascend_mag_pos_neg' can only be used for BdG vectors.

function [V, D] = sort_eigenvectors_and_eigenvalues(V, D, order)

    if strcmp(order, 'ascend') || strcmp(order, 'descend')
        [D, ind] = sort(D, order);
        V = V(:,ind);

    elseif strcmp(order, 'ascend_mag_pos_neg')
        N = length(D)/2;
        if mod(N,2) ~= 0
            error('Input order = ascend_mag_pos_neg is only for BdG eigenvectors');
        end

        [D, ind] = sort(D, 'descend');
        V = V(:,ind);

        D(1:N) = fliplr(D(1:N));
        ind_2 = N:-1:1;
        V(:,1:N) = V(:,ind_2);

        % This method works because BdG has particle-hole symmetry (PHS).
        % Done this way to ensure completely zero eigenvalues are split
        % as expected for this ordering.
        %
        % Note: This option is not used to generate any results in our paper.

    else
        error('Invalid input for ''order''');
    end

end
