% This function computes the block diagonal affinity matrix based on
% estimated parameters.
%
% For details, see:
%
% [1] A. Taştan, M. Muma and A. M. Zoubir, "Fast and Robust Sparsity-Aware
% Block Diagonal Representation," in IEEE Transactions on Signal Processing,
% 2023.
%
% Copyright (C) 2023 Aylin Tastan. All rights reserved.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% INPUTS:
%        W             :(numeric) zero diagonal symmetric affinity matrix
%                       of size (N-N_{I}) x (N-N_{I}) if Type I outlier
%                       removal has been preferred in Step 1.1.
%                       (N x N otherwise)
%        W_sim_hat     :(numeric) similarity coefficients matrix of size
%                       K_hat x K_hat. It contains similarity coefficients
%                       that the blocks are concentrated around.
%        bsize_vec_hat :(numeric) estimated block size vector of size
%                       K_hat x 1
%        K_hat         :(numeric) estimated number of blocks
%
% OUTPUTS:
%         W_block_hat  :(numeric) zero diagonal block affinity matrix of
%                       size (N-N_{o_I}) x (N-N_{o_I}) if Type I outlier
%                       removal has been preferred in Step 1.1.
%                       (N x N otherwise)
%
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W_block_hat = compute_block_diagonal_affinity_matrix(W,W_sim_hat,bsize_vec_hat,K_hat)

if (all(max(W_sim_hat(find(triu(ones(K_hat),1)))) < diag(W_sim_hat)))

    W_block_hat = W;
    W_block_hat(W_block_hat < min(diag(W_sim_hat))) = 0;

else

    block_endpoints_hat = cumsum(bsize_vec_hat);
    blocks_idx = find(diag(W_sim_hat)< max(W_sim_hat(find(triu(ones(K_hat),1)))));
    W_block_hat = W;
    W_block_hat(W_block_hat < max(W_sim_hat(find(triu(ones(K_hat),1))))) = 0;

    for i = 1:length(blocks_idx)

        block_i = blocks_idx(i);

        if(block_i > 1)

            W_block_i = W(block_endpoints_hat(block_i-1)+1:block_endpoints_hat(block_i),:);
            W_block_i(W_block_i < min(diag(W_sim_hat))) = 0;
            W_block_hat(block_endpoints_hat(block_i-1)+1:block_endpoints_hat(block_i),:) = W_block_i;
            W_block_hat(:,block_endpoints_hat(block_i-1)+1:block_endpoints_hat(block_i)) = W_block_i.';

        else

            W_block_i = W(1:block_endpoints_hat(block_i),:);
            W_block_i(W_block_i < min(diag(W_sim_hat))) = 0;
            W_block_hat(1:block_endpoints_hat(block_i),:) = W_block_i;
            W_block_hat(:,1:block_endpoints_hat(block_i)) = W_block_i.';

        end

        clear block_i W_block_i
    end

end%endif
end%endfunction