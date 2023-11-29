% This function generates matrix of candidate block sizes (Step 2.1).
%
% For details, see:
%
% [1] A. Taştan, M. Muma and A. M. Zoubir, "Fast and Robust Sparsity-Aware
% Block Diagonal Representation," in IEEE Transactions on Signal Processing,
% 2023.
%
% Copyright (C) 2023 Aylin Tastan. All rights reserved.
%
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
%
% INPUTS:
%        tau                 :(numeric) changepoint locations vector of
%                             maximum size N_cmax x 1
%        N                   :(numeric) number of samples valued as N-N_{I}
%                             if Type I outlier removal has been preferred
%                             in Step 1.1.(N otherwise)
%        K                   :(numeric) number of pieces in the linear function
%        min_num_samples     :(numeric) minimum number of samples in the blocks
%
%
% OUTPUTS:
%        cand_block_size_mat :(numeric) matrix of candidate block sizes. It
%                             has maximum size of zeta x K (default).
%
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cand_block_size_mat = compute_candidate_block_sizes(tau,N,K,min_num_samples)

%Find the candidate changepoints in the piece-wise linear function
cand_limits = nchoosek(tau,K-1);
ind_cand_changepoints = [ones(size(cand_limits,1),1),cand_limits,(N+1)*ones(size(cand_limits,1),1)];

%Compute the matrix of possible block sizes
cand_block_size_mat = diff(ind_cand_changepoints,1,2); %candidate vectors n

%Remove candidate block size vectors whose block/s includes smaller than
%minimum number of samples
if(min_num_samples)

    [ind_row,~] = find(cand_block_size_mat < min_num_samples); ind_row = unique(ind_row);
    cand_block_size_mat(ind_row,:) = [];
end

end%endfunction
