% This function estimates the vector v.
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
%        v_ord_target_hat  :(numeric)the target vector v estimate of size
%                           (N-N_{I}) x 1 if Type I outlier removal has
%                           been preferred in Step 1.1. (N x 1 otherwise)
%        W_sim_hat_K       :(numeric)similarity coefficients matrix of size
%                           K x K. It contains similarity coefficients that
%                           the blocks are concentrated around.
%        block_startpoints :(numeric)block startpoints vector of size K x 1
%        block_endpoints   :(numeric)block endpoints vector of size K x 1
%        bsizes_vec        :(numeric) block sizes vector of size K x 1
%        K                 :(numeric) number of pieces in the linear
%                           function
%
%
% OUTPUTS:
%         v_ord_hat        :(numeric) the vector v estimate of size
%                           (N-N_{I}) x 1 if Type I outlier removal has
%                           been preferred in Step 1.1. (N x 1 otherwise)
%
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v_ord_hat = estimate_vector_v(v_ord_target_hat,W_sim_hat_K,block_startpoints,block_endpoints,bsizes_vec,K)

%Compute piece-wise linear function shifts caused by group similarity
shift = 0;
for i = 2:K
    for j = (i-1):-1:1

        shift = bsizes_vec(j)*W_sim_hat_K(i,j)+shift;
    end

    groupsim_shift(i) = shift;
    shift = 0;
end

%Estimate v_ord
for i = 1:K

    v_ord_hat(block_startpoints(i):block_endpoints(i),1) = v_ord_target_hat(block_startpoints(i):block_endpoints(i))+groupsim_shift(i);
end

end%endfunction