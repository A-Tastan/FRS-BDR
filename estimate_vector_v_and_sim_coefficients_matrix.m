% This function estimates the vector v and the matrix of similarity
% coefficients (Step 2.2).
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
% INPUTS:
%       L_ord               :(numeric) block diagonal ordered Laplacian
%                            matrix of size (N-N_{I}) x (N-N_{I}) if Type I
%                            outlier removal has been preferred in Step 1.1.
%                            (N x N otherwise)
%       v_ord               :(numeric) the vector v of size(N-N_{I}) x 1 if
%                            Type I outlier removal has been preferred in
%                            Step 1.1.(N x 1 otherwise)
%       cand_block_size_mat :(numeric) matrix of candidate block sizes. It
%                             has maximum size of xi x K (default).
%       K                   :(numeric) the number of pieces in the linear
%                            function
%       fit_opt             :(string)piece-wise linear function fitting
%                            option; 'Mreg' for M-regression,
%                                    'ls' for least squares,
%                                    'plane' for the plane-based(default)
%
% OUTPUTS:
%         v_ord_hat_K     :(numeric) the vector v estimate of size
%                          (N-N_{I}) x 1 if Type I outlier removal has been
%                          preferred in Step 1.1. (N x 1 otherwise)
%         W_sim_hat_K     :(numeric) similarity coefficients matrix of size
%                           K x K. It contains similarity coefficients that
%                           the blocks are concentrated around.
%         bsize_vec_hat_K :(numeric) estimated block size vector of size
%                           K x 1
%         error_K         :(numeric) estimation error associated to K
%
% dependencies: estimate_sim_coefficients_matrix,estimate_vector_v
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[v_ord_hat_K,W_sim_hat_K,bsize_vec_hat_K,error_K] = estimate_vector_v_and_sim_coefficients_matrix(L_ord,v_ord,pblock_size_mat,K,fit_opt)

if nargin < 5 || isempty(fit_opt)
    fit_opt = 'plane';
end

%Estimate vector v
error_K = inf;
v_ord_hat_K = [];
W_sim_hat_K = [];
bsize_vec_hat_K = [];
for cand = 1:size(pblock_size_mat,1)

    %Define start and endpoints of the blocks
    block_endpoints = (cumsum(pblock_size_mat(cand,:))).';
    block_startpoints = [1;block_endpoints(1:end-1)+1];

    %Estimating matrix of similarity coefficients W_sim
    [W_sim_hat_cand,v_ord_target_hat_cand] = estimate_sim_coefficients_matrix(L_ord,K,block_startpoints,block_endpoints,pblock_size_mat(cand,:).',fit_opt);

    %Vector v estimate
    v_ord_hat_cand = estimate_vector_v(v_ord_target_hat_cand,W_sim_hat_cand,block_startpoints,block_endpoints,pblock_size_mat(cand,:).',K);

    %Calculate estimation error
    error_cand = norm(v_ord-v_ord_hat_cand);

    %Update estimate
    if(error_K > error_cand & isequal((max(W_sim_hat_cand)).',diag(W_sim_hat_cand)))

        error_K = error_cand;
        v_ord_hat_K = v_ord_hat_cand;
        W_sim_hat_K = W_sim_hat_cand;
        bsize_vec_hat_K = pblock_size_mat(cand,:).';
    end
    clear error_cand v_ord_hat_cand W_sim_hat_cand v_ord_target_hat_cand

end%endfor

end%endfunction