% This function performs FRS-BDR algorithm for unknown number of clusters.
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
% !NOTE :
% This code requires an additional function 'Mreg' and its all dependencies
% which is available in :
%
% https://github.com/RobustSP/toolbox/tree/master/codes
%
% INPUTS:
%       X               :(numeric) a data matrix of size M x N
%       K_min           :(numeric) specified minimum number of blocks
%       K_max           :(numeric) specified maximum number of blocks
%       N_cmax          :(numeric) maximum number of changepoints
%       min_num_samples :(numeric) minimum number of samples in the blocks
%       order_opt       :(string) desired ordering algorithm; 'RCM' for the
%                        reverse Cuthill-McKee algorithm, 'sBDO' for the
%                        similarity based block diagonal ordering algorithm
%                        (default='sBDO')
%       sim_func_opt    :(string) similarity measure 'AdaptiveThresholding'
%                         or 'p-nearestNeighbor' (default)
%       fit_opt         :(string)piece-wise linear function fitting option;
%                        'Mreg' for M-regression,
%                        'ls' for least squares,
%                        'plane' for the plane-based(default)
%       eig_func_opt    :(string) 'standard' or 'generalized' (default)
%                         to identify the desired eigen-decomposition
%       removal_opt     :(string) 'no' or 'yes' (default). Set as 'yes'
%                         for Type I outlier removal
%       TOL             :(numeric) tolerence level for zero eigenvalues
%                        (default=1e-3)
%
% OUTPUTS:
%       W_block_hat     :(numeric) zero diagonal block affinity matrix of
%                         size (N-N_{o_I}) x (N-N_{o_I}) if Type I outlier
%                         removal has been preferred in Step 1.1.
%                         (N x N otherwise)
%       v_dddot_hat     :(numeric) the vector v estimate of size 
%                        (N-N_{I}) x 1 if Type I outlier removal has been
%                         preferred in Step 1.1. (N x 1 otherwise)
%       W_sim_hat       :(numeric) similarity coefficients matrix of size
%                         K_hat x K_hat. It contains similarity coefficients that
%                         the blocks are concentrated around.
%       bsize_vec_hat   :(numeric) estimated block size vector of size
%                         K_hat x 1
%       error           :(numeric) estimation error associated to K_hat
%
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W_block_hat,v_dddot_hat,W_sim_hat,bsize_vec_hat,error] = perform_FRS_BDR_algorithm_unknownK(X,K_min,K_max,N_cmax,min_num_samples,order_opt,sim_func_opt,fit_opt,eig_func_opt,removal_opt,TOL)

if nargin < 11 || isempty(TOL)
    TOL = 1e-3;
end

if nargin < 10 || isempty(removal_opt)
    removal_opt = 'yes';
end

if nargin < 9 || isempty(eig_func_opt)
    eig_func_opt = 'generalized';
end

if nargin < 8 || isempty(fit_opt)
    fit_opt = 'plane';
end

if nargin < 7 || isempty(sim_func_opt)
    sim_func_opt = 'p-nearestNeighbor';
end

if nargin < 6 || isempty(order_opt)
    order_opt = 'sBDO';
end

if nargin < 5 || isempty(min_num_samples)
    min_num_samples = round(size(X,2)/K_max);
end


%% STEP 1 : ENHANCING BLOCK DIAGONAL STRUCTURE
%% Step 1.1 : Initialize Parameters and Type I Outlier Removal
[X_dot,W_dot,D_dot,L_dot,N_of,ind_typeIout,sparsity] = initialization_and_type_I_outlier_removal(X,eig_func_opt,removal_opt,TOL); %N_of = N-N_{I}

%% Step 1.2 : Estimating Block Diagonal Order based on sBDO or RCM
[W_ddot,D_ddot,L_ddot,v_ddot,order_hat] = estimate_block_diagonal_order(W_dot,D_dot,L_dot,N_of,order_opt);

%% Step 1.3 : Sparsity for Excessive Group Similarity
if~(sparsity)
    try

        [W_dddot,D_dddot,L_dddot,v_dddot] = sparsity_for_excessive_group_similarity(W_ddot,N_of,K_max,eig_func_opt,sim_func_opt);
    end
end

if(isempty(W_dddot))

    W_dddot=W_ddot; D_dddot=D_ddot; L_dddot=L_ddot; v_dddot=v_ddot;
end

%% STEP 2 : ESTIMATE VECTOR v
error = inf;
tau = estimate_changepoints(v_dddot,N_cmax); %Compute vector tau containing estimated changepoint locations

for K_cand = K_min:K_max

    %% Step 2.1 : Estimating Candidate Block Sizes
    cand_block_size_mat = compute_candidate_block_sizes(tau,N_of,K_cand,min_num_samples);

    %% Step 2.2 : Estimating W_sim
    if(~isempty(cand_block_size_mat))
        [v_dddot_hat_K_cand,W_sim_hat_K_cand,bsize_vec_hat_K_cand,error_K_cand] = estimate_vector_v_and_sim_coefficients_matrix(L_dddot,v_dddot,cand_block_size_mat,K_cand,fit_opt);
    else
        error_K_cand = inf;
        clear cand_block_size_mat
    end

    %Update estimate
    if(error > error_K_cand )
        K_hat = K_cand; error = error_K_cand;  W_sim_hat = W_sim_hat_K_cand; v_dddot_hat = v_dddot_hat_K_cand; bsize_vec_hat = bsize_vec_hat_K_cand;
    end

    clear error_K W_sim_hat_K v_dddot_hat_K bsize_vec_hat_K

end%endfor

% Compute Block Diagonal Affinity Matrix
if (~isempty(W_sim_hat))

    W_block_hat = compute_block_diagonal_affinity_matrix(W_dddot,W_sim_hat,bsize_vec_hat,K_hat);
else

    W_block_hat = W_dddot;
end

end

