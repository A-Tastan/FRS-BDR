%% This demo file runs  Fast and Robust Sparsity-Aware Block Diagonal
%% Representation (FRS-BDR) method for cluster enumeration
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
clear all
close all

%% Call data set
load fisheriris
X = meas.';
K = 3;

%% Define parameters
K_min = 2;
K_max = 2*K;
N_cmax = 2*(K_max-1);

%% FRS-BDR Algorithm for Unknown Number of Clusters
[W_block_hat,v_dddot_hat,W_sim_hat,bsize_vec_hat,error] = perform_FRS_BDR_algorithm_unknownK(X,K_min,K_max,N_cmax);







