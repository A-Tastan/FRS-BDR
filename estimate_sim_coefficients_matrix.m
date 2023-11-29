% This function estimates the matrix of similarity coefficients (Step 2.2.1
% and Step 2.2.2)
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
%
% INPUTS:
%       L_ord             :(numeric)block diagonal ordered Laplacian matrix
%                          of size (N-N_{I}) x (N-N_{I}) if Type I outlier
%                          removal has been preferred in Step 1.1.
%                          (N x N otherwise)
%       K                 :(numeric) the number of pieces in the linear
%                          function
%       block_startpoints :(numeric)block startpoints vector of size K x 1
%       block_endpoints   :(numeric)block endpoints vector of size K x 1
%       bsizes_vec        :(numeric)block sizes vector of size k x 1
%       fit_opt           :(string) piece-wise linear function fitting
%                          option; 'Mreg' for M-regression,
%                                  'ls' for least squares,
%                                  'plane' for the plane-based(default)
%
% OUTPUTS:
%        v_ord_target_hat :(numeric)the target vector v estimate of size
%                          (N-N_{I}) x 1 if Type I outlier removal has been
%                          preferred in Step 1.1. (N x 1 otherwise)
%        W_sim_hat_K      :(numeric)similarity coefficients matrix of size
%                          K x K. It contains similarity coefficients that
%                          the blocks are concentrated around.
%
% dependencies: Mreg, plane_based_piecewise_linear_fit
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W_sim_hat_K,v_ord_target_hat] = estimate_sim_coefficients_matrix(L_ord,K,block_startpoints,block_endpoints,bsizes_vec,fit_opt)

if nargin < 6 || isempty(fit_opt)
    fit_opt = 'plane';
end

%% Estimate similarity coefficients within the blocks
for i = 1:K

    %Target ith vector v piece
    v_ord_target(block_startpoints(i):block_endpoints(i),1) = sum(tril(L_ord(block_startpoints(i):block_endpoints(i),block_startpoints(i):block_endpoints(i))*(-1),-1),2);

    %Estimating ith target linear piece
    if (strcmp(fit_opt,'Mreg'))

        [v_ord_target_hat(block_startpoints(i):block_endpoints(i),1),W_sim_hat_K(i,i)] = Mreg(v_ord_target(block_startpoints(i):block_endpoints(i),1),[1:bsizes_vec(i)].');

    else if(strcmp(fit_opt,'ls'))

            W_sim_hat_K(i,i) = mldivide([1:bsizes_vec(i)].',v_ord_target(block_startpoints(i):block_endpoints(i),1));
            v_ord_target_hat(block_startpoints(i):block_endpoints(i),1) = [1:bsizes_vec(i)].'*W_sim_hat_K(i,i);
    else

        [v_ord_target_hat(block_startpoints(i):block_endpoints(i),1),W_sim_hat_K(i,i)] = plane_based_piecewise_linear_fit(v_ord_target(block_startpoints(i):block_endpoints(i),1),[1:bsizes_vec(i)]);

    end%endelseif

    end%endif

end%endfor

%% Estimate similarity coefficients between different blocks
for i = 2:K
    for j = 1:(i-1)

        v_ord_shift = v_ord_target(block_startpoints(i):block_endpoints(i),1)+sum(L_ord(block_startpoints(i):block_endpoints(i),block_startpoints(j):block_endpoints(j))*(-1),2);
        W_sim_hat_K(i,j) = median(v_ord_shift-v_ord_target_hat(block_startpoints(i):block_endpoints(i),1))/bsizes_vec(j); %Intercept/n_size(i)

    end
    clear v_ord_shift

end%endfor

W_sim_hat_K = W_sim_hat_K+(W_sim_hat_K.' - diag(diag(W_sim_hat_K))); %Symmmetrize similarity coefficients matrix

end%endfunction