% This function estimates the underlying block diagonal order based on
% desired function (Step 1.2).
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
%       W         :(numeric) zero diagonal symmetric affinity matrix of
%                  size (N-N_{I}) x (N-N_{I}) if Type I outlier removal has
%                  been preferred in Step 1.1.(N x N otherwise)
%       D         :(numeric) diagonal overall edge weight matrix of size
%                  (N-N_{I}) x (N-N_{I}) if Type I outlier removal has been
%                  preferred in Step 1.1. (N x N otherwise)
%       L         :(numeric)Laplacian matrix of size(N-N_{I})x(N-N_{I})
%                  if Type I outlier removal has been preferred in Step 1.1.
%                  (N x N otherwise)
%       N         :(numeric)number of samples valued as N-N_{I} if
%                  Type I outlier removal has been preferred in Step 1.1.
%                  (N otherwise)
%       order_opt :(string) desired ordering algorithm; 'RCM' for the
%                  reverse Cuthill-McKee algorithm, 'sBDO' for the
%                  similarity based block diagonal ordering algorithm
%                  (default='sBDO')
%
% OUTPUTS:
%        W_ord     :(numeric) block diagonal ordered zero diagonal symmetric
%                   affinity matrix of size (N-N_{I}) x (N-N_{I}) if
%                   Type I outlier removal has been preferred in Step 1.1.
%                   (N x N otherwise)
%        D_ord     :(numeric) block diagonal ordered diagonal overall edge
%                   weight matrix of size (N-N_{I}) x (N-N_{I}) if
%                   Type I outlier removal has been preferred in Step 1.1.
%                   (N x N otherwise)
%        L_ord     :(numeric) block diagonal ordered Laplacian matrix of
%                   size (N-N_{I}) x (N-N_{I}) if Type I outlier removal
%                   has been preferred in Step 1.1.(N x N otherwise)
%        v_ord     :(numeric) the vector v of size(N-N_{I}) x 1 if
%                   Type I outlier removal has been preferred in Step 1.1.
%                   (N x 1 otherwise)
%        order_hat :(numeric) block diagonal order indexes vector of
%                   size (N-N_{I}) x 1 if Type I outlier removal has been
%                   preferred in Step 1.1.(N x 1 otherwise)
%
% dependencies: sBDO
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W_ord,D_ord,L_ord,v_ord,order_hat] = estimate_block_diagonal_order(W,D,L,N,order_opt)

if nargin < 5 || isempty(order_opt)
    order_opt = 'sBDO';
end

%Estimate order
if(strcmp(order_opt,'sBDO'))

    order_hat = sBDO(W,D,N);

else if(strcmp(order_opt,'RCM'))

        order_hat = symrcm(L);
end
end

% Re-order initial parameters
W_ord = W(order_hat,order_hat);
D_ord = D(order_hat,order_hat);
L_ord = L(order_hat,order_hat);
v_ord = sum(triu(L_ord),2);
order_hat = order_hat.';

end%endfunction