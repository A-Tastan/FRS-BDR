% This function estimates the underlying block diagonal order based on
% the proposed similarity-based block diagonal ordering function (Step 1.2).
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
%         W          :(numeric) zero diagonal symmetric affinity matrix of
%                     size (N-N_{I}) x (N-N_{I}) if Type I outlier removal 
%                     has been preferred in Step 1.1.(N x N otherwise)
%         D          :(numeric) diagonal overall edge weight matrix of size
%                     (N-N_{I}) x (N-N_{I}) if Type I outlier removal 
%                     has been preferred in Step 1.1. (N x N otherwise)
%         N          :(numeric) number of samples valued as N-N_{I} if
%                     Type I outlier removal has been preferred in Step 1.1.
%                     (N otherwise)
%
% OUTPUT:
%         order_hat  :(numeric) vector of block diagonal order indexes of
%                     size (N-N_{I}) x 1 if Type I outlier removal has
%                     been preferred in Step 1.1. (N x 1 otherwise)
%
%
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function order_hat = sBDO(W,D,N)

d_vec = diag(D);
[~,order_hat(1)] = max(d_vec); %node correspond to maximum overall edge weight

for m = 2:N

    %Find Neighbors
    [~,neighbors] = find(W(order_hat,:));
    neighbors = setdiff(unique(neighbors),order_hat);

    if(~isempty(neighbors))

        %Order Neighbors
        [~,idx] = max(sum(W(order_hat,neighbors),1));
        order_hat(m) = neighbors(idx);

    else

        d_vecdiff = d_vec;   d_vecdiff(order_hat) = -inf;
        [~,order_hat(m)] = max(d_vecdiff);
    end

end%endfor

end%endfunction
