% For a given nearest neighbor number p, this function computes a sparse
% affinity matrix.
%
% For details, see :
%
% [1] D. Cai, X. He and J. Han, "Document clustering using locality
% preserving indexing” IEEE Trans. Knowl. Data Eng., vol. 17, pp. 1624-
% 1637, 2005.
%
% [2] A. Taştan, M. Muma and A. M. Zoubir, "Fast and Robust Sparsity-Aware
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
%
% INPUTS:
%       dump         :(A numeric) sorted affinity matrix of size
%                     (N-N_{I}) x (N-N_{I}) if Type I outlier removal has
%                     been preferred in Step 1.1. (N x N otherwise)
%       idx          :(A numeric) sorting indexes matrix of size
%                     (N-N_{I}) x (N-N_{I}) if Type I outlier removal has
%                     been preferred in Step 1.1. (N x N otherwise)
%       p            :(A numeric) specified number of neighbors
%       N            :(A numeric) number of samples N-N_{I} if
%                     preprocessing is preferred in the previous step.
%                     (N otherwise)
%
% OUTPUTS:
%        W_pnn       :(numeric)sparse affinity matrix corresponds to p
%                     nearest neighbor of size (N-N_{I}) x (N-N_{I})
%                     if Type I outlier removal has been preferred
%                     in Step 1.1.(N x N otherwise)
%
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W_pnn=find_affinity_p_nearest_neighbor(dump,idx,p,N)

nBlock = 4000;
G = zeros(N*(p+1),3);

for m = 1:ceil(N/nBlock)
    if m == ceil(N/nBlock)

        smpIdx = (m-1)*nBlock+1:N;
        idx = idx(:,1:p+1);
        dump = -dump(:,1:p+1);

        G((m-1)*nBlock*(p+1)+1:N*(p+1),1) = repmat(smpIdx',[p+1,1]);
        G((m-1)*nBlock*(p+1)+1:N*(p+1),2) = idx(:);
        G((m-1)*nBlock*(p+1)+1:N*(p+1),3) = dump(:);

    else

        smpIdx = (m-1)*nBlock+1:m*nBlock;
        idx = idx(:,1:p+1);
        dump = -dump(:,1:p+1);

        G((m-1)*nBlock*(p+1)+1:m*nBlock*(p+1),1) = repmat(smpIdx',[p+1,1]);
        G((m-1)*nBlock*(p+1)+1:m*nBlock*(p+1),2) = idx(:);
        G((m-1)*nBlock*(p+1)+1:m*nBlock*(p+1),3) = dump(:);
    end
end

W_pnn = sparse(G(:,1),G(:,2),G(:,3),N,N);
clear G dist
W_pnn = max(W_pnn,W_pnn');
W_pnn=full(W_pnn);  W_pnn=W_pnn-diag(diag(W_pnn));%full zero diagonal affinity matrix

end







