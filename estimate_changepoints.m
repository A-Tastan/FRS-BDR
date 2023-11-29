% This function estimates the changepoints in the vector v (Step 2.1).
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
%        v      :(numeric)the vector v of size (N-N_{I}) x 1  if Type I
%                outlier removal has been preferred in Step 1.1.
%                (N x 1 otherwise)
%        N_cmax :(numeric) maximum number of changepoints
%
% OUTPUTS:
%         tau   :(numeric) changepoint locations vector of maximum size
%                N_cmax x 1
%
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tau = estimate_changepoints(v,N_cmax)

%Estimate changepoints
[tau,~] = find(ischange(v,'linear','MaxNumChanges',N_cmax));

end%endfunction