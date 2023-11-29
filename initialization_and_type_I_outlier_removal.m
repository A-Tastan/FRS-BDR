% This function computes initial estimate affinity matrix and performs
% Type I outlier removal (Step 1.1).
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
%         X                 :(numeric) data matrix of size M x N
%         eig_func_opt      :(string) 'standard' or 'generalized' (default)
%                            to identify the desired eigen-decomposition
%         removal_opt       :(string) 'no' or 'yes' (default). Set as 'yes'
%                            for Type I outlier removal
%         TOL               :(numeric) tolerence level for zero eigenvalues
%                            (default=1e-3)
%
% OUTPUTS:
%         X                 :(numeric) data matrix of size M x (N-N_{I})
%                            if Type I outlier removal is preferred
%                            (M x N otherwise)
%         W                 :(numeric) zero diagonal symmetric affinity
%                            matrix of size (N-N_{I}) x (N-N_{I}) if
%                            Type I outlier removal is preferred.
%                            (N x N otherwise)
%         D                 :(numeric) diagonal overall edge weight matrix
%                            of size (N-N_{I}) x (N-N_{I}) if
%                            Type I outlier removal is preferred.
%                            (N x N otherwise)
%         L                 :(numeric) Laplacian matrix of size
%                            (N-N_{I}) x (N-N_{I}) if Type I outlier removal
%                            is preferred. (N x N otherwise)
%         N                 :(numeric) number of samples valued as N-N_{I}
%                            if Type I outlier removal is preferred.
%                            (N otherwise)
%         ind_typeIout      :(numeric) Type I outlier index vector of size
%                             N_{I} x 1
%         sparsity          :(logical) 1 for a sparse matrix zero otherwise
%
%
% version: 29/11/2023
% authors: Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,W,D,L,N,ind_typeIout,sparsity] = initialization_and_type_I_outlier_removal(X,eig_func_opt,removal_opt,TOL)

if nargin < 4 || isempty(TOL)
    TOL = 1e-3;
end

if nargin < 3 || isempty(removal_opt)
    removal_opt = 'yes';
end

if nargin < 2 || isempty(eig_func_opt)
    eig_func_opt = 'generalized';
end


% Initial parameters
X = normalize(X,1,'norm',2); %Normalize columns
W = X.'*X; %affinity matrix of size N x N
W = W-diag(diag(W));  %zero diagonal affinity matrix of size N x N
D = diag(sum(W)); %Diagonal weight matrix of size N x N

%Type I outlier removal
ind_typeIout = [];
N = size(X,2);
if(strcmp(removal_opt,'yes'))

    ind_typeIout=find(ismember(W,zeros(N,1).','rows')).';

    if~(isempty(ind_typeIout))
        X(:,ind_typeIout) = [];
        W(ind_typeIout,:) = []; W(:,ind_typeIout) = [];
        D = diag(sum(W));
        N = size(X,2);
    end

end
L = D - W; %Laplacian matrix of size N x N

%Computing vector of eigenvalues
if(strcmp(eig_func_opt,'standard'))

    [~,diag_Lambda] = eig(L); %standard eigen-decomposition
else

    [~,diag_Lambda] = eig(L,D); %generalized eigen-decomposition
end

[Lambda,~] = sort(diag(diag_Lambda)); %vector of eigenvalues in ascending order

%The Laplacian matrix is sparse if and only if it has at least two eigenvalues decrease to zero
if(Lambda(2)< TOL & isreal(Lambda) & all(Lambda >= 0)) %eigenvalues are real and non-zero

    sparsity = 1;
else

    sparsity = 0;
end

end%endfunction