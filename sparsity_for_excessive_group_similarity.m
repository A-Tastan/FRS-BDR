% This function designs a sparse Laplacian matrix for the case of excessive
% group similarity (Step 1.3).
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
%       W               :(numeric) zero diagonal symmetric affinity matrix
%                        of size (N-N_{I}) x (N-N_{I}) if Type I outlier
%                        removal has been preferred in Step 1.1.
%                        (N x N otherwise)
%       N               :(numeric) number of samples N-N_{I} if
%                        preprocessing if Type I outlier removal has been
%                        preferred in Step 1.1.(N otherwise)
%       K_max           :(numeric) specified maximum number of blocks
%       eig_func_opt    :(string) 'standard' or 'generalized' (default)
%                        to identify the desired eigen-decomposition
%       sim_func_opt    :(string) similarity measure 'AdaptiveThresholding'
%                        or 'p-nearestNeighbor' (default)
%       options         :(struct)
%        - options.TOL             :(numeric) tolerence level for zero
%                                   eigenvalues(default=1e-3)
%        - options.ini_thr_AdapThr :(numeric) initial threshold (T) for the
%                                   adaptive thresholding (default=0.5)
%        - options.inc_thr_AdapThr :(numeric) increasement factor for the
%                                   adaptive thresholding (default=1e-3)
%        - options.ini_p           :(numeric)initial nearest neighbor value
%                                   (p) for the p-nearest neighbor graph
%                                   (default=N-2)
%        - options.dec_pnn         :(numeric) decreasement factor for the
%                                   p-nearest neighbor graph (default=1)
%        - options.min_num_samples :(numeric) minimum number of samples
%                                   (defaut~N/k_max)
%
% OUTPUTS:
%        W_sparse_hat    :(numeric)sparse affinity matrix estimate of size
%                         (N-N_{I}) x (N-N_{I}) if if Type I outlier
%                         removal has been preferred in Step 1.1.
%                         (N x N otherwise)
%        D_sparse_hat    :(numeric) estimated overall edge weight matrix
%                         (associated with W_sparse_hat) of size
%                         (N-N_{I}) x (N-N_{I}) if if Type I outlier
%                         removal has been preferred in Step 1.1.
%                         (N x N otherwise)
%        L_sparse_hat    :(numeric)the Laplacian matrix estimate(associated
%                         with W_sparse_hat) of size(N-N_{I})x(N-N_{I})
%                         if Type I outlier removal has been preferred in
%                         Step 1.1.(N x N otherwise)
%        v_sparse_hat    :(numeric)the vector v (associated with
%                         L_sparse_hat) of size(N-N_{I}) x 1 if Type I
%                         outlier removal has been preferred in Step 1.1.
%                         (N x 1 otherwise)
%
% dependencies: find_affinity_p_nearest_neighbor
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W_sparse_hat,D_sparse_hat,L_sparse_hat,v_sparse_hat] = sparsity_for_excessive_group_similarity(W,N,K_max,eig_func_opt,sim_func_opt,options)

if nargin < 6 || isempty(options)
    options.TOL = 1e-3;               %tolerence level for zero eigenvalues
    options.ini_thr_AdapThr = 0.5;    %initial threshold for the adaptive thresholding
    options.inc_thr_AdapThr = 1e-3;   %increasement factor for the adaptive thresholding
    options.ini_p = N-2;              %initial nearest neighbor value for the p-nearest neighbor graph
    options.dec_pnn = 5;              %decreasement factor for the p-nearest neighbor graph
    options.min_num_samples = round(N/K_max); %minimum number of samples in each block
end

if nargin < 5 || isempty(sim_func_opt)
    sim_func_opt = 'p-nearestNeighbor';
end

if nargin < 4 || isempty(eig_func_opt)
    eig_func_opt = 'generalized';
end

W_cand = W;
%% Sparse Laplacian matrix design based on adaptive thresholding
if(strcmp(sim_func_opt,'AdaptiveThresholding'))

    T = options.ini_thr_AdapThr;

    while(1 > T)

        W_cand(W_cand < T) = 0; %candidate affinity matrix
        D_cand = diag(sum(W_cand)); % candidate overall edge weight matrix
        L_cand = D_cand-W_cand; % candidate Laplacian matrix

        if(strcmp(eig_func_opt,'standard'))
            [~,diag_Lambda_cand] = eig(L_cand);
        else
            [~,diag_Lambda_cand] = eig(L_cand,D_cand);
        end
        Lambda_cand = sort(diag(diag_Lambda_cand)); % candidate vector of eigenvalues

        if( Lambda_cand(2) < options.TOL & isreal(Lambda_cand) & all(Lambda_cand >= 0))

            W_sparse_hat = W_cand; D_sparse_hat = D_cand; L_sparse_hat = L_cand;
            break;

        else if ( ~(Lambda_cand(2) < options.TOL) & isreal(Lambda_cand) & all(Lambda_cand >= 0))

                W_sparse_alt = W_cand; D_sparse_alt = D_cand; L_sparse_alt = L_cand; %maximum sparse alternative that provides real nonnegative eigenvalues
                clear D_cand L_cand diag_Lambda_cand Lambda_cand
                T = T + options.inc_thr_AdapThr;
        else

            clear D_cand L_cand diag_Lambda_cand Lambda_cand
            T = T + options.inc_thr_AdapThr;
        end

        end
    end%endwhile

    if(~exist('L_sparse_hat','var') & exist('L_sparse_alt','var'))

        W_sparse_hat = W_sparse_alt; D_sparse_hat = D_sparse_alt; L_sparse_hat = L_sparse_alt;
    else if~(exist('L_sparse_hat','var') || exist('L_sparse_alt','var'))

            error('Please select a smaller initial threshold (T)');
    end
    end

    %% Sparse Laplacian matrix design based on nearest neighbor graph
else

    [dump,idx] = sort(-W_cand,2); % sort each row
    p = options.ini_p;

    while(p >= options.min_num_samples)
        W_cand = find_affinity_p_nearest_neighbor(dump,idx,p,N);
        D_cand = diag(sum(W_cand)); % candidate overall edge weight matrix
        L_cand = D_cand-W_cand; % candidate Laplacian matrix

        if(strcmp(eig_func_opt,'standard'))
            [~,diag_Lambda_cand] = eig(L_cand);
        else
            [~,diag_Lambda_cand] = eig(L_cand,D_cand);
        end
        Lambda_cand = sort(diag(diag_Lambda_cand)); % candidate vector of eigenvalues

        if( Lambda_cand(2) < options.TOL & isreal(Lambda_cand) & all(Lambda_cand >= 0))

            W_sparse_hat = W_cand; D_sparse_hat = D_cand; L_sparse_hat = L_cand;
            break;

        else if ( ~(Lambda_cand(2) < options.TOL) & isreal(Lambda_cand) & all(Lambda_cand >= 0))

                W_sparse_alt = W_cand; D_sparse_alt = D_cand; L_sparse_alt = L_cand; %maximum sparse alternative that provides real nonnegative eigenvalues
                clear W_cand D_cand L_cand diag_Lambda Lambda
                p = p - options.dec_pnn;
        else

            clear W_cand D_cand L_cand diag_Lambda Lambda
            p = p - options.dec_pnn;
        end
        end
    end%endwhile

    if(~exist('L_sparse_hat','var') & exist('L_sparse_alt','var'))

        W_sparse_hat = W_sparse_alt; D_sparse_hat = D_sparse_alt; L_sparse_hat = L_sparse_alt;

    else if~(exist('L_sparse_hat','var') || exist('L_sparse_alt','var'))

            error('Please select a greater initial nearest neighbor value (p)');
    end
    end

end%endif

%Compute the vector v
v_sparse_hat = sum(triu(L_sparse_hat),2);

end%endfunction