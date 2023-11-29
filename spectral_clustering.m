% This function performs spectral clustering.
%
% For details, see:
%
% [1] U. Von Luxburg, “A tutorial on spectral clustering,”Stat. Comput.,vol.
% 17, pp. 395-416, 2007.
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
%         W              :(numeric) zero diagonal symmetric affinity matrix
%                         of size (N-N_{I}) x (N-N_{I}) if Type I outlier
%                         removal has been preferred in Step 1.1.
%                         (N x N otherwise)
%         K              :specified number of clusters
%         options_sc     :(struct) spectral clustering options
%         - options_sc.Laplacian_opt :(string) 'ShiMalik' for a normalized
%                                     Laplacian matrix according to Shi and
%                                     Malik, 'JordanWeiss' for a normalized
%                                     Laplacian matrix according to Jordan
%                                     and Weiss, 'unnormalized' (default)
%                                     for an unnormalized Laplacian matrix
%         - options_sc.eig_func_opt :(string) 'standard' or 'generalized'
%                                    (default)to identify the desired
%                                     eigen-decomposition
%         - options_sc.partition_opt :(string) 'KmedoidsTukey' for the
%                                     K-medoids algorithm based partitioning
%                                     with Tukey's distance function,
%                                     'Kmedoids' for the K-medoids algorithm
%                                     based partitioning, 'Kmeans' (default)
%                                     for the K-means based partitioning
%
% OUTPUTS:
%         C_hat                     :(numeric) estimated labels vector of
%                                    size (N-N_{o_I}) x (N-N_{o_I}) if
%                                    Type I outlier removal has been
%                                    preferred in Step 1.1.(N x N otherwise)
%
% dependencies: distfun_tuk (if 'KmedoidsTukey' is preferred)
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C_hat = spectral_clustering(W,K,options_sc)

if nargin < 3 || isempty(options_sc)
    options_sc.Laplacian_opt = 'unnormalized';
    options_sc.eig_func_opt = 'generalized';
    options_sc.partition_opt = 'Kmeans';
end


%Compute overall edge weight matrix
d = sum(W,2);
D = diag(d);

%Compute Laplacian matrix
if(strcmp(options_sc.Laplacian_opt,'ShiMalik'))

    % calculate normalized Laplacian
    L = D-W;
    L = (D^(-1))*L;

else if(strcmp(options_sc.Laplacian_opt,'JordanWeiss'))

        % calculate normalized Laplacian
        L = D-W;
        L = (D^(0.5))*L*(D^(0.5));
else

    L = D-W; % calculate unnormalized Laplacian
end
end

%Compute eigenvectors
if(strcmp(options_sc.eig_func_opt,'standard'))

    [Eigvectors,diag_Lambda] = eig(L); %standard eigen-decomposition
else

    [Eigvectors,diag_Lambda] = eig(L,D); %generalized eigen-decomposition
end

[Lambda,idx] = sort(diag(diag_Lambda)); %vector of eigenvalues in ascending order
Y = Eigvectors(:,idx(1:K));

%Apply partitioning method
if(strcmp(options_sc.partition_opt,'KmedoidsTukey'))

    [C_hat,~] = kmedoids(Y,K,'Distance',@distfun_tuk, 'Start', 'sample');

else if(strcmp(options_sc.partition_opt,'Kmedoids'))

        [C_hat,~] = kmedoids(Y,K,'Replicates',10);
else

    [C_hat,~] = kmeans(Y,K,'Replicates',10);
end%endelseif

end%endif


end%functionend
