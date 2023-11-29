% This function performs a plane piece-wise linear fit.
%
% For details, see :
%
% [1] X. Yang, H. Yang, F. Zhang, L. Zhang, X. Fan, Q. Ye and L. Fu,
% “Piecewise linear regression based on plane clustering,” IEEE Access,
%  vol. 7, pp. 29845-29855, 2019.
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
% INPUTS:
%       v_i        :(A numeric) vector that contains vector v components
%                   associated with ith block. The dimension of v_i is
%                   N_{r_i} x 1
%       m          :(A_numeric) vector that contains sample points
%                   corresponds to v_i. The dimension of m is N_{r_i} x 1.
%
% OUTPUTS:
%       v_i_hat    :(A numeric) vector v estimate associated with ith
%                   block. The dimension is N_{r_i} x 1.
%       reg_coeff  :(A numeric) contains the estimated regression
%                   coefficient
%
% version     : 29/11/2023
% authors     : Aylin Tastan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_i_hat,reg_coeff] = plane_based_piecewise_linear_fit(v_i,m)

%Form the sample matrix
Upsilon_i = [m;v_i.'].';

%Compute covariance matrix
Phi_i = cov(Upsilon_i);    %2x2 matrix

%Compute mean vector
Mu_i = mean(Upsilon_i).';  %2x1 vector

%Estimate the normal vector and the bias using eigen-decomposition of Phi
[Eigenvec_Phi,Lambdas_Phi] = eig(Phi_i);
[~, ind_Lambdas] = sort(diag(Lambdas_Phi));
vartheta_i = Eigenvec_Phi(:,ind_Lambdas(1));
b_i = -vartheta_i.'*Mu_i;

%Linear regression
v_i_hat = -(vartheta_i(1)*m+b_i)/vartheta_i(2);
reg_coeff = -(vartheta_i(1)/vartheta_i(2));

end