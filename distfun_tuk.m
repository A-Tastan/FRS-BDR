function tuk_dist = distfun_tuk(mu,X)
% mu is row-vector of centroid (1xp)
% X is matrix of data-points (Nxp)

%if size(X,1)<size(X,2)
%X = X';
%end

c = 3;
p = length(mu);
N = size(X,1);


tuk_dist_temp = zeros(N,p);

for ii = 1:N
    for jj = 1:p
        tuk_dist_temp(ii,jj) = (c^2/3)*((1-(1-(abs(mu(jj)-X(ii,jj))/c).^2).^3).*(abs(mu(jj)-X(ii,jj))<=c) + (abs(mu(jj)-X(ii,jj))>c) ); % Tukey
    end
end

tuk_dist = mean(tuk_dist_temp,2);

end