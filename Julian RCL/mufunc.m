function f = mufunc(x_2,sigmavector)
% This function computes the non-linear part of the utility (mu_ijt in the Guide)

% inspired by Aviv Nevo, May 1998.

global ns vfull
[n k] = size(x_2);
% n is number of products I think. ->970 in our case; and k is 4
mu = zeros(n,ns);
for i = 1:ns
        %for every "sample", pick the k=4 random realizations of v.
    	v_i = vfull(:,(1+k*(i-1)):k*i);
 		mu(:,i) = ((x_2.*v_i)*sigmavector);
end
f = mu;
% this returns a 970x(ns) matrix of
end