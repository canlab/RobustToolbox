function [r,xy,v,wmean] = weighted_corrcoef(x,w)
% [r,xy,v,wmean] = weighted_corrcoef(x,w)
%
% x is m subjects by p variables data matrix
% w is a vector of case weights
%
% r = weighted correlation coeff across columns of x
% xy = weighted covariance matrix
% d = weighted variance estimates for each column of x
% wmean = weighted mean of each column of X
%
% Weights are as from IRLS robustfit, stat.w
% Assumed to have mean of 1
% To get from weights of inv var, e.g., subtract mean weight and add 1
%
% tor wager

if abs(mean(w)) - 1 > eps,
    error('weighted corrcoef: Weights should have mean of 1');
end


W = diag(w);

  % weighted covariance
  [m,n] = size(x);
  xc = x - repmat(sum(x)/m,m,1);  % Remove mean
  xy = xc' * W * xc / (m-1);          % sum of squared dev. / m-1
  
  xy = 0.5*(xy+xy');              % Remove rounding error

  
  %xy = cov(x)
  v = diag(xy);
  r = xy./sqrt(v*v');

  v = v';                       % make a row vector
  if nargout > 3
      wmean = sum(W*x) ./ m; 
  end
  
  return
  
  % statrobustfit
  %[n,p] = size(X);
  %[Q,R,perm] = qr(X,0);
  %tol = abs(R(1)) * max(n,p) * eps(class(R));
  %xrank = sum(abs(diag(R)) > tol)

  %rank(X)
  %dfe = n - xrank
  