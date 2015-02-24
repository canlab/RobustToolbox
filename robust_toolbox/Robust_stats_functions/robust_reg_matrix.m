function [b,t,p,sig,f,fp,fsig,stat] = robust_reg_matrix(X,Y,dochk)
% [b,t,p,sig,F,fp,fsig,stat] = robust_reg_matrix(X,dat,[do checks and verbose output (1/0)])
%
% takes a data matrix Y (voxels x obs.) and model matrix X ( and returns the robust reg
% coefficients and other things, and significance
%
% X should already contain an intercept.
% Empty X will perform a one-sample t-test
%
% F-values for full model (including intercept!) can be requested
% (slower...)
%
% Fastest:
% tic, [b,t] = robust_reg_matrix(X,Y,0); toc
%
% tor wager, july 06

% use typical robust regression right now; slower, but easier
%
% compare weighted_glmfit using robust weights w with output of robustfit
% betas are same, stes are larger for robustfit
%
% stat is a structure with: wts,se,resid, s from robust regression

% -------------------------------------------------------------------
% Setup and check
% -------------------------------------------------------------------
sig = []; fsig = [];

if nargin < 3, dochk = 1; end

[n2,v] = size(Y);
if isempty(X), X = ones(n2,1); end
[n,k] = size(X);

if dochk
    if n ~= n2, error('data and model sizes do not match.'); end

    if no_intercept(X), error('X should contain an intercept.'); end

    if no_variance(Y), error('Some Y vectors have no variability.  You must remove these before running.'); end
end

% -------------------------------------------------------------------
% Setup outputs
% -------------------------------------------------------------------
b = zeros(k,v);
t = b;
if nargout > 2
    p = b;
    f = zeros(1,v);
    fp = f;
end

if nargout > 7
    stat.wts = zeros(n,v);
    stat.se = zeros(k,v); 
    stat.s = zeros(1,v);
    stat.resid = zeros(n,v);
end

% -------------------------------------------------------------------
% Loop through data vectors
% -------------------------------------------------------------------
warning('off','stats:statrobustfit:IterationLimit');

if dochk
    str = sprintf('Running robust regression for %3.0f data vectors: Done  %03d%%',v,0); fprintf(1,str);
    updateiterations = 1:round(v ./ 100):v;
    updateperc = round(linspace(0,100,length(updateiterations)));
end

for i = 1:v

    % Robust version
    [b(:,i),stats] = robustfit(X,Y(:,i),'bisquare',[],'off');
    
    % OLS version
    % [b,dev,stats] = glmfit(X,Y(:,i),'normal','constant','off');
    
    t(:,i) = stats.t;
    if nargout > 2
        p(:,i) = stats.p;
    end

    if nargout > 4
        [f(:,i), fp(:,i)] = F_test_no_intercept(X,Y(:,i),stats.s);
    end

    if nargout > 7
        stat.wts(:,i) = stats.w;
        stat.se(:,i) = stats.se;
        stat.resid(:,i) = stats.resid;
        stat.s(1,i) = stats.s;
    end
        
    if dochk
        % print string
        update = updateiterations == i;
        if any(update), fprintf(1,'\b\b\b\b%03d%%',updateperc(update)); end
    end

end

if dochk
    erase_string(str);
    %fprintf(1,'\n');

    sig = p < .05;
end

if dochk && nargout > 4
    fsig = fp < .05;
end

warning('on','stats:statrobustfit:IterationLimit');


return

% -------------------------------------------------------------------
% Sub-functions
% -------------------------------------------------------------------

function val = no_intercept(X)

val = all(any(diff(X)));

return

function val = no_variance(Y)

val = ~all(any(diff(Y)));

return

function erase_string(str1)
fprintf(1,repmat('\b',1,length(str1))); % erase string
return

