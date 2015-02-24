function [t,p,maxt,maxt_by_cl, primary_uncor] = robust_pooled_weight_core(X,Y,b,resid,wts,niter,wh_intercept,clindx,varargin)
%
%  [t,p,maxt,maxt_by_cl, primary_uncor] = robust_pooled_weight_core(X,Y,[],[],[],2000,3,R.volInfo.clindx,[existing maxt,maxtbycl]);
%
%  [R.correct.weightedt,R.correct.weightedp,R.maxt,R.maxt_by_cl] = ...
%  robust_pooled_weight_core(X,Y,R.correct.b,R.correct.resid, ...
%  R.correct.wts,2000,R.wh_intercept,R.volInfo.clindx);
%
% clindx = R.volInfo.clindx
%
% Updated: Jan 2008, by Tor Wager, to add primary uncorrected thresholds

% --------------------------------------
% * Setup, and run robust reg if needed
% --------------------------------------
[n,v] = size(Y);
k = size(X,2);

t = zeros(k,v);
p = t;

% set up chunks of iterations that will fit in memory
maxmemsize = 2^28;
memneeded = prod([niter v k]) * 8;
chunks = ceil(memneeded./maxmemsize);   % how many "chunks" of iterations needed to fill the request
maxiter = floor(maxmemsize ./ (8*v*k)); % max number of iterations we can run
if maxiter == 0, error('Not enough memory. You must increase maxmemsize in robust_pooled_weight_core.'); end
maxiter = min(maxiter,niter);

iterst = 1:maxiter:niter;               % starting indices of iterations
iteren = [iterst(2:end)-1 niter];       % ending indices of iterations for each chunk
iterations = iteren - iterst + 1;       % number of iterations for each chunk
fprintf(1,'%3.0f total iterations.  %3.0f chunks, each of %3.0f iterations.\n',niter,chunks,maxiter)
fprintf(1,'Adding iterations to analysis : %3.0f through %3.0f\n',iterst(1),iteren(end));

if isempty(b) || isempty(resid) || isempty(wts)
    [b,t,p,sig,f,fp,fsig,stat] = robust_reg_matrix(X,Y);

    resid = stat.resid;
    wts = stat.wts;

end

if isempty(clindx), clindx = ones(1,v);  end

% pooling se is a bad idea w/multiple columns, because of predictors
% account for lots of variance, this will produce low se, but it'll be
% averaged with surrounding noise.

% Set up primary uncorrected thresholds
primary_uncor.primary_pthr = [.0005 .001 .0025 .005 .01 .025 .05];
primary_uncor.primary_pthr_descrip = 'one-tailed primary p-value thresholds';

primary_tthr = tinv(1 - primary_uncor.primary_pthr, size(X, 1) - size(X, 2));  % one-tailed p-vals
primary_uncor.primary_tthr = primary_tthr;


% --------------------------------------
% * Do iterations
% --------------------------------------


nclust = max(clindx);
maxt_by_cl = cell(1,nclust);

% replace old estimates with pooled-weight least squares for output
fprintf(1,'\nPooling weights: crunching by robust_pooled_weight_core.m\n');

% break into "chunks"

for thischunk = 1:chunks

    cstr = sprintf('Chunk %3.0f. ',thischunk); fprintf(1,cstr);

    % Run a set of iterations (a chunk)
    % -------------------------------------------------------------------
    [t,p,tnull] = pooled_weight_lsq(X, b, resid, wts, wh_intercept, iterations(thischunk));

    % Save counts at uncorrected primary threshold for these iterations
    % ------------------------------------------------------------------- 
    % tnull is voxels x iterations x regressors
    for pt = 1:length(primary_tthr)
        sumt = squeeze(sum(tnull > primary_tthr(pt)))';  % a matrix of counts, regressors x iterations
    
        primary_uncor.null_sig{pt}(:,iterst(thischunk):iteren(thischunk)) = sumt;
    end
    
    % Save maxima for these iterations
    % -------------------------------------------------------------------
    % mapwise max
    maxtnull = squeeze(max(tnull));
    
    try
        maxt(:,iterst(thischunk):iteren(thischunk)) = maxtnull;
    catch
        disp('transpose issue');
        maxt(:,iterst(thischunk):iteren(thischunk)) = maxtnull';
    end

    % tnull: voxels x iterations (x regressors)

    for c = 1:max(clindx)

        wh = clindx == c;
        
        try
            maxt_by_cl{c}(:,iterst(thischunk):iteren(thischunk)) = squeeze(max(tnull(wh,:,:),[],1));  %tmax(:,wh);
        catch
            disp('transpose issue');
            maxt_by_cl{c}(:,iterst(thischunk):iteren(thischunk)) = squeeze(max(tnull(wh,:,:),[],1))';
        end
    end

    clear tnull

    erase_string(cstr)

end   % chunk

return


% -------------------------------------------------------------------
%
%
% Sub-functions
%
%
% -------------------------------------------------------------------

% -------------------------------------------------------------------
% Run a set of iterations (a chunk)
% -------------------------------------------------------------------

function [t,p,tnull] = pooled_weight_lsq(X,b,resid,wts,wh_intercept,niter)
% weighted least squares for each effect
% pooling weights over all voxels

sstr = sprintf('Allocating memory. '); fprintf(1,sstr);


[n,v] = size(resid);
k = size(X,2);

tnull = zeros(v,niter,k);

dfe = n - k;

t = zeros(k,v);
p = zeros(k,v);
%tmax = zeros(k,v);
%nulltmax = zeros(k,niter);

erase_string(sstr)

for wh_reg = 1:k

    kstr = sprintf('Reg. %02d: ',wh_reg); fprintf(1,kstr);



    % --------------------------------------
    % * Y (data) adjusted for effects of
    %  no (current) interest
    % --------------------------------------
    Yadj = X(:,wh_reg) * b(wh_reg,:) + resid;

    % --------------------------------------
    % * Get betas, t, std errs based on pooled weights
    %   and robust standard errors
    % --------------------------------------

    [t(wh_reg,:),p(wh_reg,:)] = weighted_glmfit_reddf(X(:,wh_reg),Yadj,mean(wts,2),dfe);

    % --------------------------------------
    % * permutations
    % --------------------------------------
    str = sprintf('Iteration %05d',0); fprintf(1,str);

    for i = 1:niter

        if wh_reg == wh_intercept
            % permute signs
            signmtx = repmat(sign(randn(n,1)),1,v);
            Yp = Yadj .* signmtx;
        else
            % permute order of rows
            Yp = Yadj(randperm(n)',:);
        end

        % tnull: voxels x iterations (x regressors)
        tnull(:,i,wh_reg) = weighted_glmfit_reddf(X(:,wh_reg),Yp,mean(wts,2),dfe)';
        %nulltmax(wh_reg,i) = max(tnull);

        if mod(i,100) == 0
            fprintf(1,'\b\b\b\b\b%05d',i)
        end

    end

    %tmax(wh_reg,:) = max(tnull);

    erase_string(str)
    erase_string(kstr)
end

%%****allocate tmax
% squeeze down to k x voxels
%tmax = squeeze(max(tnull))';
%clear tnull



return



function [t,p,b,v] = weighted_glmfit_reddf(X,Y,wts,dfe)

n = size(Y,2);


W = diag(wts);                    % Weight matrix
invxwx = inv(X'*W*X);
bform = invxwx * X'* W;         % beta-forming matrix.  hat = X * bform
%
% % rows are columns of design (X), cols are Y variables
b = bform*Y;
k = size(b,1);

e = Y - X * b;         % residuals

% % standard error
MSE = zeros(1,n);
for i=1:n, MSE(i) = e(:,i)'*W*e(:,i); end, MSE = MSE/dfe;

v = repmat(diag(invxwx),1,n) .* repmat(MSE,k,1);

t = b ./ sqrt(v);

if nargout > 1, p = 2 * ( 1 - tcdf(abs(t),dfe) ); end

return


function erase_string(str1)
fprintf(1,repmat('\b',1,length(str1))); % erase string
return
