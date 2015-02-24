function robust_nonpar_set_level_svc(R,wh_cl)
% robust_nonpar_set_level_svc(R,wh_cl)
%
% Computes number of significant ROIs needed to reject omnibus hypothesis
% of no significant activation in any ROI.
%
% numc is the number of ROIs.  Optional input wh_cl (vector of which clusters) overrides the total
% number of ROIs saved in the R structure, in case you want to correct over
% a different number (i.e., a subset of ROIs of primary interest).

if nargin < 2
    numc = R.volInfo.c; 
    wh_cl = 1:numc;
else
    numc = length(wh_cl);
end

    
niter = size(R.maxt,2);

% Number of significant ROIs expected (mapwise correction for # ROIS)
% -------------------------------------------------------------------
ntoprint = min(8,numc);
pvals = 1 - binomcdf(1:ntoprint,numc,.05);
fprintf(1,'\nSet-level p-values for number of significant ROIs')
fprintf(1,'\nThis provides a correction based on the number of') 
fprintf('\n significant contiguous ROIs in your search mask')
fprintf('\n for each column in your design matrix.')
fprintf('\n significant p-values indicate more ROIs show effects')
fprintf('\n with small-volume correction than would be expected by chance.')
fprintf(1,'\n-------------------------------------------------\n')
fprintf(1,'%3.0f\t',1:ntoprint)
fprintf(1,'\nBinomial model: One-tailed:\n')
fprintf(1,'%3.4f\t',pvals)
fprintf(1,'\nBinomial model: Two-tailed:\n')
fprintf(1,'%3.4f\t',min(1,2*pvals))

n_sig_needed = min(find(pvals <= .05));
n_sig_needed2 = min(find(pvals <= .025));
fprintf(1,'\n\n# sig. ROIs needed to reject mapwise omnibus of no activation\n');
fprintf(1,'%3.0f one-tailed, %3.0f two-tailed\n',n_sig_needed,n_sig_needed2);

% nonparametric version
fprintf(1,'\nNonparametric: \n')
k = size(R.X,2);

for col = 1:k
    % t values for this regressor in this cluster
    dat = zeros(numc, niter);
    
    % get list of which ROIs activated in the Monte Carlo simulation for
    % each iteration and for this regressor
    for i = 1:length(wh_cl)
        dat(i,:) = R.maxt_by_cl{wh_cl(i)}(col,:);
    end
    thr = repmat(R.svc.tcrit(wh_cl, col), 1, niter);
    dat = sum(dat > thr);

    % get p-values for small numbers of sig. ROIs
    for i = 1:ntoprint
        pvals(i) = sum(dat >= i) ./ niter;
    end
    
    fprintf(1,'\n%s, One-tailed nonparametric p-values:\n',R.names{col})
    fprintf(1,'%3.4f\t',pvals)
    
    % get number of sig. ROIs needed
    n_sig_needed = prctile(dat,95);
    n_sig_needed2 = prctile(dat,97.5);
    fprintf(1,'\nCritical number: %3.0f one-tailed, %3.0f two-tailed\n',n_sig_needed,n_sig_needed2);
    fprintf(1,'\n')
    
end

return