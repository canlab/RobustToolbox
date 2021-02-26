function OUT = robfit_parcelwise(imgs)
%
% OUT = robfit_parcelwise(imgs)
% 
% Example:
% imgs = load_image_set('emotionreg');
% 
% Xinfo = importdata(which('Wager_2008_emotionreg_behavioral_data.txt'));
% Xnames = Xinfo.textdata(2);
% imgs.X = zscore(Xinfo.data(:, 2));
% 
% % Voxel-wise 
% out = regress(imgs);
% t = threshold(out.t, .05, 'fdr');
% 
% o2 = canlab_results_fmridisplay(get_wh_image(t, 2), 'montagetype', 'full', 'noverbose');
% o2 = removeblobs(o2);
% o2 = addblobs(o2, get_wh_image(t, 2));
% 
% o2 = addblobs(o2, get_wh_image(t_obj, 1));


% --------------------------------------------
% Defaults and inputs
% --------------------------------------------
names = {}; 



% --------------------------------------------
% set up images and parcels
% --------------------------------------------

% Checks for required functions
if isempty(which('mafdr')), error('Sorry, you need the Matlab bioinformatics toolbox on your Matlab path to use the Storey 2002 mafdr function called here.'); end



% Load parcel atlas and initalize object
b = brainpathway();

% Make image space match atlas space
imgs = resample_space(imgs, b.region_atlas);

% Trigger extraction and averaging of parcel-wise data
b.voxel_dat = imgs.dat;

%%
datmatrix = double(b.region_dat);
[t, v] = size(datmatrix);

% --------------------------------------------
% set up and check design
% --------------------------------------------
X = imgs.X;

if isempty(X)
    X = ones(t, 1);  % for robust reg
else
    % Make sure intercept is last column
    X = intercept(X, 'end');
end

% Check for non-centered predictors (excluding the intercept)
if any(abs(mean(X)) > (100 * eps) & any(abs(X - mean(X)) > 100 * eps, 1))
    warning('SOME PREDICTORS ARE NOT CENTERED - THE INTERCEPT MAP WILL NOT BE INTERPRETABLE AS THE GROUP MEAN');
    pause(2)
end

k = size(X, 2); % number of maps to estimate - one per regressor

for i = 1:k-1, names{i} = sprintf('Predictor %d', i); end, names{k} = 'Intercept (Group avg)';

% --------------------------------------------
% Initialize output arrays
% --------------------------------------------
% Parcels are rows, images are columns

[nsubjects, maskvol, dfe] = deal(zeros(v, 1) .* NaN);     % save numbers of missing values

weights = zeros(v, t) .* NaN;                        % robust reg weights

[betas, tscores, pvalues] = deal(zeros(v, k) .* NaN);

% --------------------------------------------
% Run robust reg for each parcel
% --------------------------------------------

for i = 1:v
    
    % Remove NaNs and run only if we have enough observations
    [wasnan, Xvox, tvox] = nanremove(X, datmatrix(:, i));
    
    isok = sum(~wasnan) > 4;
    
    if isok
        % OK to run
        [bb,stats]=robustfit(Xvox, tvox, 'bisquare', [], 'off');
        
        w = naninsert(wasnan, stats.w);
        
    else
        % empty voxel: not enough data
        disp('THIS CODE SHOULD NEVER RUN B/C EXCLUDING BAD VOXELS ABOVE');
        bb = NaN .* ones(k, 1);
        stats = struct('t', bb, 'p', ones(k, 1));
        w = NaN .* ones(N, 1);
        
    end
    
    % Redistribute to maps
        betas(i, :) = bb';
        tscores(i, :) = stats.t';
        pvalues(i, :) = stats.p';
        
        dfe(i, 1) = stats.dfe;
        nsubjects(i, 1) = sum(~wasnan);
        maskvol(i, 1) = isok;
        
        weights(i, :) = w;
    
end % loop through nodes

OUT = struct('betas', betas, 'tscores', tscores, 'pvalues', pvalues, 'nsubjects', nsubjects, 'maskvol', maskvol, 'weights', weights, 'dfe', dfe);

%% --------------------------------------------
% FDR correction
% --------------------------------------------
for i = 1:k
    % for each map
    FDRq(:, i) = mafdr(pvalues(:, i));
    
    % P-threshold for FDR q < 0.05 for each map
    pthr(i) = max(pvalues(FDRq(:, i) < 0.05, i));
end

sig_q05 = FDRq < 0.05;

% Other interesting metrics to save

OUT.pthr_FDRq05 = pthr;
OUT.pthr_FDRq05_descrip = 'Highest P-value significant at FDR q < 0.05; threshold at p <= pthr';
OUT.sig_q05 = sig_q05;

maxn = max(OUT.nsubjects);
OUT.cohens_d_fdr05 = tinv(1 - OUT.pthr_FDRq05, maxn - k) ./ sqrt(maxn - k);
OUT.cohens_d_fdr05_descrip = 'Min Cohen''s d detectable at FDR q < 0.05';

% sum(sig_q05)  Sig parcels


% --------------------------------------------
% Transform back to voxelwise output maps
% --------------------------------------------

OUT.t_obj = parcel_stats2statistic_image(b.region_atlas, tscores, pvalues, dfe, sig_q05);


% --------------------------------------------
% Create mask (could do this from nsubjects later)
% --------------------------------------------

mask = get_wh_image(OUT.t_obj, 1);
mask = threshold(mask, 1 - eps, 'unc', 'noverbose');
mask = remove_empty(mask);
mask.dat = ones(size(mask.dat));
mask.p = .001 * ones(size(mask.dat)); % kludge to avoid region() using p-values
OUT.mask = mask;



% --------------------------------------------
% Print report & summary results table
% --------------------------------------------

resultstable = table;
resultstable.Properties.Description = 'Parcel-wise robust regression';
resultstable.maxT = max(OUT.tscores)';
resultstable.minP = min(OUT.pvalues)';


resultstable.sig05 = sum(OUT.pvalues < 0.05)';
resultstable.sig005 = sum(OUT.pvalues < 0.005)';
resultstable.sig001 = sum(OUT.pvalues < 0.001)';
resultstable.sigFDR05 = sum(OUT.sig_q05)';
resultstable.p_thr_FDR05 = OUT.pthr_FDRq05';
resultstable.min_d_FDR05 = OUT.cohens_d_fdr05';

resultstable.Properties.RowNames = names;
OUT.resultstable = resultstable;


% Print to screen

dashes = '__________________________________________________________________';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

printhdr(resultstable.Properties.Description);
disp(resultstable);
disp(' ')
disp('sig*: Significant parcels at given threshold (p < 0.05 two-tailed, q < 0.05 FDR, etc.)');
disp('p_thr_FDR05: P-value threshold to achieve q < 0.05 FDR-corrected for each predictor');
fprintf('min_d_FDR05: %s', OUT.cohens_d_fdr05_descrip);
disp('dashes')
disp(' ')

% --------------------------------------------
% Print diagnostic info
% --------------------------------------------
printhdr('Input image diagnostic information');

[OUT.group_metrics OUT.individual_metrics values gwcsf gwcsfmean gwcsfl2norm] = qc_metrics_second_level(imgs);

% mahalanobis distance
[ds, expectedds, p, wh_outlier_uncorr, wh_outlier_corr] = mahal(imgs, 'noplot', 'corr');
[dscov, expectedds, p, wh_outlier_uncorr_cov, wh_outlier_corr_cov] = mahal(imgs, 'noplot');

% save quality dat

ind_quality_dat = [mean(OUT.weights)' OUT.individual_metrics.csf_to_gm_signal_ratio' OUT.individual_metrics.gm_L1norm OUT.individual_metrics.csf_L1norm];

%  OUT.individual_metrics.
% corr()

%%
% --------------------------------------------
% Print montages and tables of regions at q < 0.05 FDR
% --------------------------------------------
printhdr('Tables of regions at q < 0.05 FDR');

for i = 1:k
    
    printhdr(sprintf('Predictor %d: %s', i, names{i}));
    
    o2 = canlab_results_fmridisplay(get_wh_image(OUT.t_obj, i), 'montagetype', 'full', 'noverbose');
    
    r = region(get_wh_image(OUT.t_obj, i), 'noverbose');
    
    [posr, negr, OUT.contrast_tables_FDR05{i}] = table(r, 'nolegend'); 
    
    disp(' ')
end

%% 
% --------------------------------------------
% Mask Figure
% --------------------------------------------

create_figure('mask'); axis off; montage(OUT.mask, 'color', [0 .5 0], 'trans');

drawnow, snapnow;

% --------------------------------------------
% Data Figure
% --------------------------------------------
create_figure('data'); axis off; plot(imgs);

drawnow, snapnow;

% --------------------------------------------
% Weights and diagnostics figure
% --------------------------------------------

create_figure('weights and metrics', 2, 2);
xlabel('Image'); ylabel('Weights');
errorbar(mean(OUT.weights), std(OUT.weights), 'bo', 'MarkerFaceColor', [0 0 .5])
title('Mean weights across parcels (s.d. error bars) per image');

subplot(2, 2, 2);
imagesc(OUT.weights);
xlabel('Image'); ylabel('Parcel');
title('Weights by parcel');
colorbar;
axis tight; set(gca, 'YDir', 'Reverse');
        

subplot(2, 2, 3);
xlabel('Image'); ylabel('Z(Weights)');
errorbar(zscore(mean(OUT.weights)), ste(OUT.weights), 'bo-', 'MarkerFaceColor', [0 0 .5], 'LineWidth', 2)
title('Mean weights (s.e. error bars) and quality metrics');
plot(zscore(OUT.individual_metrics.gm_L1norm), 'LineWidth', 2);
plot(zscore(OUT.individual_metrics.csf_L1norm), 'LineWidth', 2);
plot(zscore(ds), 'LineWidth', 2);
plot(zscore(dscov), 'LineWidth', 2);
legend({'Z(Weights)' 'Z(GM L1 norm)' 'Z(CSF L1 norm)' 'Mahal corr dist' 'Mahal cov dist'});

% ***mark off who are outliers

subplot(2, 2, 4)
plot_correlation_matrix(datmatrix, 'dofigure', false);
title('inter-parcel correlations across images');
drawnow, snapnow;


% distribution of weights

end % function

