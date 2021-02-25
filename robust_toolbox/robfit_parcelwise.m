imgs = load_image_set('emotionreg');

X = importdata(which('Wager_2008_emotionreg_behavioral_data.txt'));
Xnames = X.textdata(2);
X = X.data(:, 2);

%%


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
if isempty(X)
    X = ones(t, 1);  % for robust reg
else
    % Make sure intercept is last column
    X = intercept(X, 'end');
end

k = size(X, 2); % number of maps to estimate - one per regressor

% --------------------------------------------
% Initialize output arrays
% --------------------------------------------
% Parcels are rows, images are columns

[nsubjects, maskvol] = deal(zeros(v, 1) .* NaN);     % save numbers of missing values

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
        
        nsubjects(i, 1) = sum(~wasnan);
        maskvol(i, 1) = isok;
        
        weights(i, :) = w;
    
end % loop through nodes

OUT = struct('betas', betas, 'tscores', tscores, 'pvalues', pvalues, 'nsubjects', nsubjects, 'maskvol', maskvol, 'weights', weights);


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

OUT.pthr_FDRq05 = pthr;
OUT.pthr_FDRq05_descrip = 'Highest P-value significant at FDR q < 0.05; threshold at p <= pthr';
OUT.sig_q05 = sig_q05;

% sum(sig_q05)  Sig parcels


% --------------------------------------------
% Transform back to voxelwise output maps
% --------------------------------------------

