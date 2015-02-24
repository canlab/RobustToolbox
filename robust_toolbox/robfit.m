% EXPT = robfit(EXPT, [which cons], [do ols], [mask], [nocenter])
%
% This function does a robust fit on contrasts specified in EXPT.SNPM.connames,
% using images in EXPT.SNPM.P
% The model should be stored (without intercept) in EXPT.cov
% Empty EXPT.cov will result in a one-sample t-test
%
% For each contrast, it saves con, t, p, filtered t imgs, and clusters.mat
% p-value images are 2-tailed!
%
% All predictors should be in a matrix of column vectors, stored in EXPT.cov.
% 
% All covariates are scaled to mean=0 and sd=1, UNLESS they are composed of 0's, 1's, and -1's,
% in which case they are assumed to be contrast codes and are not scaled.
% the 'nocenter' flag prevents automatic scaling.  This behavior was
% slightly modified on 5/29/2014 -- see programmer's notes below. 
%
% Intercept is added as 1st predictor
%
% fields should be:
% EXPT.SNPM.P           containing image names of individual subjects
% EXPT.SNPM.connames    str mtx of contrast names for each P string matrix
% EXPT.SNPM.connums     contrast numbers for each P
% EXPT.cov              empty for 1-sample ttest, or containing covariates
% [mask] is the name of the mask file to use
% [nocenter] is a flag to be used if you don't want to center some
% variables (e.g., you have three different groups, two sets of contrasts)
%
%
% [Which cons] optional vector of which arrays in EXPT.SNPM.P to use
% e.g., [1 2 4] runs the first, second, and 4th sets of image names
%
% Example:
% %load('subjs');
% %cd into results directory first!
% EXPT.cov = [];
% load(sprintf('%s/SPM.mat', subjs{1}));
% wh_spm_cons = strmatch('T', char(SPM.xCon.STAT));
% EXPT.SNPM.connames = strvcat(SPM.xCon(wh_spm_cons).name);
% EXPT.SNPM.connums = 1:size(EXPT.SNPM.connames, 1);
% EXPT.SNPM.P = {};
% for i=1:length(wh_spm_cons)
%       EXPT.SNPM.P{i} = filenames(sprintf('%s/*/con_%04d.img', pwd(), wh_spm_cons(i)), 'char');
% end
%
% EXPT.mask = fullfile(study_root, 'rscalped_avg152T1_graymatter.img');
% %or
% EXPT.mask = fullfile(study_root, 'scalped_avg152T1_graymatter.img');
% %or if you need to make a new mask:
% mask = which('scalped_avg152T1_graymatter_smoothed.img');
% [dummy, EXPT.mask] = reslice_imgs(EXPT.SNPM.P{1}(1,:),mask);
%
% EXPT = robfit(EXPT, 1:length(EXPT.SNPM.connums), 0, EXPT.mask);
%
% To do OLS, cons 1-4, and use gray matter mask:
% EXPT = robfit(EXPT, 1:4, 1, EXPT.mask);
%
% Programmer's notes:
%
% On 5/19/2014, covariates are automatically scaled unless composed solely
% of 0s, 1s, and -1s.  Previously, if ANY covariate was not composed soley
% of 0s, 1s, and -1s, then NO covariate was centered.  Now, this decision
% is made for each covariate separately.

function EXPT = robfit(EXPT,varargin)

linestr = '____________________________________________________________';

% -----------------------------------------------------
% center and scale predictors
% -----------------------------------------------------
if ~isfield(EXPT,'cov'), EXPT.cov = []; end
covt = EXPT.cov;

if length(varargin)>3, nocenter = varargin{4}; else nocenter = 0; end

if ~isempty(covt)
    
    fprintf(1,'Robfit.m\n%s\nPreparing covariates:\n',linestr);
    
    if length(covt) > size(covt,1), covt = covt'; fprintf(1,'Transposed covariate matrix.\n'); end
    
    if nocenter == 0
        for i = 1:size(covt,2)
            if all( covt(:,i) == 1 | covt(:,i) == -1 | covt(:,i) == 0)
                fprintf(1,'Cov %3.0f: Contrast codes; no centering. Will be saved in rob_tmap_%04d\n',i,i+1);
                % Contrast codes; do not center; this gives Type III SS
            else
                fprintf(1,'Cov %3.0f: Centering and scaling to unit variance. Will be saved in rob_tmap_%04d\n',i,i+1);
                covt(:,i) = scale(covt(:,i));
            end
            
        end
    else
        for i = 1:size(covt,2)
            fprintf(1,'Cov %3.0f: Requested no centering. This will change the interpretation of the intercept!\nWill be saved in rob_tmap_%04d\n',i,i+1);
        end
    end
    
else
    fprintf(1,'Robfit.m\n%s\nIntercept only (no covariates)\n',linestr);
end

EXPT.cov = covt;

fprintf(1,'%s\n',linestr);

% intercept: automatically added if doing robustfit.m
% covt = [ones(1,size(covt,1)) covt];
%covt(:,end+1) = 1;



if ~isempty(varargin), wh=varargin{1}; else wh=1:length(EXPT.SNPM.P); end
if length(varargin)>1, dools = varargin{2}; else dools = 1; end
if length(varargin)>2, domask = varargin{3}; else domask = 0; end

% ----------------------------------------------------
% check some things
% ----------------------------------------------------

if ~isfield(EXPT, 'SNPM')
    disp('robfit: Enter EXPT.SNPM substructure with file names, contrast names, contrast numbers!');
    disp('Skipping analysis.')
    return
elseif ~isfield(EXPT.SNPM, 'P')
    disp('robfit: Enter EXPT.SNPM.P cell array with file names, one cell per robust analysis');
    disp('Skipping analysis.')
    return
elseif ~isfield(EXPT.SNPM, 'connames')
    disp('robfit: Enter EXPT.SNPM.connames string matrix with analysis/contrast name for each image set');
    disp('Creating dummy names and proceeding')
    tmp = {};
    for i = 1:length(EXPT.SNPM.P)
        tmp{i} = sprintf('analysis%3.0f', i);
    end
    EXPT.SNPM.connames = char(tmp{:});
elseif size(EXPT.SNPM.connames, 1) ~= length(EXPT.SNPM.P)
    if size(EXPT.SNPM.connames, 2) == length(EXPT.SNPM.P)
        EXPT.SNPM.connames = EXPT.SNPM.connames'; %trivial: sometimes need to invert to get in expected format
    else
        disp('robfit: Contrast names and image file list are different lengths!! Check this!!');
        error('Fix names/image lists and re-run.');
    end
elseif ~isfield(EXPT.SNPM, 'connums')
    disp('robfit: Enter EXPT.SNPM.connums vector with analysis number for each robust analysis');
    disp('Creating default vector and proceeding. Fix to avoid this warning.')
    EXPT.SNPM.connums = 1:length(EXPT.SNPM.P);
end

% ----------------------------------------------------
% * run robust regression for selected sets
% ----------------------------------------------------

basedir = pwd;
fprintf('robfit: Will write to robust???? directories in: %s\n', basedir);

for i = wh
    % * run robust regression
    % ----------------------------------------------------
    warning off
    
    disp(''); disp(['Robfit.m - working on ' EXPT.SNPM.connames(i,:)])
    
    if dools, disp('Running OLS and IRLS comparison (slower) - to turn this off, use 0 as 3rd input argument.'),end
    disp(linestr)
    
    if domask
        fprintf(1,'Using mask image: \n%s\n', domask);
    else
        fprintf(1,'No mask found. Analyzing all voxels.\n');
    end
    
    if dools
        [EXPT.SNPM.rob_betas{i},EXPT.SNPM.ols_betas{i}] = rob_fit(EXPT.SNPM.P{i}, covt, EXPT.SNPM.connums(i), dools, domask, EXPT.SNPM.connames(i,:), basedir);
    else
        [EXPT.SNPM.rob_betas{i}] = rob_fit(EXPT.SNPM.P{i}, covt, EXPT.SNPM.connums(i), dools, domask, EXPT.SNPM.connames(i,:), basedir);
    end
    warning on
    
    % * compare
    % ----------------------------------------------------
    
    
end



end % main function



function [newP,newP2] = rob_fit(P, covt, index, dools, domask, analysisname, basedir)

% --------------------------------------------
% Set up output directory
% --------------------------------------------
if index < 10, myz = '000';  else myz = '00';  end
mydir = ['robust' myz num2str(index)];

fprintf(1,'\nSaving results in : %s ', mydir);

cd(basedir)
cwd = fullfile(basedir, mydir);

if ~exist(cwd,'dir')
    fprintf(1,'(Creating directory)');
    eval(['mkdir ' cwd])
else
    fprintf(1,'(Directory already exists. Will overwrite.)');
end
fprintf(1,'\n');
cd(cwd)

% initialize SETUP
SETUP = struct('dir', mydir, 'name', analysisname, 'files', P, 'covariates', covt, 'dools', dools, 'domask', domask);

% --------------------------------------------
% get data and some var sizes
% --------------------------------------------
V = spm_vol(P);
v = spm_read_vols(V);
nvols = size(v, 4);

n_out_imgs = size(covt,2) + 1;

% --------------------------------------------
% find the in-analysis voxels
% --------------------------------------------

% mask, if specified
if domask
    % sample mask data in space of input images
    vm = scn_map_image(domask, P(1,:));
    
    % make sure mask is 1's or 0's
    vm = double(vm > 0);
    
    for i = 1:nvols
        v(:, :, :, i) = v(:, :, :, i) .* vm;
    end
    
end

N = size(P, 1);

% minimum allowed observations.  robust requires 3 df
if isempty(covt), k = 1;  %only estimate intercept
else k = size(covt,2); end
Ncrit = N-k-3;
fprintf('Minimum allowed observations per voxel:%d\n', Ncrit);

wh = N - sum(isnan(v) | v == 0, 4); %how many observations at this voxel
wh = wh >= Ncrit; % voxels to keep

[x,y,z] = ind2sub(size(wh),find(wh));
tmp = sum(wh(:)); if tmp == 0, disp('No voxels in analysis!'), end
whplanes = unique(z);

% --------------------------------------------
% set up output arrays
% --------------------------------------------

% F-maps for omnibus Act + Covs
%     fmap = zeros(size(wh)) .* NaN;     % save F vals
%     fpmap = zeros(size(wh)) .* NaN;     % save p vals

[nsubjects, maskvol] = deal(zeros(size(wh)) .* NaN);     % save numbers of missing values

weights = zeros([size(wh) N]) .* NaN; % robust reg weights

[vo, vot, vop, vo2, vot2, vop2, voz3, vop3] = deal(cell(1, n_out_imgs));

for j = 1 : n_out_imgs
    
    %IRLS
    vo{j} = zeros(size(wh)) .* NaN;     % save betas
    vot{j} = zeros(size(wh)) .* NaN;    % save t values
    vop{j} = zeros(size(wh)) .* NaN;    % save p values
    
    % OLS
    vo2{j} = zeros(size(wh)) .* NaN;     % save betas
    vot2{j} = zeros(size(wh)) .* NaN;    % save t values
    vop2{j} = zeros(size(wh)) .* NaN;    % save p values
    
    % comparison between IRLS and OLS
    voz3{j} = zeros(size(wh)) .* NaN;    % save Z scores for difference
    vop3{j} = zeros(size(wh)) .* NaN;    % save p values
    
end
fprintf(1,'\nImage name for first image:\n %s\nNumber of images: %3.0f\n\n%6.0f voxels, %3.0f planes in analysis\n\t',P(1,:),size(P,1),sum(wh(:)),length(whplanes))
fprintf(1,'Done: 000%%');
savewh = sum(wh(:));

% --------------------------------------------
% set up design
% --------------------------------------------
if isempty(covt)
    X = ones(nvols,1);  % for robust reg
else
    X = [ones(size(covt,1),1) covt]; % for robust reg, run with no intercept!
end

k = size(X, 2);

% --------------------------------------------
% perform regression
% --------------------------------------------

for i = 1:length(x)
    
    t = squeeze(v(x(i),y(i),z(i),:));
    
    %   b = pinv(covt) * t;         % regular linear model fit
    
    % alternative way: uses Rousseeuw's algorithm - much slower
    % The problem is that this is biased towards finding results!  Improper
    % control of false positive rate.
    %Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
    %Journal of the American Statistical Association, Vol. 79, pp. 871-881.
    
    %[res]=fastmcd_noplot([covt t]);
    % remove n most extreme outliers and recompute correlation
    %wh = res.flag==0; tnew = t; xnew = covt;
    %tnew(wh) = []; xnew(wh,:) = [];
    %[bb,dev,stats] = glmfit(xnew,tnew);
    
    % * fit models for IRLS and OLS
    % ----------------------------------------------------
    
    % robustfit, IRLS,
    doirls = 1; %dools = 1;
    if doirls
        
        % Remove NaNs and run only if we have enough observations
        [wasnan, Xvox, tvox] = nanremove(X, t);
        
        isok = sum(~wasnan) > 4; % THIS CHECK LIKELY UNNECESSARY B/C EXCLUDING BAD VOXELS ABOVE - YONI
        
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
    end % if robustfit
    
    for j = 1:k
        vo{j}(x(i),y(i),z(i)) = bb(j);
        vot{j}(x(i),y(i),z(i)) = stats.t(j);
        vop{j}(x(i),y(i),z(i)) = stats.p(j);
        
        nsubjects(x(i),y(i),z(i)) = sum(~wasnan);
        maskvol(x(i),y(i),z(i)) = isok;
        
        weights(x(i),y(i),z(i), :) = w;
    end
    
    % save F-stats
    % Tor removed Aug 2012 - not very useful overall, so save the time.
    %         [Fobs, p, dfb, dfe] = F_test_no_intercept(X,t,stats.s);
    %         fmap(x(i),y(i),z(i)) = Fobs;
    %         fpmap(x(i),y(i),z(i)) = p;
    
    
    if dools
        if sum(~wasnan) > 4
            % OK to run
            [bb2,dev,stats2]=glmfit(Xvox,tvox,[],[],'off',[],[],'off');
            
            % * calculate Z-test for comparison
            % ----------------------------------------------------
            if i == 1, df = stats2.dfe;  [mn,vv] = tstat(df); end
            
            zdiff = stats.t ./ sqrt(vv) - stats2.t ./ sqrt(vv); % z-scores by dividing by t-distribution variance
            zdiff = zdiff ./ sqrt(2);                           % diff between z-scores is distributed with var=sum of var(z1) + var(z2)
            % ...thus,sigma
            % (z1 - z2) = sqrt(2)
            % e.g., - http://www.mathpages.com/home/kmath046.htm
            pdiff = 2 * (1 - normcdf(abs(zdiff)));              % 2-tailed
            
            
        else
            % empty voxel: not enough data
            bb2 = NaN .* ones(k, 1);
            stats2 = struct('t', bb2, 'p', ones(k, 1), 'dfe', NaN);
            zdiff = NaN;
            pdiff = 1;
            
        end
        
        % save in matrix for output
        for j = 1:length(bb)
            vo2{j}(x(i),y(i),z(i)) = bb2(j);
            vot2{j}(x(i),y(i),z(i)) = stats2.t(j);
            vop2{j}(x(i),y(i),z(i)) = stats2.p(j);
            
            voz3{j}(x(i),y(i),z(i)) = zdiff(j);
            vop3{j}(x(i),y(i),z(i)) = pdiff(j);
        end
        
    end
    
    % display progress
    if rem(i,100) == 0
        fprintf(1,'\b\b\b\b%03d%%',round(100*i / length(x)));
    end
    
end
fprintf(1,' done. ')
fprintf(1,'\twriting volumes.\n')


%cd(mydir)

% save setup stuff
try
    %         SETUP.files = P;
    %         SETUP.covariates = covt;
    SETUP.X = X;
    SETUP.V = V;
    SETUP.planes = whplanes;
    SETUP.nvoxels = savewh;
    SETUP.sample_robust_res = stats;
    SETUP.df = stats.dfe;
%     SETUP.descr = 'F-test df follow.';
%     SETUP.dfb = dfb;
%     SETUP.dfe = dfe;
    save SETUP SETUP
catch
    warning('Error creating SETUP file');
end

V = V(1); [d,f,e] = fileparts(V.fname);

% set data type to float
switch(spm('Ver'))
    case 'SPM2'
        V.dim(4) = spm_type('float');
    case 'SPM5'
        V.dt(1) = spm_type('float32');
    case 'SPM8'
        V.dt(1) = spm_type('float32');
    otherwise
        error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
end

%cwd = pwd;

%     fprintf(1,'Writing F- and p-images for full model. \n')
%     fprintf(1,'Note: Significance of F includes intercept, so overall activation contributes to significance\n');

% * write output images for fmap
% ----------------------------------------------------
%     V.fname = fullfile(cwd,['rob_fmap_full' e]);
%     V.descrip = 'IRLS robust regression F-map for full model, including intercept';
%     disp(['Writing ' V.fname])
%     spm_write_vol(V,fmap);
%
%     V.fname = fullfile(cwd,['rob_pmap_full' e]);
%     V.descrip = 'IRLS robust regression p-map for full model, including intercept';
%     disp(['Writing ' V.fname])
%     spm_write_vol(V,fpmap);

V.fname = fullfile(cwd,['nsubjects' e]);
V.descrip = 'Number of observations (subjects) with valid data';
disp(['Writing ' V.fname])
spm_write_vol(V, nsubjects);

V.fname = fullfile(cwd,['mask' e]);
V.descrip = 'Mask of in-analysis voxels (at least 4 valid observations)';
disp(['Writing ' V.fname])
spm_write_vol(V, maskvol);

% Weights
% V.fname = fullfile(cwd,['weights' e]);
%     V.descrip = 'Robust regression weights for each subject';
% disp(['Writing ' V.fname])
% spm_write_vol(V, squeeze(weights(:, :, :, 1)));

maskvol(isnan(maskvol)) = 0;

try
    
    V2 = V; % copy to avoid problems later with 4-D files where they should be 3-D
    V2.fname = fullfile(cwd,['weights' e]);
    disp(['Writing 4-D weight file: ' V2.fname])
    for i = 1:N
        V2.n(1) = i; % volume number
        V2.descrip = 'Robust regression weights for each subject';
        
        myw = squeeze(weights(:, :, :, i));
        myw(~maskvol) = NaN;
        spm_write_vol(V2, myw);
    end
catch
    disp('Error writing weights.');
end


fprintf(1,'Writing t- and p-images for univariate effects in model. \n')
for i = 1:length(vo)

    % * write output images for IRLS
    % ----------------------------------------------------
    
    if i < 10, myz = '000';  else myz = '00';  end
    V.fname = fullfile(cwd,['rob_beta_' myz num2str(i) e]);
    V.descrip = 'IRLS robust regression betas, beta_0001 is intercept';
    
    if i == 1, newP = V.fname;
    else
        newP = str2mat(newP,V.fname);
    end
    
    disp(['Writing ' V.fname])
    spm_write_vol(V,vo{i});
    
    V.fname = fullfile(cwd,['rob_tmap_' myz num2str(i) e]);
    V.descrip = 'IRLS robust regression t-scores, tmap_0001 is intercept';
    disp(['Writing ' V.fname])
    spm_write_vol(V,vot{i});
    
    V.fname = fullfile(cwd,['rob_p_' myz num2str(i) e]);
    V.descrip = 'IRLS robust regression p values, rob_p_0001 is intercept';
    disp(['Writing ' V.fname])
    spm_write_vol(V,vop{i});
    
    % * write output images for OLS
    % ----------------------------------------------------
    
    if dools
        V.descrip = 'OLS regression betas from robfit.m, beta_0001 is intercept';
        V.fname = fullfile(cwd,['ols_beta_' myz num2str(i) e]);
        
        if i == 1, newP2 = V.fname;
        else
            newP2 = str2mat(newP,V.fname);
        end
        
        disp(['Writing ' V.fname])
        spm_write_vol(V,vo2{i});
        
        V.fname = fullfile(cwd,['ols_tmap_' myz num2str(i) e]);
        V.descrip = 'OLS regression t-scores from robfit.m, tmap_0001 is intercept';
        disp(['Writing ' V.fname])
        spm_write_vol(V,vot2{i});
        
        V.fname = fullfile(cwd,['ols_p_' myz num2str(i) e]);
        V.descrip = 'OLS regression p values from robfit.m, rob_p_0001 is intercept';
        disp(['Writing ' V.fname])
        spm_write_vol(V,vop2{i});
        
        % * write output images for comparison
        % ----------------------------------------------------
        V.descrip = 'IRLS-OLS difference z-scores from robfit.m, beta_0001 is intercept';
        V.fname = fullfile(cwd,['irls-ols_z_' myz num2str(i) e]);
        
        if i == 1, newP2 = V.fname;
        else
            newP2 = str2mat(newP,V.fname);
        end
        
        disp(['Writing ' V.fname])
        spm_write_vol(V,voz3{i});
        
        V.fname = fullfile(cwd,['irls-ols_p_' myz num2str(i) e]);
        V.descrip = 'IRLS-OLS 2-tailed p-values from robfit.m, tmap_0001 is intercept';
        disp(['Writing ' V.fname])
        spm_write_vol(V,vop3{i});
        
    end % do ols
    
end % contrasts (vo)

cd(basedir)

end % subfunction