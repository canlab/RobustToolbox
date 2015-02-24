function EXPT = robseed(EXPT,cl,varargin)
% EXPT = robseed(EXPT,cl,[which cons],[dools],[mask],[doglobal])
%
% Tor Wager, last modified 4/12/04 to add OLS and difference maps
%       Modified may 21 by tor to force output data type to float
%
% This function does a robust GLM 
% that regresses each of the images in EXPT.SNPM.P
% (voxel by voxel) on each of the 
% seed timecourses in cl.indiv_timecourse or cl.timecourse.
% Each seed is run in a separate analysis, and is saved in a separate set
% of results images. 
% Seed data is observations x seeds, and comes from one of the following
% sources (searched for in order):
% 1) EXPT.seeds, obs x seeds, will run analysis for each column sep.
%       In this case, cl input can be empty
%
% 2) cl.indiv_timeseries.  Assumes each of cl(x).indiv_timeseries contains
% a vector of observations for one seed region (cl(x)).  Vectors for all
% seed regions should contain the same number of observations.
%
% 3) cl(x).timeseries.  Same as above.
% A copy of the data used as seeds is stored in SETUP.covariates.
% Seeds are centered and scaled (coverted to z-scores) in this function.
% Intercept is added as 1st predictor
%
% In the output:
% Each output directory robseed00** is a contrast
% Each covariate image (after the intercept) is a new seed
%   e.g., rob_tmap_0001.img is the intercept
%         rob_tmap_0002.img is the first seed
%         rob_tmap_0003.img is the second seed, etc.
%
% For each correlation with a seed, it saves con, t, p, filtered t imgs, and clusters.mat 
% p-value images are 2-tailed.
%
%
% see get_expt_info.m for how to make the EXPT structure
% fields should be:
% EXPT.SNPM.P           containing image names of individual subjects
% EXPT.SNPM.connames    str mtx of contrast names for each P string matrix
% EXPT.SNPM.connums     contrast numbers for each P
% EXPT.cov              empty for 1-sample ttest, or containing covariates
% [mask] is the name of the mask file to use
%
% [Which cons] optional vector of which arrays in EXPT.SNPM.P to use
% e.g., [1 2 4] runs the first, second, and 4th sets of image names
%
% Example:
% EXPT.seeds = rand(20, 3); % 20 observations, 3 seeds, made-up example.
% dools = 0;    % do not write OLS images in addition to robust images
%  %(default is 1, or "Yes")
% mask = 'my_brain_mask.img';  % should be in same space (dims, voxel
% sizes, origin) as images in analysis
% doglobal = 0; %  include global image intensity as covariate.
% EXPT = robseed(EXPT, [], dools, mask, doglobal);   
%
% See robfit.m for additional formatting help and general "how it works."

% -----------------------------------------------------
% center and scale predictors
% -----------------------------------------------------
if isfield(EXPT,'seeds')
    covt = EXPT.seeds;
else
    try 
        covt = cat(2,cl.indiv_timeseries);
        fprintf(1,'robseed.  using individual timeseries from cl.indiv_timeseries.\n');
    catch
        try
            covt = cat(2,cl.timeseries);
            fprintf(1,'robseed.  using average timeseries from cl.timeseries as covariate.\n');
        catch
            error('Enter seed region data EXPT.seeds, or in cl.indiv_timeseries or cl.timeseries.');
        end
    end
end

if ~isempty(covt)
    if length(covt) > size(covt,1), covt = covt';  end
    
    covt = scale(covt);
    
    %covt = covt - repmat(mean(covt),size(covt,1),1);
    %for i = 1:size(covt,2)
    %    covt(:,i) = covt(:,i) ./ var(covt(:,i));
    %end
end

% intercept: automatically added if doing robustfit.m
% covt = [ones(1,size(covt,1)) covt];
%covt(:,end+1) = 1;


if length(varargin)>0,wh=varargin{1}; else wh=1:length(EXPT.SNPM.P); end
if length(varargin)>1,dools = varargin{2};  else  dools = 1;  end
if length(varargin)>2,domask = varargin{3};  else  domask = 0;  end
if length(varargin)>3,doglobal = varargin{4};  else  doglobal = 0;  end

% ----------------------------------------------------
% save seeds, etc., in current dir.
% ----------------------------------------------------
robseed_seeds = covt;
robseed_whole_brain_imgs = EXPT.SNPM.P(wh);
if isfield(EXPT,'seednames'), robseed_seed_names = EXPT.seednames;  end
save robseed_setup_info robseed*
clear robseed*



% ----------------------------------------------------
% * set up the gui figure
% ----------------------------------------------------
f = [];
%f = figure('Color','w');
%tmp = get(gcf,'Position') .* [1 1 .5 .1];
%set(gcf,'Position',tmp)
%set(gcf,'MenuBar','none','NumberTitle','off')
%figure(f), set(gca,'Xlim',[0 100])

% ----------------------------------------------------
% * run robust regression for selected sets
% ----------------------------------------------------

for i = wh
    % * run robust regression 
    % pass in ONLY the covt for the current seed region
    % corresponding to this run/directory.
    % one run/dir per seed region (i.e., cluster)
    % ----------------------------------------------------
    warning off
    %str = ['Cl ' num2str(i)];
    %if isfield(cl,'shorttitle'), str = [str ': ' cl(i).shorttitle]; end
    str = EXPT.SNPM.connames(i,:);
    try disp(['Robseed.m - working on ' str]), catch end
    if dools, disp('Running OLS and IRLS comparison (slower) - to turn this off, use 0 as 3rd input argument.'),end
    disp('____________________________________________________________')
    
    % ----------------------------------------------------
    % get globals, if entered
    % ----------------------------------------------------
    if doglobal
        %glob = spm_global(EXPT.SNPM.P{i});  % crashing!
        disp(['Entering global contrast values as covariate.']);
        glob = tor_global(EXPT.SNPM.P{i},domask);
        
    else
        glob = [];
        disp(['No global contrast covariate.']); 
    end
    
    if dools
        [EXPT.SNPM.rob_betas{i},EXPT.SNPM.ols_betas{i}] = rob_fit(EXPT.SNPM.P{i},covt,EXPT.SNPM.connums(i),f,dools,domask,glob);
    else
        [EXPT.SNPM.rob_betas{i}] = rob_fit(EXPT.SNPM.P{i},covt,EXPT.SNPM.connums(i),f,dools,domask,glob);
    end
    warning on
    
    save seed_clusters cl
end



return








function [newP,newP2] = rob_fit(P,covt,index,f,dools,domask,glob)

    
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

% mask, if spc
if domask
    
    % sample mask data in space of input images
        vm = scn_map_image(domask, P(1,:));

        % make sure mask is 1's or 0's
        vm = double(vm > 0);
        
        for i = 1:nvols
            v(:, :, :, i) = v(:, :, :, i) .* vm;
        end

end

wh=sum(v,4);
wh = ~isnan(wh) & wh ~= 0; 
%wh(all(v == 0,4) = 0;
[x,y,z] = ind2sub(size(wh),find(wh));

tmp = sum(wh(:)); if tmp == 0, disp('No voxels in analysis!'), end
whplanes = unique(z);

% --------------------------------------------
% set up output arrays
% --------------------------------------------

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
fprintf(1,'\nworking on %s: %6.0f voxels, %3.0f planes in analysis\n\t',P(1,:),sum(wh(:)),length(whplanes))
savewh = sum(wh(:));



fprintf('Done: 00%');

% --------------------------------------------
% perform regression
% --------------------------------------------

for i = 1:length(x)
  
    t = squeeze(v(x(i),y(i),z(i),:));
    %   b = pinv(covt) * t;         % regular linear model fit
    
    % NEW way: uses Rousseeuw's algorithm - much slower
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
        if isempty(covt)
            % intercept only (one-sample t-test)
            [bb,stats]=robustfit([ones(length(t),1) glob],t,'bisquare',[],'off');   
        else
            bbsave = []; tsave = []; psave = [];
            for seed = 1:size(covt,2)
                [bb,stats]=robustfit([covt(:,seed) glob],t);  
                if seed == 1
                    bbsave = [bb(1:2)];  
                    tsave = [tsave; stats.t(1:2)]; 
                    psave = [psave; stats.p(1:2)];  
                else  
                    bbsave(end+1) = bb(2); 
                    tsave(end+1) = stats.t(2);
                    psave(end+1) = stats.p(2);
                end
            end
        end
      
    end % if robustfit
    
    for j = 1:length(bbsave)
        vo{j}(x(i),y(i),z(i)) = bbsave(j);
        vot{j}(x(i),y(i),z(i)) = tsave(j);
        vop{j}(x(i),y(i),z(i)) = psave(j);
    end
    
    
    if dools
        if isempty(covt)
            % intercept only (one-sample t-test)
            [bb2,dev,stats2]=glmfit([ones(length(t),1) glob],t,[],[],'off',[],[],'off');
        else
            bb2save = []; t2save = []; p2save = [];
            for seed = 1:size(covt,2)
                [bb2,dev,stats2]=glmfit([covt(:,seed) glob],t);
                if seed == 1, 
                    bb2save = [bb2(1:2)];  
                    t2save = [t2save; stats2.t(1:2)]; 
                    p2save = [p2save; stats2.p(1:2)];  
                else  
                    bb2save(end+1) = bb2(2); 
                    t2save(end+1) = stats2.t(2);
                    p2save(end+1) = stats2.p(2);
                end
            end
        end
        
        for j = 1:length(bb2save)
            vo2{j}(x(i),y(i),z(i)) = bb2save(j);
            vot2{j}(x(i),y(i),z(i)) = t2save(j);
            vop2{j}(x(i),y(i),z(i)) = p2save(j);
        end
        
        % * calculate Z-test for comparison
        % ----------------------------------------------------
        if i == 1, df = stats2.dfe;  [mn,vv] = tstat(df); end
        
        zdiff = tsave ./ sqrt(vv) - t2save ./ sqrt(vv); % z-scores by dividing by t-distribution variance
        zdiff = zdiff ./ sqrt(2);                           % diff between z-scores is distributed with var=sum of var(z1) + var(z2)
                                                            % ...thus,sigma
                                                            % (z1 - z2) = sqrt(2) 
                                                            % e.g., - http://www.mathpages.com/home/kmath046.htm
        pdiff = 2 * (1 - normcdf(abs(zdiff)));              % 2-tailed
         
        for j = 1:length(bb)
            voz3{j}(x(i),y(i),z(i)) = zdiff(j);
            vop3{j}(x(i),y(i),z(i)) = pdiff(j);
        end
    end



    if rem(i,10) == 0
        %figure(f), try,barh(100*i / length(x)),catch,end
        %set(gca,'Xlim',[0 100]),set(gca,'YTickLabel',i),drawnow
        %text(5,.5,['Voxel ' num2str(i) ' of ' num2str(length(x))],'Color','r')
        fprintf(1,'\b\b\b%02d%%',round(100*i / length(x)));
    end
   
end
fprintf(1,' done. ')
fprintf(1,'\twriting volumes.\n')

if index < 10, myz = '000';  else  myz = '00';  end
mydir = ['robseed' myz num2str(index)];
eval(['mkdir ' mydir])
cd(mydir)


% save setup stuff
if ~(exist('df') == 1), df = [];  end
if ~(exist('stats') == 1), stats = [];  end
SETUP.df = df;
SETUP.sample_robust_res = stats;

try
    SETUP.files = P;
    SETUP.covariates = covt;
    SETUP.V = V;
    SETUP.planes = whplanes;
    SETUP.nvoxels = savewh;
    SETUP.glob = glob;
    save SETUP SETUP
catch
    disp('Error creating SETUP file');
end



V = V(1); [d,f,e] = fileparts(V.fname);
cwd = pwd;

% set data type to float
    switch(spm('Ver'))
        case 'SPM2'
            V.dim(4) = spm_type('float');
        case 'SPM5'
            V.dt(1) = spm_type('float32');
        otherwise
            error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
    end


for i = 1:length(vo)
    

    % * write output images for IRLS
    % ----------------------------------------------------
    
    if i < 10, myz = '000';  else  myz = '00';  end
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



cd ..

return