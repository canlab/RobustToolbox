function R = robust_reg_nonparam(X,niter,varargin)
%
% Robust regression with nonparametric correction for multiple comparisons
% Correction is returned both mapwise (all voxels) and contiguous
% cluster-by-cluster, for small volume correction (svc)
%
% Overview: takes a data matrix Y (voxels x obs.) and
% model matrix X (subj. x predictors) and returns the robust reg
% coefficients and other things, and significance
%
% X should already contain an intercept.
% Empty X will perform a one-sample t-test.
%
% Usage:
% Results = robust_reg_nonparam(X,niter,'mask',maskname,'data',image_names,'names',x_names,'file',mysavefile)
% OR
% Results = robust_reg_nonparam(Results,niter)
% R =
% robust_reg_nonparam(R,3500,'mask','pag_roi.img','file','nparm_grp3_olstest', ...
% 'names',{'Placebo' 'Order' 'Intercept'},'startover');

%
% F-values for full model (including intercept!) can be requested
% (slower...)
% image_names = EXPT.SNPM.P{3}
% X = EXPT.cov;
% X(:,end+1) = 1;
%
% tor wager, july 06
% Updated: Jan 2008, by Tor Wager, to add primary uncorrected thresholds
%
% Examples:
% Take an old Results (R) struct, add a new mask, keep same input images,
% new output file
% ----------------------------------------------------------------------
% R = robust_reg_nonparam(R,niter,'savefile','nparm_results_pag')

% use typical robust regression right now; slower, but easier


% banner
fprintf(1,'\n* ======================================================================== *')
fprintf(1,'\n*                                                                          *')
fprintf(1,'\n*                 Running function: robust_reg_nonparam.m                  *')
fprintf(1,'\n* Robust regression with nonparametric correction for multiple comparisons *')
fprintf(1,'\n*                                                                          *')
fprintf(1,'\n*                        Tor Wager, version July 2006                      *')
fprintf(1,'\n*                                                                          *')
fprintf(1,'\n* ======================================================================== *\n')
fprintf(1,'Setup:\n')

global mytail
mytail = 'twotail'; % or 'twotail'
savefile = 'nparm_results';
image_names = [];
maskname = [];
names = [];
Y = [];

if isstruct(X)
    % we already have R
    R = X;

    % remove critical field if any input is 'startover'
    R = check_startover(R,varargin);

    [X,R.maskname,R.image_names,R.names,mytail] = setup_variables(R); % somewhat redundant, but does some display and err checking
    if ~isfield(R,'savefile'), R.savefile = savefile; end
    R.mytail = mytail;

else
    % create structure
    R = struct('image_names',image_names,'Y',Y,'maskname',maskname,'X',X,'mytail',mytail,'savefile',savefile);

    % any cell array will cause R to be dealt into multiple structures, so add
    % names here
    R.names = names;
end

% inputs
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'mask', R.maskname = varargin{i+1};
            case 'data'
                if ischar(varargin{i+1}), R.image_names = varargin{i+1};
                else Y = varargin{i+1};
                end
            case 'tail', mytail = varargin{i+1}; R.mytail = mytail;
            case 'file', R.savefile = varargin{i+1};
            case 'names', R.names = varargin{i+1};
                if ~iscell(R.names), error('Enter names in { } cell array.'); end
                %case 'volinfo', volInfo = varargin{i+1};
        end
    end
end

[gook, X, R.maskname, R.names] = check_required_inputs(X,R.maskname,R.image_names,Y, R.names);
if ~gook, return, end

fprintf(1,'\n')
fprintf(1,'Mask of in-analysis voxels is: \n%s\n', R.maskname)
fprintf(1,'First data image name is: \n%s\n', R.image_names(1,:))
fprintf(1,'\n')
fprintf(1,'Tails for max T: %s\n', mytail)
fprintf(1,'\n')
disp(['Saving results periodically in ' R.savefile '.mat'])
disp(['In current working dir: ' pwd])
fprintf(1,'\n')


% -------------------------------------------------------------------
% Get data
% -------------------------------------------------------------------
[R,Y] = load_data_subfunction(R);


% -------------------------------------------------------------------
% Setup permutations and other variables
% -------------------------------------------------------------------

display_str('Setting variables and permutations.');

[n2,v] = size(Y);
if isempty(X), X = ones(n2,1); end
[n,k] = size(X);

% we need to treat the intercept specially later
wh_intercept = find_intercept(X);
R.wh_intercept = wh_intercept;

if n ~= n2, error('data and model sizes do not match.'); end
if no_variance(Y)
    % Check which voxels:
    volInfo = iimg_read_img(R.maskname, 2);
    wh = ~(any(diff(Y)));
    iimg_reconstruct_vols(double(wh)', volInfo, 'outname', 'no_variance_vox.img');
    cl = mask2clusters('no_variance_vox.img');
    cluster_orthviews();
    cluster_orthviews(cl, {[0 0 0]}, 'add', 'solid');
    
    disp('Your mask includes some voxels for which all values are the same.')
    disp('You need to create a mask first that has only valid voxels in it.')
    disp('e.g., use: mask_create_from_image_set(imgs, ''nonparam_mask.img'', size(imgs, 1));');
    disp(' ...where imgs are your input images.');
    error('Some Y vectors have no variability.  You must remove these before running.'); 
end

% names, if empty
if isempty(R.names), for i = 1:k, R.names{i} = ['Reg' num2str(i)]; end, end

% clusters, for by-cluster small volume correction (svc)
[clindx,nvox] = iimg_cluster_index(ones(v,1),R.volInfo.xyzlist');
nvox = nvox{1}';
c = length(nvox);   % number of clusters
R.volInfo.clindx = clindx;
%R.volInfo.nvox = nvox;  % bad; volInfo.nvox is the total num vox.  this is
%in-mask
R.volInfo.c = c;

% for covariates: niter permutations in columns, each to be applied to Y
% could apply to X, but just easier this way
perms = setup_perms(n,niter);

% for intercept, we must permute signs
signmtx = sign(randn(n,niter));

% -------------------------------------------------------------------
% Setup null-hypothesis stat outputs; append to existing if we have them
% -------------------------------------------------------------------
% R contains [maxt,maxf,maxt_by_cl,maxf_by_cl] with empty iterations
% removed and  niter new ones added

[R,permindx] = setup_iterations(R,k,niter,c);

fprintf(1,'Setup done.\n-------------------------------------------------------------------\n');

% End Setup
% -------------------------------------------------------------------




% -------------------------------------------------------------------
% Correct permutation
% -------------------------------------------------------------------
display_str('Correct permutation:');

if isfield(R,'correct')
    fprintf(1,' Already exists in input structure. Using existing.\n');
else
    display_str('Computing.');

    [R.correct.b,R.correct.t,R.correct.p, ...
        R.correct.sig,R.correct.f,R.correct.fp,R.correct.fsig,stat] = robust_reg_matrix(X,Y);

    R.correct.resid = stat.resid;
    R.correct.wts = stat.wts;
    %R.correct.se = stat.se;
    %R.correct.s = stat.s;
end

fprintf(1,'\nSaving results structure R in %s.mat\n',R.savefile);
eval(['save ' R.savefile ' R']);

% print table here
print_maxt(R.correct.t,R.names,wh_intercept)



% -------------------------------------------------------------------
% Update correct stats and Loop through iterations
% -------------------------------------------------------------------
disp('Permuting data to test image-wise and cluster-wise maxima');
fprintf(1,'\n-------------------------------------------------------------------\n');

[R.correct.t,R.correct.p,R.maxt,R.maxt_by_cl, R.primary_uncor] = ...
    robust_pooled_weight_core(X,Y,R.correct.b,R.correct.resid, ...
    R.correct.wts,niter,R.wh_intercept,R.volInfo.clindx);


eval(['save ' R.savefile ' R']);

%R = robust_old_iteration_core(R,niter,perms,X,Y,permindx,clindx,c,signmtx,wh_intercept);


% -------------------------------------------------------------------
% Map-wise and cluster-wise (svc) stats: Print tables and save
% mapwise and svc fields
% -------------------------------------------------------------------

R = robust_nonparam_results(R);


return












% -------------------------------------------------------------------
%
%
% Sub-functions
%
%
% -------------------------------------------------------------------

% -------------------------------------------------------------------
% SETUP: start over by removing critical fields, if we have them
% -------------------------------------------------------------------

function R = check_startover(R,varinputs)
for i = 1:length(varinputs)
    arg = varinputs{i};
    if ischar(arg)
        switch lower(arg)
            case 'startover'

                disp('* Note: Removing all results,data, and permutations for clean start.')

                if isfield(R,'Y'),R = rmfield(R,'Y');end
                if isfield(R,'correct'), R = rmfield(R,'correct'); end
                if isfield(R,'maxt_by_cl'),R = rmfield(R,'maxt_by_cl');end
                if isfield(R,'maxf_by_cl'),R = rmfield(R,'maxf_by_cl');end
                if isfield(R,'maxf'),R = rmfield(R,'maxf');end
                if isfield(R,'maxt'),R = rmfield(R,'maxt');end
                if isfield(R,'volInfo'),R = rmfield(R,'volInfo');end
                if isfield(R,'savefile'),R = rmfield(R,'savefile');end
                if isfield(R,'mapwise'),R = rmfield(R,'mapwise');end
                if isfield(R,'svc'),R = rmfield(R,'svc');end
        end
    end
end
return

% -------------------------------------------------------------------
% SETUP: check inputs
% -------------------------------------------------------------------
function [gook, X, maskname, names] = check_required_inputs(X,maskname,image_names,Y,names)

gook = 1;

% check required
if isempty(X), disp('Must enter valid design matrix in X.'), gook = 0; return, end

wh = intercept(X, 'which');
if isempty(wh)
    X = intercept(X, 'add');
    disp('Check_reqired_inputs: No intercept found; adding one.')
    names = [names {'Intercept'}];
end

if length(names) ~= size(X, 2)
    disp('Check_reqired_inputs: names must be same length as columns of X, including intercept');
end

% mask
if isempty(maskname) && ~isempty(image_names)
    disp('No mask entered.  Using first image as mask.')
    maskname = image_names(1,:);
elseif isempty(maskname)
    disp('Neither mask nor image names entered.  Cannot establish analysis mask.')
    gook = 0;
    return
end

if isempty(Y) && isempty(image_names)
    disp('Neither mask nor image names entered.  Cannot establish analysis mask.')
    gook = 0;
    return
end

return


% -------------------------------------------------------------------
% SETUP: variables from structure
% -------------------------------------------------------------------


function [X,maskname,image_names,names,mytail] = setup_variables(R)
% We have an input R struct
% deal variables out so we can add missing ones

image_names = []; maskname = []; X = []; names = []; mytail = 'twotail';

display_str('Loading variables from input Results structure.'); fprintf(1,'\n');
n = fieldnames(R);
for i = 1:length(n)
    eval([n{i} ' = R.' n{i} ';']);
end

if isempty(X) || isempty(maskname) || isempty(image_names)
    error('Results structure is missing a required field, and is not valid.');
end
return



% -------------------------------------------------------------------
% SETUP: iterations
% -------------------------------------------------------------------

function [R,permindx] = setup_iterations(R,k,niter,c)

if isfield(R,'maxt') && isfield(R,'maxt_by_cl')
    % find empty
    whomit = find(any(R.maxt == 0));
    R.maxt(:,whomit) = [];
    for i = 1:c
        whomit = find(any(R.maxt_by_cl{c} == 0));
        try R.maxt_by_cl{c}(:,whomit) = []; catch disp('Warning: maxt_by_cl is wrong length.'), end
    end
else
    R.maxt = [];
    R.maxt_by_cl = cell(1,c);
end

if isfield(R,'maxf') && isfield(R,'maxf_by_cl')
    % find empty
    whomit = find(R.maxf == 0);
    R.maxf(whomit) = [];
    for i = 1:c
        whomit = find(R.maxf_by_cl{c} == 0);
        try R.maxf_by_cl{c}(whomit) = []; catch disp('Warning: maxf_by_cl is wrong length.'),end
    end
else
    R.maxf = [];
    R.maxf_by_cl = cell(1,c);
end

startat = size(R.maxt,2) + 1;

newt = zeros(k,niter);
newf = zeros(1,niter);
R.maxt = [R.maxt newt];     % max abs t (2-tailed);
R.maxf = [R.maxf newf];

% by-cluster
for i = 1:c
    R.maxt_by_cl{c} = [R.maxt_by_cl{c} newt];
    R.maxf_by_cl{c} = [R.maxf_by_cl{c} newf];
end

endat = size(R.maxt,2);
permindx = startat:endat;   % which indices in output arrays for each permutation

fprintf(1,'\nExisting iterations: %3.0f.  Adding new: %3.0f\n',startat-1,niter);

return


% -------------------------------------------------------------------
% DATA: check for R.Y and load if necessary.  Check image names
% -------------------------------------------------------------------

function [R,Y] = load_data_subfunction(R)
sprintf('Loading data. ');

if isfield(R,'Y') && ~isempty(R.Y)
    disp('R.Y data field found.  Using it.');
    Y = R.Y;

    if ~exist(deblank(R.maskname),'file')
        fprintf(1,'\n* Error: R.maskname is not a valid file.\nImage: %s\nPlease select.\n',R.maskname)
        R.maskname = spm_get(1,'*','Select new mask image');
    end
    
else
    % check image names and get data

    if ~exist(deblank(R.image_names(1,:)),'file')
        fprintf(1,'\n* Error: R.image_names are not valid files.  Assuming names are correct but path is wrong!\n')
        R.image_names = filename_get_new_root_dir(R.image_names, ...
            spm_get(-1,'*','Select root dir for files'),2);
    end

    if ~exist(deblank(R.maskname),'file')
        fprintf(1,'\n* Error: R.maskname is not a valid file.\nImage: %s\nPlease select.\n',R.maskname)
        R.maskname = spm_get(1,'*','Select new mask image');
    end

    % load data and save in R
    Y = iimg_get_data(R.maskname,R.image_names);
    R.Y = Y;
end

if ~isfield(R,'volInfo') || ~strcmp(R.volInfo.fname,R.maskname)
    R.volInfo = iimg_read_img(R.maskname,1);
else
    fprintf(1,'\n* Using existing volInfo structure. Remove before running if mask info has changed.');
end

if size(R.Y,2) ~= R.volInfo.n_inmask
    fprintf(1,'\n* Existing data Y has wrong number of voxels.  Re-loading data from images.');
    Y = iimg_get_data(R.maskname,R.image_names);
    R.Y = Y;
end
    
return


% % % % -------------------------------------------------------------------
% % % % COMPUTE: old iteration core, runs robfit twice on the set
% % % % F-test is still valid, but the t-values are not
% % % % replaced by robust_pooled_weight_core
% % % % -------------------------------------------------------------------
% % % 
% % % function R = robust_old_iteration_core(R,niter,perms,X,Y,permindx,clindx,c,signmtx,wh_intercept)
% % % 
% % % global mytail
% % % 
% % % str = sprintf('Iteration %05d :',0); fprintf(1,str);
% % % 
% % % % timing info
% % % tic
% % % estr = 'Iteration '; % should be same length as iteration text
% % % 
% % % for i = 1:niter
% % % 
% % %     Yperm = Y(perms(:,i),:);
% % % 
% % % 
% % %     % first permute data, for covariates
% % %     % retain conditionality on the marginal distributions of the data &
% % %     % model
% % %     % ----------------------------------------
% % %     if i < 10, str2 = display_str('covariate t-values.');  end
% % %     [bprm,tprm,pprm] = robust_reg_matrix(X,Yperm,0);
% % % 
% % %     if strcmp(mytail,'twotail')
% % %         maxt_i = max(abs(tprm'))';   % for each column    %slower: max([max(tprm'); max(-tprm')]);
% % %     elseif strcmp(mytail,'onetail')
% % %         maxt_i = max(tprm')';
% % %     end
% % % 
% % %     R.maxt(:,permindx(i)) = maxt_i;
% % %     R.maxt_by_cl = update_maxstat_by_cluster(clindx,c,permindx(i),tprm,R.maxt_by_cl);
% % % 
% % % 
% % %     if i < 10, erase_string(str2); end
% % % 
% % %     % then permute signs, for intercept and F
% % %     % ----------------------------------------
% % %     if i < 10, str2 = display_str('intercept and f-values.'); end
% % %     sgns = repmat(signmtx(:,i),1,v);
% % %     [bprm,tprm,pprm,sigprm,fprm] = robust_reg_matrix(X,Yperm .* sgns,0);
% % % 
% % %     if strcmp(mytail,'twotail')
% % %         maxt_i = max(abs(tprm(wh_intercept,:)));
% % %     elseif strcmp(mytail,'onetail')
% % %         maxt_i = max(tprm(wh_intercept,:));
% % %     end
% % % 
% % %     R.maxt(wh_intercept,permindx(i)) = maxt_i;
% % %     R.maxt_by_cl = update_maxstat_by_cluster(clindx,c,permindx(i),tprm,R.maxt_by_cl,wh_intercept);
% % % 
% % %     R.maxf(permindx(i)) = max(fprm);
% % %     R.maxf_by_cl = update_maxstat_by_cluster(clindx,c,permindx(i),fprm,R.maxf_by_cl);
% % % 
% % %     if i < 10, erase_string(str2); end          % erase intercept/fvals
% % %     fprintf(1,'\b\b\b\b\b\b\b');  % erase iteration number
% % % 
% % % 
% % %     % save periodically in case of crash
% % %     if mod(i,10) == 0
% % %         erase_string(estr);
% % %         remhrs = (niter - i) * (toc ./ 10 ./ 3600);
% % %         estr = sprintf('Avg time per iteration: %3.1f s, Est. remaining: %3.1f hrs ',toc ./ 10,remhrs); fprintf(1,estr);
% % % 
% % %         drawnow   % try to fix problem with display in nojvm mode
% % %         pause(.1)
% % % 
% % %         eval(['save ' R.savefile ' R']);
% % %         tic
% % %     end
% % % 
% % %     % print iteration number
% % %     fprintf(1,'%05d :',i);
% % %     drawnow
% % % end

toc;    % just stop the timer


% -------------------------------------------------------------------
% COMPUTE: misc functions
% -------------------------------------------------------------------

% -------------------------------------------------------------------

function val = find_intercept(X)

val = find(~any(diff(X)));

if isempty(val), error('X should contain an intercept.'); end
return

% -------------------------------------------------------------------

function val = no_variance(Y)

val = ~all(any(diff(Y)));

return

% -------------------------------------------------------------------

function erase_string(str1)
fprintf(1,repmat('\b',1,length(str1))); % erase string
return

% -------------------------------------------------------------------

function str = display_str(str)
fprintf(1,str);
return

% -------------------------------------------------------------------

function print_maxt(t,names,wh_intercept)

k = size(t,1);
tmax = max(t');  tmin = min(t');
fprintf(1,'\tMax t\tMin t\t\n')
for i = 1:k, fprintf(1,'%s\t%3.2f\t%3.2f\t\n',names{i},tmax(i),tmin(i)); end
fprintf(1,'\nIntercept is column %3.0f, called %s\n\n',wh_intercept,names{wh_intercept})

return

% -------------------------------------------------------------------

function perms = setup_perms(n,niter)
perms = zeros(n,niter);
for i = 1:niter
    perms(:,i) = randperm(n)';
end
return

% -------------------------------------------------------------------

function maxstat = update_maxstat_by_cluster(clindx,c,i,nullstat,maxstat,pprm,varargin)
% nullstat is k x v, predictors x voxels, matrix of test statistics

global mytail
doint = 0;
if length(varargin) > 1, wh_intercept = varargin{1}; doint = 1; end

for cl = 1:c

    wh = find(clindx == cl);

    if ~doint
        if strcmp(mytail,'twotail')
            mymax = max(abs(nullstat(:,wh)'))';     % this cluster only
        elseif strcmp(mytail,'onetail')
            mymax = max(nullstat(:,wh)')';     % this cluster only
        end
        maxstat{cl}(:,i) = mymax;

        % number of voxels above threshold, p<.05 2-tailed
        % only for t contrast, not F (use empty input for F)
        %if ~isempty(pprm)
        %    sum(pprm(:,wh)' < .05);
        %end

    else
        if strcmp(mytail,'twotail')
            mymax = max(abs(nullstat(wh_intercept,wh))); % this cluster only
        elseif strcmp(mytail,'onetail')
            mymax = max(nullstat(wh_intercept,wh)); % this cluster only
        end
        maxstat{cl}(wh_intercept,i) = mymax;
    end

end

return
