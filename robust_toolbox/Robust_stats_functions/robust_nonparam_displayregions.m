function [A,R,cl_extent,dat_extent,cl_extent_pos,cl_extent_neg] = robust_nonparam_displayregions(R,wh_regressor,wh_field,cortype,varargin)
% [A,R,cl_extent,dat_extent,cl_extent_pos,cl_extent_neg] = robust_nonparam_displayregions(R,wh_regressor,wh_field,cortype,[optional args])
%
% wh_regressor = 1;       % number of regressor to get results for
% wh_field = 't';         % t or f %%% Only works now for t %%%
% cortype = 'svc';        % mapwise or svc
%
% optional:
% 'overlay', followed by overlay
% 'extended', get partial correlations in output tables
%             and save data from significant voxels in A.cl_sig for plots
% 'data', followed by data matrix
%             (data in R.Y, R.image_names, or entered here is necessary for
%             extended output
%
% tor wager, july 2006
%
% TO-DO: fill in other cortype and wh_field options
% step-down test option
%
% Example: Get regressor 3 with extended correlation output
% [A,R] = robust_nonparam_displayregions(R,3,'t','svc','overlay',EXPT.overlay,'extended');

% -------------------------------------------------------------------
% Set up input arguments
% -------------------------------------------------------------------
doextended = 0;
overlay = which('scalped_single_subj_T1.img');
cl_extent = [];

% inputs
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'overlay', overlay = varargin{i+1};
            case 'data', R.Y = varargin{i+1};

            case 'extended', doextended = 1;
        end
    end
end

% -------------------------------------------------------------------
% Sizes and data, print banner
% -------------------------------------------------------------------

% numc = R.volInfo.c;         % number of clusters
% n = R.volInfo.nvox;
[k,v] = size(R.correct.t);

% all voxels in cluster
allvox = ones(v,1); % for actual t-values : R.correct.t(wh_regressor,:)';
A.cl_all = iimg_indx2clusters(allvox,R.volInfo);

% banner -- and check fields
% This runs robust_nonparam_results if the critical fields are missing.
% If so, that function will print a results table.
[gook,R] = print_banner(R,wh_regressor,wh_field,cortype);
if ~gook, return, end

% -------------------------------------------------------------------
% Get significant clusters
% -------------------------------------------------------------------
[A,R] = robust_nonparam_get_sig_clusters(R,wh_regressor,doextended);


% -------------------------------------------------------------------
% Display images on orthviews
% Setup to print tables
%
% Create cl_extent, dat_extent (pos or neg sig. image vectors)
% dat, data at all thresholds in columns of image vector (for tables)
%
% Other useful outputs not used directly: cl_uncor, cl_uncor_neg
% -------------------------------------------------------------------
showrois = 1;

[dat,dat_pos,dat_neg,A] = robust_nonparam_orthviews(R,wh_regressor,overlay,showrois,A);

dat_extent = sum(dat,2);
cl_extent = iimg_indx2clusters(dat_extent,R.volInfo);

cl_extent_pos = iimg_indx2clusters(sum(dat_pos,2),R.volInfo);
cl_extent_neg = iimg_indx2clusters(sum(dat_neg,2),R.volInfo);



% ---------------------------------------------------------------------
% Make Tables
% ---------------------------------------------------------------------
uthr = [.005 .01 .05];   % uncorrected thresholds for display
dat = [A.p_sig>0 dat];
robust_nonpar_cluster_tables(R,doextended,wh_regressor,A.cl_all,A.cl_sig,A.cl_pos,A.cl_neg,A.rfieldname,A.pfieldname,uthr,dat,cl_extent_pos,cl_extent_neg)



domont = input('Make montage? ');
if domont
    montage_clusters(overlay,A.cl_all,cl_extent_pos,cl_extent_neg,A.cl_pos,A.cl_neg,{'g' 'r' 'c' 'y' 'b'},'nooverlap');
end


return







function [gook,R] = print_banner(R,wh_regressor,wh_field,cortype)
gook = 1;

if ~isfield(R,'svc') || ~isfield(R,'mapwise')
    R = robust_nonparam_results(R);
end

fprintf(1,'\n* ========================================================================')
fprintf(1,'\n* robust_nonparam_displayregions                         ')

switch wh_field
    case 't'
        fprintf(1,'\n* t-test for Regressor: %s, Correction: %s',R.names{wh_regressor},cortype)
    case 'f'
        fprintf(1,'\n* t-test for Regressor: %s, Correction: %s',R.names{wh_regressor},cortype)
    otherwise
        fprintf(1,'\n* Error: wh_field must be t or f\n')
        gook = 0;
end

fprintf(1,'\n* Tail is: %s',R.mytail)
fprintf(1,'\n* Corrected alpha level is: %3.2f',R.mapwise.alpha)
fprintf(1,'\n* ========================================================================\n')

if ~( strcmp(cortype,'svc') || strcmp(cortype,'mapwise') )
    fprintf(1,'\n* Error: cortype must be svc or mapwise\n')
    gook = 0;
end

return






% these subfunctions are no longer needed because they're in
% % robust_nonparam_get_sig_clusters
% 
% % -------------------------------------------------------------------
% % DATA: check for R.Y and load if necessary.  Check image names
% % -------------------------------------------------------------------
% 
% function [R,Y] = load_data_subfunction(R)
% sprintf('Loading data. ');
% 
% if isfield(R,'Y') && ~isempty(R.Y)
%     disp('R.Y data field found.  Using it.');
%     Y = R.Y;
% 
%     if ~exist(deblank(R.maskname),'file')
%         fprintf(1,'\n* Error: R.maskname is not a valid file.\nImage: %s\nPlease select.\n',R.maskname)
%         R.maskname = spm_get(1,'*','Select new mask image');
%     end
% 
% else
%     % check image names and get data
% 
%     if ~exist(deblank(R.image_names(1,:)),'file')
%         fprintf(1,'\n* Error: R.image_names are not valid files.  Assuming names are correct but path is wrong!\n')
%         R.image_names = filename_get_new_root_dir(R.image_names, ...
%             spm_get(-1,'*','Select root dir for files'),2);
%     end
% 
%     if ~exist(deblank(R.maskname),'file')
%         fprintf(1,'\n* Error: R.maskname is not a valid file.\nImage: %s\nPlease select.\n',R.maskname)
%         R.maskname = spm_get(1,'*','Select new mask image');
%     end
% 
%     % load data and save in R
%     Y = iimg_get_data(R.maskname,R.image_names);
%     R.Y = Y;
% end
% if ~isfield(R,'volInfo'), R.volInfo = iimg_read_img(R.maskname); end
% 
% return
% 
% % -------------------------------------------------------------------
% % DATA: max partial correlation in a cluster (or set of voxels generally)
% % -------------------------------------------------------------------
% 
% function [p,r] = get_max_partial_cor(X,Y,vox_index,wh_regressor)
% 
% warning('off','stats:statrobustfit:IterationLimit');
% 
% dat = Y(:,vox_index); % for mean of region: mean(R.Y(:,in_cluster == j), 2);
% 
% % get max partial correlation and min p
% % partialcor uses robust by default
% nvox = size(dat,2);
% rrob = zeros(nvox,1);
% prob = rrob;
% for vox = 1:nvox
%     [x,y,rrob(vox),prob(vox)] = partialcor(X,dat(:,vox),wh_regressor);
% end
% 
% p = min(prob);
% p = p(1);
% r = rrob(prob == p);
% r = r(1);
% 
% warning on
% return
% 
% 

