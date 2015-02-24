function [A,R] = robust_nonparam_get_sig_clusters(R,wh_regressor,doextended)
% function [A,R] = robust_nonparam_get_sig_clusters(R,wh_regressor,doextended)
%
% A has:
%A.cl_all,A.cl_sig,A.wh_sig,A.p_sig,R,A.cl_pos,A.cl_neg,A.p_pos,A.p_neg,
%A.pfieldname,A.rfieldname
%
% Get cluster significance information: indices of significant values and
% clusters structures
% from robust_nonpar structure R, given contrast and effect of interest
%
% tor wager

% make sure we have data, if we're doing partial correlations and stuff
if doextended
    R = load_data_subfunction(R);
end

% -------------------------------------------------------------------
% Sizes and data
% -------------------------------------------------------------------

numc = R.volInfo.c;         % number of clusters
n = R.volInfo.nvox;
[k,v] = size(R.correct.t);

% all voxels in cluster
allvox = ones(v,1); % for actual t-values : R.correct.t(wh_regressor,:)';
A.cl_all = iimg_indx2clusters(allvox,R.volInfo);


A.cl_sig = cell(1,numc);            % significant clusters
A.wh_sig = logical(zeros(1,numc));  % which clusters
A.p_sig = (zeros(v,1));
A.p_pos = A.p_sig;
A.p_neg = A.p_sig;


% -------------------------------------------------------------------
% Loop through ROIs
% -------------------------------------------------------------------

for c = 1:numc

    % this is signed, pos or neg.
    sigvox = R.svc.tsig{c}(:,wh_regressor);

    whpos = find(sigvox > 0);
    whneg = find(sigvox < 0);

    % all sig. voxels
    whvox = abs(sigvox) > 0;

    % for t-values (signed)
    sigvox = R.correct.t(wh_regressor,:)' .* whvox;

    % save vector of all significant results for visualization (p-values)
    A.p_sig(whvox) = R.correct.p(wh_regressor,whvox);

    % significant voxels only
    [A.cl_sig{c},in_cluster] = iimg_indx2clusters(sigvox,R.volInfo,R.svc.tcrit(c,wh_regressor));

    % positive and negative clusters and p-values
    sigpos = sigvox;
    sigpos(sigpos < 0) = 0;

    A.p_pos(whpos) = R.correct.p(wh_regressor,whpos);

    A.cl_pos{c} = iimg_indx2clusters(sigpos,R.volInfo,R.svc.tcrit(c,wh_regressor));

    signeg = sigvox;
    signeg(signeg > 0) = 0;
    A.p_neg(whneg) = R.correct.p(wh_regressor,whneg);

    A.cl_neg{c} = iimg_indx2clusters(signeg,R.volInfo,R.svc.tcrit(c,wh_regressor));


    % save additional info: from_cluster, and partial corrs if requested
    % -------------------------------------------------------------------
    if isfield(R,'Y'), str = sprintf('Getting partial correlations for cluster %03d',c); fprintf(1,str); end

    [A.cl_sig,A.pfieldname,A.rfieldname] = get_partial_correlations(A.cl_sig,c,in_cluster,R,wh_regressor);

    A.cl_pos = get_partial_correlations(A.cl_pos,c,in_cluster,R,wh_regressor);
    A.cl_neg = get_partial_correlations(A.cl_neg,c,in_cluster,R,wh_regressor);


    A.wh_sig(c) = ~isempty(A.cl_sig{c});

    if isfield(R,'Y'), erase_string(str); end
end

A.cl_all = A.cl_all(A.wh_sig);
A.cl_sig = cat(2,A.cl_sig{:});

A.cl_pos = cat(2,A.cl_pos{:});
A.cl_neg = cat(2,A.cl_neg{:});

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
if ~isfield(R,'volInfo'), R.volInfo = iimg_read_img(R.maskname); end

return


% -------------------------------------------------------------------
% DATA: from_cluster and max partial corrs for a series of clusters in cell
% c
% -------------------------------------------------------------------

function [cl_sig,pfieldname,rfieldname] = get_partial_correlations(cl_sig,c,in_cluster,R,wh_regressor)

% names
k = size(R.X,2);
pfieldname = cell(1,k);
rfieldname = cell(1,k);

for myregressor = 1:k

    pfieldname{myregressor} = [R.names{myregressor} '_p'];
    rfieldname{myregressor} = [R.names{myregressor} '_partialr_or_mean'];
    if any(pfieldname{myregressor} == ' '), error('No spaces or special chars allowed in R.names!!'); end
end

% correlations
for j = 1:length(cl_sig{c})
    cl_sig{c}(j).from_cluster = c;

    cl_all(c).threshold = R.svc.tcrit(c,wh_regressor);
    cl_sig{c}(j).threshold = R.svc.tcrit(c,wh_regressor);

    if isfield(R,'Y')

        % get max effect or partial correlation for the effect of
        % interest
        %[cl_sig{c}(j).p,cl_sig{c}(j).correl] = ...
        %    get_max_partial_cor(R.X,R.Y,in_cluster == j,wh_regressor);

        % get max effect or partial correlation for the other
        % predictors
%         for myregressor = 1:size(R.X,2)
% 
%             [cl_sig{c}(j).(pfieldname{myregressor}),cl_sig{c}(j).(rfieldname{myregressor})] = ...
%                 get_max_partial_cor(R.X,R.Y,in_cluster == j,myregressor);
% 
%             
%         end
        
        % get max effect or partial correlation for all
        % predictors
        [rrob,prob,Ymaxeffect,whY,Xadj,Yadj] = robust_max_partial_corr(R.X,R.Y(:,in_cluster==j),0);
        for myregressor = 1:size(R.X,2)
            cl_sig{c}(j).(rfieldname{myregressor}) = rrob(myregressor);
            cl_sig{c}(j).(pfieldname{myregressor}) = prob(myregressor);
        end

        cl_sig{c}(j).Ymaxeffect = Ymaxeffect;
        cl_sig{c}(j).whY = whY;
        
        % get max effect or partial correlation for the effect of
        % interest
        cl_sig{c}(j).p = prob(wh_regressor);
        cl_sig{c}(j).correl = rrob(wh_regressor);
    end
end

return


% -------------------------------------------------------------------
% DATA: max partial correlation in a cluster (or set of voxels generally)
% -------------------------------------------------------------------

% function [p,r] = get_max_partial_cor(X,Y,vox_index,wh_regressor)
% 
% warning('off','stats:statrobustfit:IterationLimit');
% 
% dat = Y(:,vox_index); % for mean of region: mean(R.Y(:,in_cluster == j), 2);
% 
% 
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



function erase_string(str1)
fprintf(1,repmat('\b',1,length(str1))); % erase string
return
