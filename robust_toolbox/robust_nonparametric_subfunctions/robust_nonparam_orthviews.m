function [dat,dat_pos,dat_neg,A] = robust_nonparam_orthviews(R,wh_regressor,overlay,showrois,A)
% [dat,dat_pos,dat_neg,A] = robust_nonparam_orthviews(R,wh_regressor,overlay,showrois,[A])
%
% p_sig,cl_pos, and cl_neg inputs can be empty or missing, in which case
% robust_nonparam_get_sig_clusters is run to get it.
%
% A is "activation" structure with cl_pos and cl_neg in it; see
% robust_nonparam_get_sig_clusters.  If not input,
% robust_nonparam_get_sig_clusters is run here.
%
% Examples:
%  To run using saved clusters from previous
%  robust_nonparam_displayregions (get 3rd regressor stats):
%  [dat,dat_pos,dat_neg] = robust_nonparam_orthviews(R,3,EXPT.overlay,1,A);
%
%  To run without saved clusters, creating them here:
%  [dat,dat_pos,dat_neg,A] = robust_nonparam_orthviews(R,wh_regressor,overlay,showrois);
%
% tor wager

% wh_regressor = 3;
% overlay = EXPT.overlay;
% showrois = 1;

pvals = R.correct.p(wh_regressor,:)';
tvals = R.correct.t(wh_regressor,:)';


uthr = [.005 .01 .05];   % uncorrected thresholds for display
poscolors = [1 1 0; 1 .5 0; 1 .3 .3; .8 0 0; 0 1 0]; % corrected first; last is color for ROIs if transseed
negcolors = [0 0 1; 0 .5 1; .3 .5 1; .3 .3 1; 0 1 0];

% need all voxels even at lowest threshold
% dat needed as output for main program
% now handled below
%[cl_uncor,dat] = iimg_multi_threshold(pvals,'p','thresh',uthr,'size',[1 1 1], ...
%'pruneseed',p_sig,'volInfo',R.volInfo,'overlay',overlay,'colors',poscolors);


% just need p_sig
if nargin < 5
    A = robust_nonparam_get_sig_clusters(R,wh_regressor,0);
end


fprintf(1,'\n* Orthviews: All results within ROIs                  ')
% ---------------------------------------------------------------------
% Orthviews: Sig. within ROIs (all regions)
% ---------------------------------------------------------------------
% only in extended sig region; to get DAT
% this returns only areas within ROIs that are significant, because of
% the implicit mask

% Without regions of interest (already done above)
%[cl_uncor,dat] = iimg_multi_threshold(pvals,'p','thresh',uthr,'size',[1 1 1], ...
%'pruneseed',A.p_sig,'volInfo',R.volInfo,'overlay',overlay,'colors',poscolors);

fprintf('\nShowing results for uncorrected thresholds in orange/red for positive results');
fprintf('\nand blue/purple for negative results.')
fprintf('\nResults with corrected significance will be shown in yellow or dark blue.');
fprintf('\n')

% ---------------------------------------------------------------------

% Orthviews: Positive
% ---------------------------------------------------------------------
%
% anywhere within with ROIs
% this is enforced implicitly because of the ROI mask applied earlier
[A.cl_uncor_pos,dat_pos] = iimg_multi_threshold(pvals,'p','thresh',uthr,'size',[1 1 1], ...
    'volInfo',R.volInfo,'overlay',overlay, ...
    'colors',poscolors(2:end,:),'pos',tvals);

addstr = 'add';
if isempty(A.cl_uncor_pos{end}), addstr = 'noadd'; end

% ---------------------------------------------------------------------

% Orthviews: Negative
% ---------------------------------------------------------------------
[A.cl_uncor_neg,dat_neg] = iimg_multi_threshold(pvals,'p','thresh',uthr,'size',[1 1 1], ...
    'volInfo',R.volInfo,'overlay',overlay, ...
    'colors',negcolors(2:end,:),'neg',tvals,addstr);

% save dat for either pos or neg for later use
dat = dat_pos | dat_neg;

% ---------------------------------------------------------------------

% Orthviews: Significant within ROIs x pos. or neg. results
% ---------------------------------------------------------------------
fprintf('\nCorrected positive results: %3.0f separate regions', length(A.cl_pos));
fprintf('\nCorrected negative results: %3.0f separate regions', length(A.cl_neg));
fprintf('\n');

% significant with correction
if ~isempty(A.cl_pos)
    cluster_orthviews(A.cl_pos,{poscolors(1,:)},'solid','add');
end
if ~isempty(A.cl_neg)
    cluster_orthviews(A.cl_neg,{negcolors(1,:)},'solid','add');
end

if showrois
    fprintf(1,'Showing extent of ROIs with voxels sig. at all thresholds.');
    allsigseed = sum(dat,2);
    roivec = ones(size(pvals)); %~allsigseed;  % % all vox in mask  %R.volInfo.image_indx;
    
    roivec = iimg_cluster_prune(roivec,allsigseed,R.volInfo);
    roivec(allsigseed>0) = 0;
    cl_all2 = iimg_indx2clusters(roivec,R.volInfo);
    
    cluster_orthviews(cl_all2,{poscolors(end,:)},'trans','add');  %'overlay',overlay);
    
    makelegend({'Corrected (+)' 'p < .005' 'p < .01' 'p < .05' 'Corrected (-)' 'p < .005' 'p < .01' 'p < .05' 'ROI extent'},[poscolors(1:4,:); negcolors(1:4,:); poscolors(end,:)],1);
else
    makelegend({'Corrected (+)' 'p < .005' 'p < .01' 'p < .05' 'Corrected (-)' 'p < .005' 'p < .01' 'p < .05'},[poscolors(1:4,:); negcolors(1:4,:)],1);
end

% ---------------------------------------------------------------------

% Option to save figures here
% ---------------------------------------------------------------------

savefigs = input('Save PNG images of these clusters? (1/0) ');

if savefigs
    cl = merge_clusters(A.cl_pos,A.cl_neg); % preserve names if these are input
    cl = cluster_export_pngs(cl,1,overlay); % save names if entered here
    
    % transfer names
    if isfield(cl,'shorttitle')
        npos = length(A.cl_pos);
        nneg = length(A.cl_neg);
        for i = 1:npos
            A.cl_pos(i).shorttitle = cl(i).shorttitle;
        end
        for i = 1:nneg
            A.cl_neg(i).shorttitle = cl(i+npos).shorttitle;
        end
    end
end

% ---------------------------------------------------------------------

% Option to add extent of activation outside ROIs
% ---------------------------------------------------------------------

add_extent_outside_rois(R,uthr,overlay,poscolors,negcolors);

% ---------------------------------------------------------------------

% Option to save figures again with extent
% ---------------------------------------------------------------------

savefigs = input('Save PNG images of these clusters now? (1/0) ');

if savefigs
    cl = merge_clusters(A.cl_pos,A.cl_neg); % preserve names if these are input
    cl = cluster_export_pngs(cl,1,overlay);
    
    % transfer names
    if isfield(cl,'shorttitle')
        npos = length(A.cl_pos);
        nneg = length(A.cl_neg);
        for i = 1:npos
            A.cl_pos(i).shorttitle = cl(i).shorttitle;
        end
        for i = 1:nneg
            A.A.cl_neg(i).shorttitle = cl(i+npos).shorttitle;
        end
    end
end

return




% -------------------------------------------------------------------
% DISPLAY: add regions contiguous with (but outside) ROIs to color mapped
% orthviews
% -------------------------------------------------------------------


function add_extent_outside_rois(R,uthr,overlay,poscolors,negcolors)
doextent = input('Show extent of regions contiguous with but outside ROIs? (1/0) ');
if doextent
    z = zeros(R.volInfo.nvox,1); %whos z
    %z(R.volInfo.wh_inmask) = dat_extent;   % this line for getting only
    % clusters contiguous with sig. nonparam clusters
    
    z(R.volInfo.wh_inmask) = 1;             % to get entire ROI extent
    seeddat = z;
    
    switch spm('Ver')
        case 'SPM2'
            pimg = spm_get(1,'*img','Select brain-wise p-map to use');
            timg = spm_get(1,'*img','Select matching t-map to use');
            
        case 'SPM5'
            % spm_defaults is a function
            spm_defaults()
            pimg = spm_select(1,'image','Select brain-wise p-map to use');
            timg = spm_select(1,'image','Select matching t-map to use');
            
        case 'SPM8'
            % spm_defaults is a function
            spm_defaults()
            pimg = spm_select(1,'image','Select brain-wise p-map to use');
            timg = spm_select(1,'image','Select matching t-map to use');
            
        otherwise
            % unknown SPM
            disp('Unknown version of SPM!');
            spm_defaults()
    end
    
    
    
    % ROIs have smoothed weights.  Stats in p-map may not match those from
    % Take previous orthviews and add extent around it (anything contig.
    % with extended activation in ROI)
    
    iimg_multi_threshold(pimg,'thresh',uthr,'size',[1 1 1],'p','pruneseed',seeddat,'overlay',overlay, ...
        'pos',timg,'colors',poscolors,'add','hideseed');
    
    %cluster_orthviews(cl_sig,{[1 1 0]},'solid','add');
    
    iimg_multi_threshold(pimg,'thresh',uthr,'size',[1 1 1],'p','pruneseed',seeddat, ...
        'neg',timg,'colors',negcolors,'add','hideseed');
    
end

return

