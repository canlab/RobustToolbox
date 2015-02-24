function [cl_me_pos, cl_me_neg, ...
    cl_mepos_covpos, cl_mepos_covneg, cl_meneg_covpos, cl_meneg_covneg] ...
    = robust_act_corr_intersection_2008(me_pthresh, cov_pthresh, size_thresh)
%robust_act_corr_intersection_2008(me_pthresh, cov_pthresh, size_thresh)
%
% Documentation not complete. 
% Tor Wager, Nov 2008
%
% Example:
% [cl_me_pos, cl_me_neg, ...
%     cl_mepos_covpos, cl_mepos_covneg, cl_meneg_covpos, cl_meneg_covneg] ...
%     = robust_act_corr_intersection_2008(me_pthresh, cov_pthresh)


me_image = 'rob_p_0001.img';       % the input image used to define the mask
me_sign_image = 'rob_tmap_0001.img';
%me_pthresh = .005;               % the p-value threshold used to define the mask
me_mask = 'main_effect_mask.img';  % the output image

cov_image = 'rob_p_0002.img';      % the covariate image we want to mask
%cov_pthresh = .05;               % the p-value threshold used to define the mask
cov_sign_image = 'rob_tmap_0002.img';

% me_mask = 'main_effect_mask.img';  % the output image

%size_thresh = 10;

cov_pvals_me_mask = mask_image(cov_image, me_image, me_mask, 'maxmask', me_pthresh);
% These are only covariate p-values with significant group effects
cov_pvals_me_mask(cov_pvals_me_mask == 0) = NaN;

% Cov p-values with significant effects in both me and cov
cov_pvals_masked = cov_pvals_me_mask;
cov_pvals_masked(cov_pvals_masked > cov_pthresh) = NaN;

% now we need to get information about the signs of each effect
me_tvals = spm_read_vols(spm_vol(me_sign_image));
cov_tvals = spm_read_vols(spm_vol(cov_sign_image)); 

V = spm_vol(me_mask);

% main effects clusters

cl_me_pos = mask2clusters(double(cov_pvals_me_mask > 0 & me_tvals > 0), V);
cl_me_neg = mask2clusters(double(cov_pvals_me_mask > 0 & me_tvals < 0), V);

% covariate clusters
cl_mepos_covpos = mask2clusters(cov_pvals_masked .* (me_tvals > 0) .* (cov_tvals > 0), V);
cl_mepos_covneg = mask2clusters(cov_pvals_masked .* (me_tvals > 0) .* (cov_tvals < 0), V);

cl_meneg_covpos = mask2clusters(cov_pvals_masked .* (me_tvals < 0) .* (cov_tvals > 0), V);
cl_meneg_covneg = mask2clusters(cov_pvals_masked .* (me_tvals < 0) .* (cov_tvals < 0), V);

% 
% % impose size threshold
cl_me_pos(cat(1, cl_me_pos.numVox) < size_thresh) = [];
cl_me_neg(cat(1, cl_me_neg.numVox) < size_thresh) = [];

cl_mepos_covpos(cat(1, cl_mepos_covpos.numVox) < size_thresh) = [];
cl_mepos_covneg(cat(1, cl_mepos_covneg.numVox) < size_thresh) = [];

cl_meneg_covpos(cat(1, cl_meneg_covpos.numVox) < size_thresh) = [];
cl_meneg_covneg(cat(1, cl_meneg_covneg.numVox) < size_thresh) = [];

if length(cl_me_pos) == 0, cl_me_pos = []; end
if length(cl_me_neg) == 0, cl_me_neg = []; end
if length(cl_mepos_covpos) == 0, cl_mepos_covpos = []; end
if length(cl_mepos_covneg) == 0, cl_mepos_covneg = []; end
if length(cl_meneg_covpos) == 0, cl_meneg_covpos = []; end
if length(cl_meneg_covneg) == 0, cl_meneg_covneg = []; end

fprintf('Group avg Positive (Top Orthviews panel): %3.0f regions\n', length(cl_me_pos));
fprintf('Group avg Negative (Bottom Orthviews panel): %3.0f regions\n', length(cl_me_neg));

fprintf('Group avg Positive, Covariate Positive (Red): %3.0f regions\n', length(cl_mepos_covpos));
fprintf('These regions show A > B overall, and greater A > B with positive covariate values.\n\n');

fprintf('Group avg Positive, Covariate Negative (Blue): %3.0f regions\n', length(cl_mepos_covneg));
fprintf('These regions show A > B overall, and reduced A > B with positive covariate values.\n\n');

fprintf('Group avg Negative, Covariate Positive (Blue): %3.0f regions\n', length(cl_meneg_covpos));
fprintf('These regions show B > A overall, and reduced B > A with positive covariate values.\n\n');

fprintf('Group avg Negative, Covariate Negative (Red): %3.0f regions\n', length(cl_meneg_covneg));
fprintf('These regions show B > A overall, and greater B > A with positive covariate values.\n\n');

%%
overlay = which('spm2_single_subj_T1_scalped.img');
nimgs = 2;
spm_check_registration(repmat(overlay, nimgs, 1));

cluster_orthviews(cl_me_pos, {[0 1 0]}, 'solid', 'handle', 1, 'add');
cluster_orthviews(cl_me_neg, {[0 1 0]}, 'solid', 'handle', 2, 'add');

% kludge: do this many times so colors come out
for i = 1:8
    cluster_orthviews(cl_mepos_covpos, {[1 0 0]}, 'solid', 'handle', 1, 'add');
    cluster_orthviews(cl_mepos_covneg, {[0 .3 1]}, 'solid', 'handle', 1, 'add');

    cluster_orthviews(cl_meneg_covpos, {[0 .3 1]}, 'solid', 'handle', 2, 'add');
    cluster_orthviews(cl_meneg_covneg, {[1 0 0]}, 'solid', 'handle', 2, 'add');

end

spm_orthviews_name_axis('Group avg (A>B, Pos)', 1);
spm_orthviews_name_axis('Group avg (B>A, Neg)', 2);
