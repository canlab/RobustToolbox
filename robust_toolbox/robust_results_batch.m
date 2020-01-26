%% Robust regression report

%% About this report
%
% This script displays a summary of robust regression results
% Stored in a CANLab robust regression folder
% 
% It can be used with the publish() command to generate a published html
% report.  e.g.,
%
%%
%  cd(my robust results directory)
%  publish('robust_results_batch')
% 

%%
% publish.m is a Matlab report-generating function that can print to PDF,
% HTML, Powerpoint, and other formats. See |help publish| for more.
%
% To create reports and save them in the current robust results dir, try:
%
%  publish_robust_regression_report
%
% This script can also be customized to create different types of reports

% Load files from current directory
% -----------------------------------------------------------------------

[trob, names, mask_obj, nsubjects, weights, SETUP] = robust_reg_load_files_to_objects(pwd);

% Load files from a CANLab robust regression directory into objects
% Objects contain information needed for results display and tables.
%
% trob      statistic_image object with one image per contrast (the model intercept is the first image)
% names     cell array of names of each image (intercept, regressor 1, etc.)
% mask_obj  image defining the set of voxels analyzed
% nsubjects image summarizing number of subjects with valid data in each voxel
% weights   fmri_data object with weight maps for each subject (i.e., input image)
% SETUP     struct saved in directory, including design matrix SETUP.X

% Preliminaries
% -----------------------------------------------------------------------

% Number of images, including intercept
k = size(trob.dat, 2);      

% Display helpers
dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

% Set up figure

create_figure('slices'); axis off
o2 = canlab_results_fmridisplay;

%% Mask of in-analysis voxels and number of participants
% This montage shows the voxels analyzed in green

o2 = addblobs(o2, region(mask_obj), 'color', [0 .7 0], 'trans');
o2 = o2.title_montage(5, 'Analysis mask');
drawnow, snapnow

o2 = removeblobs(o2);
o2 = addblobs(o2, region(nsubjects), 'mincolor', [1 1 1], 'maxcolor', [0 0 1], 'trans');
o2 = o2.title_montage(5, 'Coverage (number of participants)');
o2 = legend(o2);

drawnow, snapnow

fprintf('Modal value for number of participants: %d\n', mode(nsubjects.dat));

%% Plot of subject weights 
% This plot shows a summary of the weight maps for each subject
% Weights of less than one indicate that a subject is down-weighted for
% the low-weight voxel. Weights of zero indcate the subject was dropped entirely for
% that voxel.
%
% A subject with low weights across large areas of the brain is an outlier
% across many brain areas.
% 
% See |help fmri_data.plot| and the basic plot walkthrough on canlab.github,io 
% for more information about this plot
%
% In this plot, the 'carpet plot' of weights across all subjects is
% currently mean-zeroed across the whole dataset. The relative weights 
% provide the key information.

plot(weights);
drawnow, snapnow

% slices(weights, 'orientation', 'axial');


%% Unthresholded results maps
% This plot shows the unthresholded t-maps
%
% The first image is the intercept in a standard CANlab robust regression analysis
% If regressors are mean-centered, the intercept can be interpreted as the
% activity for the average subject (when regressor values are all 0)
%
% Controlling for regressors (of interest or nuisance covariates) can be
% very helpful in reducing known sources of error when assessing group-mean
% activation.

create_figure('robust t maps'); axis off
o2 = canlab_results_fmridisplay([], 'multirow', k);

for i = 1:k
    
    t = get_wh_image(trob, i);
    
    o2 = addblobs(o2, region(t), 'trans', 'wh_montages', 2*i - 1: 2*i);  % 2 montages in registry per slice display
    o2 = title_montage(o2, 2*i, names{i});
end

drawnow, snapnow

%% FDR-corrected results maps and tables for each regressor
%
% This plot shows False Discovery Rate corrected t-maps for each regressor
% Apply gray matter mask before FDR correction

trob = apply_mask(trob, which('gray_matter_mask.img'));

o2 = removeblobs(o2);

for i = 1:k
    
%     disp(' ')
%     printhdr(sprintf('FDR-corrected: %s', names{i}))
%     disp(' ')
%     
    t = get_wh_image(trob, i);
    
    t = threshold(t, .05, 'fdr');
    o2 = addblobs(o2, region(t), 'trans', 'wh_montages', 2*i - 1: 2*i);
    o2 = title_montage(o2, 2*i, names{i});
    
end

drawnow, snapnow

%% FDR-corrected results surfaces and tables for each regressor
%
% We apply a gray matter mask before FDR correction and correct within
% gray-matter voxels. The mask is fairly liberal to avoid excluding brain
% tissue of interest. (It assumes data are in MNI space).
%
% Correcting within gray matter applies the multiple comparisons correction
% only within areas where we would plausibly expect signai. i.e., we are
% not penalized for searching for activity in the ventricles.
%
% The table() function generates a table with region and network names
% labeled using a standard CANlab atlas object. See the walkthrough on
% atlas objects in canlab.github.io for more information. 
% The default atlas at time of writing this is the canlab2018_2mm atlas.

for i = 1:k
    
    disp(' ')
    printhdr(sprintf('FDR-corrected: %s', names{i}))
    disp(' ')
    
    t = get_wh_image(trob, i);
    
    t = threshold(t, .05, 'fdr');
    
    % Table
    r = region(t);
    r = table(r, 'nolegend');
    
    % Surfaces and slices
    
    figure; 
    surface_handles = surface(t);
    
    o2_surf = canlab_results_fmridisplay(r, 'montagetype', 'full', 'nooutline');
    
    drawnow, snapnow
    
end



