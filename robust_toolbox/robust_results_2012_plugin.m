function robust_results_2012_plugin(R)
% robust_results_2012_plugin(R)
%
% For robust nonparametric analysis:
% Load the variable R from the saved .mat file and run
%
% This script:
% 1) shows the mask
% 2) shows histograms of the observed and null max T-values
% 3) prints information about the whole-mask threshold
% 4) writes thresholded (whole-mask) .img files to disk in a subfolder
% 
% It is designed to be called by a wrapper script invoked by the publish
% command, for saving results in HTML archives.

%% Mask used in analysis

% display mask and set up fmridisplay o2 object
% ---------------------------------------------------------
r = region(R.maskname);
o2 = canlab_results_fmridisplay(r);

set(o2.activation_maps{1}.blobhandles, 'FaceColor', [0 1 0], 'FaceAlpha', .4)
set(o2.activation_maps{2}.blobhandles, 'Color', [0 .4 0])

%% Null and observed t-value histograms

% threshold for whole-mask correction and histogram
% ---------------------------------------------------------

k = size(R.maxt, 1);
create_figure('Histograms', 1, k);

for i = 1:k
    
    create_figure('Histograms', 1, k, 1);
    subplot(1, k, i);
    
    t = R.correct.t(i, :)';
    nullt = R.maxt(i, :)';
    nulltn = -nullt; % for both + and - tails
    
    x = linspace(min([nullt; nulltn; t]), max([nullt; nulltn; t]), 100);
    
    h1 = hist([nullt; nulltn], x);
    
    h2 = hist(t, x);
    
    han1 = bar(x, h1, 'FaceColor', [.5 .5 .5], 'EdgeColor', 'none');
    han2 = bar(x, h2, 'FaceColor', [0 0 1], 'EdgeColor', 'none');

    plot_vertical_line(prctile(nulltn, 5));
    plot_vertical_line(prctile(nullt, 95));
    
    patchhan = get(han2, 'Children');
    set(patchhan, 'FaceAlpha', .7);
    
    if i == 1
    legend([han1 han2], {'Null (max)' 'Observed'})
    end
    
    xlabel('T-value');
    ylabel('Frequency');
    
    title(R.names{i});
    
end

drawnow;
snapnow;

%% Whole-mask thresholded maps (FWE p < .05)

% threshold and reconstruct maps
% -------------------------------------------------------

% set up for writing images
dat = statistic_image('image_names', R.maskname);

outputdir = fullfile(pwd, 'thresholded_images');
if ~exist(outputdir, 'dir'), mkdir(outputdir); end

z = '------------------------------------------------------';
fprintf('%s\n%s\n', z, z);

fprintf('Thresholds (t-values) for p < .05 corrected across mask\n');

for i = 1:k
    
    t = R.correct.t(i, :)';
    nullt = R.maxt(i, :)';
    nulltn = -nullt; % for both + and - tails
    
    negthr = prctile(nulltn, 5);
    posthr = prctile(nullt, 95);
    
    % Table of thresholds and sig vox
    % -----------------------------------------------------------
    fprintf('%s\n%s\n%s\n', z, R.names{i}, z);

    fprintf('Param %3.0f - %s: t < %3.3f or t > %3.3f\t', i, R.names{i}, negthr, posthr)
    
    fprintf('Sig -:\t%3.0f\tvoxels\tSig +:\t%3.0f\tvoxels\n', sum(t < negthr), sum(t > posthr));
    

    % write map with sig. voxels
    tmap = dat;
    tmap.type = 'T';
    tmap.dat = t;
    tmap.sig = (t < negthr) | (t > posthr);
    
    outname = fullfile(outputdir, ['Thresholded_FWE05_' R.names{i} '.img']);
    
    tmap.fullpath = outname;
    write(tmap);
    
    % display on fmridisplay
    o2 = removeblobs(o2);
    
    if sum(tmap.sig > 0)
        cl = region(tmap);
        o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]});
    end
    
    o2 = addblobs(o2, r, 'maxcolor', [0 .2 0], 'outline');

    drawnow;
    snapnow;
end

fprintf('%s\n%s\n', z, z);

end % function


