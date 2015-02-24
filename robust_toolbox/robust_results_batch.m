function robust_results_batch(varargin)
%robust_results_batch(varargin)
%
% This function runs and saves tables, images, and clusters
% for all effects of a robust regression model
%
% You need Tor's mediation toolbox on your path.
%
% You also need to enter optional inputs as you would to
% mediation_brain_results.
% Enter any arguments to mediation_brain_results:
% e.g., ... 'thresh', [.005 .01 .05], 'size', [3 1 1], 'prune', 'overlay',
% EXPT.overlay, 'mask', EXPT.mask);
%
% You MUST enter a mask image name. robfit does not currently save a
% default one, and mediation_brain_results requires one.
%
% Tor Wager, March 2010
% Updated: Aug 2012
%       - compatibility with publish
%       - new default slice views, compatible with canlab_results_fmridisplay slice format
%       - added some diagnostic plots (qchist) as standard output
%
% Optional inputs:
% -------------------------------------------------------------------------
%             case 'mask', mask = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
%             case 'overlay', overlay = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
%                 
%             case 'thresh', pthresh = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
%             case 'size', kthresh = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
%                 
%             case 'reversecolors', reversecolors = 1;
%                 
%             case 'skipmask', skipmask = 1;
%                 
%             case {'rob', 'rob0', 'rob1', 'rob2', 'rob3', 'rob4'}
%                 
% Examples:
% -------------------------------------------------------------------------
% robust_results_batch('overlay', EXPT.overlay, 'mask', EXPT.mask, 'thresh', [.001 .005 .05], 'size', [5 1 1], 'prune');
% follow up with surface plots and a surface movie:
% load cl_rob_p_0002_001_005_05_k5_1_1_prune.mat  % this is for the covariate effect
% mediation_brain_surface_figs(clpos, clneg);

mask = 'rob_beta_0001.img';
overlay = [];
pthresh = [.001 .01 .05];  %[.005 .01 .05]; %[.001 .005 .01]; %'fdr';            % specify 'fdr' or a series of p thresholds
kthresh = [3 1 1];          % specify a series of extent thresholds
reversecolors = 0;
orig_varargin = varargin;
overlay = which('SPM8_colin27T1_seg.img');
skipmask = 0;

% special inputs that we do extra things with;
% all optional inputs are passed to mediation_brain_results
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'mask', mask = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
            case 'overlay', overlay = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
                
            case 'thresh', pthresh = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
            case 'size', kthresh = varargin{i+1}; varargin{i} = []; varargin{i+1} = [];
                
            case 'reversecolors', reversecolors = 1;
                
            case 'skipmask', skipmask = 1;
                
            case {'rob', 'rob0', 'rob1', 'rob2', 'rob3', 'rob4'}
                warning('Do not enter rob* inputs to robust_results_batch');
                varargin{i} = []; orig_varargin{i} = [];
                
                % omit warning because we may have other valid inputs to pass forward
                %otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isempty(mask), error('You must enter a mask image or omit mask for the default: rob_beta_0001.img.'); end
if isempty(overlay), error('You must enter an overlay (anatomical) image, or omit overlay for the default: SPM8_colin27T1_seg.img.'); end

currdir = pwd;
setupfile = fullfile(currdir, 'SETUP.mat');
if ~exist(setupfile, 'file'); error('SETUP.mat does not exist in current dir; go to valid robust results directory.'); end


% ---------------------------------------------------------
% Report header
% ---------------------------------------------------------

load(setupfile)
ncols = size(SETUP.X, 2);

reportfile = fullfile(pwd, 'robust_results_batch_report.txt');
diary(reportfile)

hdrline = '==========================';
sepline = '________________________________';

fprintf('%s\n=\n= robust_results_batch.m\n=\n%s\n', hdrline, hdrline)
fprintf('Directory: %s\nDate report generated: %s\n%s\n', currdir, date, sepline);

%pause(1)
fprintf('\n\nModel info:\n%s\nImages and regressors (intercept + %3.0f covs)\n', sepline, ncols - 1);
imgs = mat2cell(SETUP.files, ones(size(SETUP.files, 1), 1), size(SETUP.files, 2));
print_matrix(SETUP.X, [], imgs);

fprintf('\n%s\n', hdrline)

fprintf('\n\nResults info:%s\nInput options:\n', sepline);
for i = 1:length(orig_varargin)
    if ischar(orig_varargin{i})
        fprintf('%s\t', orig_varargin{i});
    elseif ishandle(orig_varargin{i})
    else
        print_matrix(orig_varargin{i})
    end
end

fprintf('\nResults inputs\n%s\n', sepline);
fprintf('\nMask: %s\nOverlay: %s\n', mask, overlay);
fprintf('\nP-value thresholds: ');  print_matrix(pthresh);
fprintf('\nExtent thresholds: ');  print_matrix(kthresh);

% ---------------------------------------------------------
% Mask figure
% ---------------------------------------------------------

% clear existing windows. No longer necessary, but good anyway.
figh = findobj('Tag', 'montage_axial');
if ishandle(figh), clf(figh); end

figh = findobj('Tag', 'montage_coronal');
if ishandle(figh), clf(figh); end

figh = findobj('Tag', 'montage_sagittal');
if ishandle(figh), clf(figh); end

if ~skipmask
    
    cltmp = mask2clusters(mask);
    cluster_orthviews(cltmp, {[0 1 0]}, 'trans');
    cluster_orthviews_montage(6, 'axial', overlay, 'onerow');
    %snapnow
    %mask_fig_handle = montage_clusters(mask);
    mask_fig_handle = findobj('Tag', 'montage_axial');
    set(mask_fig_handle, 'Name', 'Mask image');
    scn_export_papersetup(500); saveas(mask_fig_handle, fullfile(currdir, 'Results_Mask.png'));
    
    cluster_orthviews_montage(10, 'sagittal', overlay, 'onerow');
    
    f = findobj('Tag', 'Graphics');
    set(f, 'Visible', 'off');
    
    snapnow
end

% ---------------------------------------------------------
% Histogram+ on input data
% ---------------------------------------------------------
if ~isempty(check_valid_imagename(SETUP.files, 0))
    
    fprintf('%s\nImage dataset visualization and diagnostics\n%s\n', sepline, sepline);
    qchist(SETUP.files);
    sfig = findobj('Tag', 'Graphics');
    set(sfig, 'Visible', 'off');
    snapnow
else
    disp('Cannot find images in SETUP.files. Skipping QC histograms.')
end

% ---------------------------------------------------------
% Case (subject) weights
% ---------------------------------------------------------

fprintf('%s\nCase (subjects) weights: Low indicates likely outliers\n%s\n', sepline, sepline);
dat = fmri_data('weights.img');
plot(dat)
sfig = findobj('Tag', 'Graphics');
set(sfig, 'Visible', 'off');
snapnow


% ---------------------------------------------------------
% Report results info, including FDR correction
% ---------------------------------------------------------

fprintf('\nFDR correction info (FDR not necessarily used here; just for info):%s\n', sepline);
SETUP = robust_results_fdr_threshold('fdr', 'mask', mask);

fprintf('\n%s\nColors: ', sepline);
if reversecolors, fprintf('Reversed\n%s\n', hdrline);
else fprintf('Normal\n%s\n', hdrline);
end

if reversecolors
    negcolors = { [1 1 0] [1 .5 0] [1 .3 .3] };
    poscolors = { [0 0 1] [0 .5 1] [.3 .3 1] };
else
    poscolors = { [1 1 0] [1 .5 0] [1 .3 .3] };
    negcolors = { [0 0 1] [0 .5 1] [.3 .3 1] };
end

% ---------------------------------------------------------
% Run for each effect
% ---------------------------------------------------------

for i = 1:ncols
    
    robstring = ['rob' num2str(i - 1)];
    
    fprintf('%s\n=\n= Results for %s\n=\n%s\n', hdrline, robstring, hdrline)
    
    [clpos, clneg, clpos_data, clneg_data, clpp2, clnn2] = mediation_brain_results ...
        (robstring, 'thresh', pthresh, 'size', kthresh,  ...
        'tables', 'slices', 'save', 'overlay', overlay, 'mask', mask, 'poscolors', poscolors, 'negcolors', negcolors, varargin{:});
    
    f = findobj('Tag', 'Graphics');
    set(f, 'Visible', 'off');
    
    snapnow
end

diary off

end % function


