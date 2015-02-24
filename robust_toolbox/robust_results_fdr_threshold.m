function SETUP = robust_results_fdr_threshold(corr_type, varargin)
% SETUP = robust_results_fdr_threshold(corr_type, ['mask', mask image defining regions], ['images', p-value images])
%
% Multiple comparisons correction tool for mediation_brain
% Calculates corrected threshold and saves in mediation_SETUP.mat, SETUP variable
% Returns SETUP with corrected threshold.
% You can then use this threshold in mediation_brain_results to get results
% at corrected p-value thresholds.
%
% input: corr_type
% Must be one of the following:
% -------------------------------
% 'fdr'
%  False Discovery Rate correction across multiple effects
% (multiple images).  The idea is that you can find the threshold that
% controls the overall FDR in a robust regression analysis, including all covariate images.  
% This threshold is a single threshold that provides control across the
% different tests.
%
% P-VALUE IMAGES TO USE
% ------------------------
% This function automatically uses all p-value images from a robust results
% directory
% You can enter your own image or set of images as well by using the 
% optional inputs:
% 'images', followed by a string matrix of p-value images
%
% SEARCH AREA
% ------------------------
% The search area is defined by default as the area in mask.img in the
% mediation results directory.  This image is written automatically when
% mediation_brain analyses are run, but you have to specify one for robust reg.
% You can enter your own mask by using the optional inputs:
% 'mask', followed by the mask image defining the regions
% This provides facility for ROI-based correction.
%
% Examples:
% ------------------------
% In a mediation directory, type:
% SETUP = robust_results_fdr_threshold('fdr');
%
% To calculate FDR for only a single effect, e.g., the Path b effect in the
% example below, enter the p-value image name you'd like to use.
% You can actually do this for ANY p-value image, not just mediation analysis images:
% SETUP = mediation_brain_corrected_threshold('fdr', 'mask', mask, 'imgs', 'M-Y_pvals.img');
%
% Tor Wager, Dec 3, 2008

% Set up optional inputs and defaults
mask = 'mask.img';
imgs = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            % functional commands
            case 'mask', mask = varargin{i+1}; varargin{i + 1} = [];
            case {'images', 'imgs'}, imgs = varargin{i+1}; varargin{i + 1} = [];

            otherwise, warning('robfit:badInput', ['Unknown input string option:' varargin{i}]);
        end
    end
end

switch lower(corr_type)
    case 'fdr'
        SETUP = run_fdr();

    otherwise
        error('Unknown correction type. See help.');
end

% do not save, because we may want to run results with different masks,
% etc.
%save mediation_SETUP -append SETUP


% --------------------------
% *
% INLINE FUNCTIONS
% *
% --------------------------

    function SETUP = run_fdr

        % Load SETUP
        % Choose sensible images to get FDR correction over
        % (Ignore images that do not vary across mediation tests.)

        load(fullfile(pwd, 'SETUP.mat'),  'SETUP');

        if isempty(imgs)

           imgs = dir('rob_p_*.img'); 
           imgs = char(imgs.name);

        end

        disp('Calculating FDR threshold across family of tests in these images:')
        disp(imgs)

        if ~exist(mask, 'file'), error('Cannot find mask image file.'); end
        maskInfo = iimg_read_img(mask, 2);

        pvals = iimg_get_data(maskInfo, imgs);
        fdr_p_thresh = FDR(pvals(:), .05);
        
        if isempty(fdr_p_thresh), fdr_p_thresh = -Inf; end
        
        fprintf('Total p-values: %7d\n', length(pvals(:)));

        fprintf('Combined FDR threshold is %7f\n', fdr_p_thresh);

        disp('Returning output in SETUP.fdr_p_thresh (you must save this yourself in the .mat file)');

        SETUP.fdr_p_thresh = fdr_p_thresh;
        
        disp(' ');
        disp('Getting FDR thresholds for individual images:');

        for i = 1:size(imgs, 1)
            fdr_p_thresh = FDR(pvals(i, :), .05);
            
            if isempty(fdr_p_thresh), fdr_p_thresh = Inf; end
            
             fprintf('%s\t FDR (q < .05) = %7f\n',imgs(i,:),  fdr_p_thresh);
             
             SETUP.fdr_p_thresh_indiv(i) = fdr_p_thresh;
        end
        
        disp('Returning output in SETUP.fdr_p_thresh_indiv (you must save this yourself in the .mat file)');
        
%         create_figure('FDR'); plot(sort(pvals))
%         plot_horizontal_line(.001)
%         plot( .05 * (1:length(pvals(:))) ./ length(pvals(:)), 'r')

    end


end
