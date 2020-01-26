% [clpos, clneg, dat, volInfo] = robust_results_threshold(pthr, kthr, [cmd strings])
%
% Threshold a p-value image and return values in a second image (designed
% for t-images, but can be anything.)
%
% Command strings
%  'mask', followed by mask image name
%  'overlay', followed by overlay image name
%  'p', followed by p-value image name
%  't', followed by t-value image name
%  'intercept', pimg = 'rob_p_0001.img'; timg = 'rob_tmap_0001.img';
%  'cov1', pimg = 'rob_p_0002.img'; timg = 'rob_tmap_0002.img';
%  'cov2', pimg = 'rob_p_0003.img'; timg = 'rob_tmap_0003.img';
%  'cov3', pimg = 'rob_p_0004.img'; timg = 'rob_tmap_0004.img';
%  'fmap', pimg = 'rob_pmap_full.img'; timg = 'rob_fmap_full.img';
%  'write', write thresholded output images to disk
%  'nodisplay', do not display results on brain
%  'fdr', do FDR correction
%   'handle', followed by which orthviews plot to use
%   case 'add', add to existing orthviews
%
% Examples
% -------------------------------------------------------------------------
% Default: use rob_p_0001.img and t_0001.img
% [clpos, clneg, dat, volInfo] = robust_results_threshold(.005, 10);
%
% Do not display orthviews, but write masks
% [clpos, clneg, dat, volInfo] = robust_results_threshold(.001, 10, 'nodisplay', 'write');
%
% Threshold correlation map (rob_p_0002.img) and mask with pos intercept activation
% [clpos, clneg, dat, volInfo] = robust_results_threshold(.005, 10, 'cov1', 'mask', '../mask_images/hr_corr_pos_fdr05.img');
%
% robust_results_threshold(.01, 5, 'cov1', 'overlay', EXPT.overlay);
%
% %% put results in second orthviews pane:
% robust_results_threshold(.01, 5, 'cov1', 'overlay', EXPT.overlay, 'handle', 2);

function [clpos, clneg, dat, volInfo] = robust_results_threshold(pthr, kthr, varargin)
    % --------------------------------------
    % * Set up arguments
    % --------------------------------------

    % defaults
    
    switch spm('Ver')
        case 'SPM2'
            % spm_defaults is a script
            disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');
            
        case 'SPM5'
            % spm_defaults is a function
            spm_defaults()
    end
    
    pimg = 'rob_p_0001.img';
    timg = 'rob_tmap_0001.img';
    writeoutput = 0;
    mask  = [];
    display = 1;
    clpos = []; clneg = [];
    dofdr = 0;
    overlay = [];
    orth_handle = 1;
    addstr = 'do not add.';
    
    % inputs
    for i = 1:length(varargin)
        arg = varargin{i};
        if ischar(arg)
            switch lower(arg)
                case 'mask', mask = varargin{i+1}; varargin{i+1} = [];
                case 'overlay', overlay = varargin{i+1}; varargin{i+1} = [];
                case 'p', pimg = varargin{i+1}; varargin{i+1} = [];
                case 't', timg = varargin{i+1}; varargin{i+1} = [];
                case 'intercept', pimg = 'rob_p_0001.img'; timg = 'rob_tmap_0001.img';
                case 'cov1', pimg = 'rob_p_0002.img'; timg = 'rob_tmap_0002.img';
                case 'cov2', pimg = 'rob_p_0003.img'; timg = 'rob_tmap_0003.img';
                case 'cov3', pimg = 'rob_p_0004.img'; timg = 'rob_tmap_0004.img';
                case 'cov4', pimg = 'rob_p_0005.img'; timg = 'rob_tmap_0005.img';
                        
                case 'fmap', pimg = 'rob_pmap_full.img'; timg = 'rob_fmap_full.img';
                case 'write', writeoutput = 1;
                case 'nodisplay', display = 0;
                case 'fdr', dofdr = 1;
                    
                case 'handle', orth_handle = varargin{i + 1};
                case 'add', addstr = 'add';
                    
                otherwise
                    warning(['Unrecognized command string : ' arg])
                    
            end
        end
    end

    % --------------------------------------
    % * Read data and apply mask, if specified
    % --------------------------------------

    % check for images in  current dir
    if ~exist(fullfile(pwd, pimg)) || ~exist(fullfile(pwd, timg))
        disp('Cannot find p- or t-image.  You must go to a valid robust0??? results directory.'); 
    end
    
    % return p-values for in-analysis voxels only
    if isempty(mask)
        [volInfo, dat] = iimg_read_img(pimg, 1);
    else
        [dat, volInfo] = iimg_mask(mask, pimg);
    end
    dat = dat(volInfo.wh_inmask);
    
    eligible_vox = ~isnan(dat) & (dat > 0);
    mp = min(dat(eligible_vox));


    % --------------------------------------
    % * Begin output report
    % --------------------------------------

    fprintf(1, '\n\n---------------------------------------------------------------------\n')
    fprintf(1, 'Robust Regression output report')
    fprintf(1, '\n---------------------------------------------------------------------\n\n')

    fprintf(1, 'p-value image thresholded: %s\n', pimg)
    fprintf(1, '\nIn-analysis voxels: %3.0f total\n', sum(eligible_vox));


    fprintf(1, '\nMin p-value = %3.6f, z = %3.2f\n', mp, norminv(1-mp));

    % --------------------------------------
    % * Threshold p-values
    % --------------------------------------
    if isinf(pthr)
        dofdr = 1;
        pthr = .05;
    end

    if dofdr
        pt = FDR(dat(eligible_vox), pthr);       % get threshold for FDR
        if isempty(pt) || isinf(pt), fprintf(1, 'No voxels pass FDR threshold.\n'); return, end
        fprintf(1, 'FDR: threshold is p < %3.4f\n', pt);
    else
        % uncorrected threshold
        pt = pthr;
        % overkill  %[dat, volInfo] = iimg_threshold(dat, 'thr', [0 pthr], 'k', kthr);
    end

    dat(dat>pt) = 0;

    fprintf(1, '\nSignificant voxels at p < %3.6f:\nBefore size threshold: %3.0f ', pt, sum(dat>0));

    % height threshold based on p-value
    [dat, nvox] = iimg_cluster_extent(dat, volInfo, kthr);

    %fprintf(1, '\nAfter size threshold: %3.0f, in %3.0f clusters of size %3.0f to %3.0f voxels\n', sum(dat>0), length(nvox{1}), min(nvox{1}), max(nvox{1}));
    fprintf(1, '  After size threshold: %3.0f\n', sum(dat>0));

    % --------------------------------------
    % * get t-values
    % --------------------------------------

    % get t-values, masked with significant p- voxels
    dat = iimg_mask(dat, timg, volInfo);

    fprintf(1, ' + %3.0f significant voxels\n ', sum(dat>0)); fprintf(1, '- %3.0f significant voxels\n', sum(dat<0));

    % --------------------------------------
    % * write output images
    % --------------------------------------

    if writeoutput
        % write output images
        if dofdr, fdrstr = '_fdr'; else fdrstr = []; end
        pstr = num2str(pthr); wh = find(pstr == '.'); wh = wh(1); pstr = pstr(wh+1:end);
        toutname = [timg(1:end-4) '_thr_p' fdrstr pstr timg(end-3:end)];
        tposname = [timg(1:end-4) '_pos_mask_p' fdrstr pstr timg(end-3:end)];
        tnegname = [timg(1:end-4) '_neg_mask_p' fdrstr pstr timg(end-3:end)];
    end

    if writeoutput
        voldata = iimg_reconstruct_3dvol(dat, volInfo, 'outname', toutname);
        iimg_reconstruct_3dvol(dat>0, volInfo, 'outname', tposname);
        iimg_reconstruct_3dvol(dat<0, volInfo, 'outname', tnegname);
    else
        %voldata = iimg_reconstruct_3dvol(dat, volInfo);
    end

    % --------------------------------------
    % * make clusters
    % --------------------------------------

    % make clusters
    posdat = dat; posdat(posdat<0) = 0;
    negdat = dat; negdat(negdat>0) = 0;

    % still use size threshold in case a cluster above is both + and -
    clpos = iimg_indx2clusters(posdat, volInfo, eps, kthr);         % save vals > eps
    clneg = iimg_indx2clusters(negdat, volInfo, [-Inf eps], kthr);  % save vals between -Inf and eps

    for i = 1:length(clpos)
        clpos(i).Z_descrip = 't-values';
    end
    
    for i = 1:length(clneg)
        clneg(i).Z_descrip = 't-values';
    end
    
    if ~isempty(clpos), szp = cat(1, clpos.numVox); else szp = 0; end
    if ~isempty(clneg), szn = cat(1, clneg.numVox);  else szn = 0; end

    fprintf(1, 'Suprathreshold clusters: \n + %3.0f (%3.0f to %3.0f voxels)   ', length(clpos), min(szp), max(szp));
    fprintf(1, '\n - %3.0f (%3.0f to %3.0f voxels)\n', length(clneg), min(szn), max(szn));

    % --------------------------------------
    % * display
    % --------------------------------------

    if display && (~isempty(clpos) || ~isempty(clneg))
        % cluster_orthviews([clpos clneg], 'bivalent', 'overlay', overlay, 'handle', orth_handle, addstr);
        cluster_orthviews([clpos clneg], 'overlay', overlay, 'handle', orth_handle, addstr);
    end

    % --------------------------------------
    % * table
    % --------------------------------------
    fprintf(1, '\nPositive results\n----------------------------\n')
    cluster_table(clpos, 1, 0);
    fprintf(1, '\nNegative results\n----------------------------\n')
    cluster_table(clneg, 1, 0);
end