% function robust_results3(EXPT, u, k, ['display thresholds', thresholds], ['overlay', ovl], ['result names', resultnames], ...
%   ['write files', 0|1], ['contrast name', contrastname], ['cov names', covariatenames], ['pause for display', 0|1])
%
% Displays the results of a random effects analysis contained in a robust* directory. By default, saves .png
% files in axial and medial views.
%
% u: uncorrected threshold (e.g., .001) or 'FDR' % only for rob*img files that are produced, NOT for display
% k: extent threshold (minimum size of cluster)
%
% E.g.:
% % After loading an EXPT object...
% robdirs = filenames('robust00[0-9][0-9]');
% for i=1:length(robdirs)
%     cd(robdirs{i});
%     robust_results3(EXPT, .005, 10, 'result names', EXPT.SNPM.connames{i});
%     close all;
%     cd('..');
% end
%
% To not display FDR, run:
% for i=1:length(robdirs)
%     cd(robdirs{i});
%     robust_results3(EXPT, .005, 10, 'result names', EXPT.SNPM.connames{i}, 'display thresholds', [.001 .005]);
%     close all;
%     cd('..');
% end
%
% To pause after each one, run:
% for i=1:length(robdirs)
%     cd(robdirs{i});
%     robust_results3(EXPT, .005, 10, 'result names', EXPT.SNPM.connames{i}, 'pause for display', 1);
%     close all;
%     cd('..');
% end
%
% To not write out result files (cls and imgs), run:
% for i=1:length(robdirs)
%     cd(robdirs{i});
%     robust_results3(EXPT, .005, 10, 'result names', EXPT.SNPM.connames{i}, 'write files', 0);
%     close all;
%     cd('..');
% end

function robust_results3(EXPT, u, k, varargin)
    [tmp currentdir] = fileparts(pwd()); %#ok - this is a signal to mlint to ignore the fact that tmp is never used

    result_names = {'Intercept' 'Cov 1' 'Cov 2' 'Cov 3' 'Cov 4' 'Cov 5' 'Cov 6' 'Cov 7' 'Cov 8' 'Cov 9' 'Cov 10' 'Cov 11' 'Cov 12'};
    pause_for_display = 0;
    write_files = 1;
    thresholds = [Inf .001 .005];
    ovl = [];


    for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case {'contrast name', 'contrastname'}
                    result_names{1} = varargin{i+1};
                case {'cov names', 'covnames'}
                    result_names = {'Intercept' varargin{i+1}{:}};
                case {'result names', 'resultnames'}
                    result_names = cellstr(varargin{i+1});
                case 'display thresholds'
                    thresholds = varargin{i+1};
                case {'overlay', 'ovl'}
                    ovl = varargin{i+1};
                case {'pausefordisplay' 'pause for display'}
                    pause_for_display = varargin{i+1};
                case {'writefiles' 'write files'}
                    write_files = varargin{i+1};
            end
        end
    end

    warning off

    % convert u to t-threshold
    df = length(EXPT.subjects) - length(dir('rob_p_0*.img'));
    fprintf(1, '\n\n\nrobust_results3 thinks there are %d degrees of freedom.\n', df);

    uorig = u;

    % if uncorrected threshold
    if ~ischar(u)
        t = tinv(1-u, df);
        fprintf(1, 'Height thresh: t = %3.2f (%3.0f Ss, %3.0f df @ p < %3.4f, extent = %3.0f\n', t(1), length(EXPT.subjects), df, uorig(1), k(1));
    else
        uorig = 0.05;
    end

    if isempty(ovl) && isfield(EXPT, 'overlay')
        disp(['Using overlay image: ' EXPT.overlay]);
        ovl = EXPT.overlay;
    end

    if(write_files)
        save_cl_results();
    end
    display_activations();
    

    warning on

    
    
    %------------------%
    % Nested functions %
    %------------------%
    
    function save_cl_results()
        for i=1:length(result_names)
            tmap_file = sprintf('rob_tmap_%04d.img', i);
            if exist(tmap_file, 'file')
                switch i
                    case 1 % intercept results
                        headers = {'Overall pos activation' 'Overall neg activation'};
                        mat_files = {[currentdir '_intercept_pos'] [currentdir '_intercept_neg']};
                    otherwise
                        headers = {sprintf('%s pos effect', result_names{i}) sprintf('%s neg effect', result_names{i})};
                        mat_files = {sprintf('%s_%04d_pos', currentdir, i) sprintf('%s_%04d_neg', currentdir, i)};
                end

                fprintf(1, '\n---------------------------------------------\n%s\n---------------------------------------------\n', headers{1});
                %FDR threshold
                if (ischar(u) && strcmp(u, 'FDR')) || isinf(u)
                    t = spm_uc_FDR(.05, [1 df], 'T', 1, spm_vol(tmap_file), 0);
                    fprintf(1, 'Height thresh FDR-corr: t = %3.2f (%3.0f Ss, %3.0f df @ p < %3.4f, extent = %3.0f\n', ...
                        t(1), length(EXPT.subjects), df, uorig, k(1));
                end

                pos_thresh_imgs = threshold_imgs(tmap_file, t(1), k(1), 'pos');
                cl = mask2clusters(pos_thresh_imgs); %#ok
                save(mat_files{1}, 'cl');

                fprintf(1, '\n---------------------------------------------\n%s\n---------------------------------------------\n', headers{2})
                neg_thresh_imgs = threshold_imgs(tmap_file, t(1), k(1), 'neg');
                cl = mask2clusters(neg_thresh_imgs); %#ok
                save(mat_files{2}, 'cl');
            end
        end
    end

    function display_activations()
        result_names = cellstr(result_names);
        thr_string = strtrim(num2str(thresholds, '%0.3g '));
        thr_string = strrep(thr_string, 'Inf', 'FDR');
        for i=1:length(result_names)
            tmap_file = sprintf('rob_tmap_%04d.img', i);
            if exist(tmap_file, 'file')
                fprintf(1, '\n---------------------------------------------\nDisplaying %s\n---------------------------------------------\n', result_names{i});
                fig_title = sprintf('%s - %s - %s', currentdir, result_names{i}, thr_string);
                if(~isempty(thresholds))
                    cl = multi_threshold2(tmap_file, 'T', df, 'overlay', ovl, 'title', fig_title, ...
                        'save images', 1, 'thresholds', thresholds); %#ok
                else
                    cl = multi_threshold2(tmap_file, 'T', df, 'overlay', ovl, 'title', fig_title, ...
                        'save images', 1); %#ok
                end
                if(pause_for_display)
                    input('Press return/enter to continue');
                end
            end
        end
    end
end


