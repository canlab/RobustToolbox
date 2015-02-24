function varargout = robust_toolbox_gui(Action, varargin)
    % varargout = robust_toolbox_gui(Action, varargin)
    %
    % Type >>robust_toolbox_gui to get started
    % See ROBUST_REGRESSION_HELP.pdf for documentation.
    % See also robfit.m
    %
    % by Tor Wager
    %
    % Thanks to Tom Nichols for the excellent GUI shell!

    %-----------------------------functions-called------------------------
    %
    %-----------------------------functions-called------------------------


    %-Format arguments
    %-----------------------------------------------------------------------
    if nargin == 0, Action='Init'; end



    % Load EXPT and cl in 'nonverbose' mode
    % - from workspace (preferred) or from EXPT.mat if it exists on path
    % - unless Action is 'Init', in which case it will be loaded in verbose
    % mode
    
    if ~strcmpi(Action, 'init')
    
        [EXPT, cl] = load_expt_variable(0);

    end


    switch lower(Action)

        case lower('Init')
            %=======================================================================

            clc
            %BrainVowager_defaults;

            robust_toolbox_gui('AsciiWelcome')
            

            [EXPT, cl] = load_expt_variable(1);

            robust_toolbox_gui('CreateMenuWin')
            varargout{1} = EXPT;

        case lower('AsciiWelcome')
            %=======================================================================
            disp('==================================================')
            disp( 'Welcome to robust_toolbox_gui: Robust Analyses.')
            disp('==================================================')
            fprintf('\n')


        case lower('Ver')
            %=======================================================================
            varargout = {'SCNlab Robust Menu'};



        case lower('CreateMenuWin')
            %=======================================================================
            %close(findobj(get(0, 'Children'), 'Tag', 'robust_toolbox_gui Menu'))


            %-Initialize robust_toolbox_gui menu window
            %-----------------------------------------------------------------------
            [F, winwid, winh] = robust_toolbox_gui('initFigure');


            % default button sizes and positions, etc.

            topbutton = winh-100;        % y location of top button
            butspace = 35;               % y spacing of buttons

            fullbutxy = [160 30];       % full-length button width and height
            halfbutxy = [80 30];        % (left-hand) half-width button w and h
            rightbutxy = halfbutxy;
            rightbutx = 110+halfbutxy(1)+5;  % right-hand button start x


            %-Frames and text
            %-----------------------------------------------------------------------
            axes('Position', [0 0 80/winwid winh/winh], 'Visible', 'Off')
            text(0.5, 0.475, 'Robust Toolbox', ...
                'FontName', 'Times', 'FontSize', 36, ...
                'Rotation', 90, ...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
                'Color', [1 1 1]*.6);

            text(0.2, 0.96, 'SCN Lab', ...
                'FontName', 'Times', 'FontSize', 16, 'FontAngle', 'Italic', ...
                'FontWeight', 'Bold', ...
                'Color', [1 1 1]*.6);

            uicontrol(F, 'Style', 'Frame', 'Position', [095 005 winwid-100 winh - 30], ...
                'BackgroundColor', robust_toolbox_gui('Color'));  % colored frame
            uicontrol(F, 'Style', 'Frame', 'Position', [105 015 winwid-120 winh - 50]);  % inner gray frame

            %-Buttons to launch robust_toolbox_gui functions
            %-----------------------------------------------------------------------

            % -------------------------------------------
            % Section - Random Effects Analysis
            uicontrol(F, 'Style', 'Text', ...
                'String', 'Random Effects', 'FontSize', 14, ...
                'HorizontalAlignment', 'Center', ...
                'Position', [115 topbutton+30 fullbutxy], ...
                'ForegroundColor', 'y', 'FontWeight', 'b');
            % -------------------------------------------

            % Robust setup
            str = 'EXPT = CreateExpt(''robfit'');';                            % callback function
            str = robust_toolbox_gui('ExpandString', str);               % display then execute
            uicontrol(F, 'String', 'Setup', ...
                'Position', [110 topbutton-(1-1)*butspace fullbutxy], ...
                'CallBack', str, ...
                'Interruptible', 'on', ...
                'ForegroundColor', 'k', 'FontWeight', 'b');

            % Robust RFX Analysis
            str = 'EXPT = robfit(EXPT, 1:length(EXPT.SNPM.P), 0, EXPT.mask);';            % callback function
            str = robust_toolbox_gui('ExpandString', str);               % display then execute
            uicontrol(F, 'String', 'Robust RFX', ...
                'Position', [110 topbutton-(2-1)*butspace fullbutxy], ...
                'CallBack', str, ...
                'Interruptible', 'on', ...
                'ForegroundColor', 'k', 'FontWeight', 'b');

            % Robust Results
            str1 = 'hthr = input(''Enter p-value threshold for saving clusters: '');';
            str2 = 'sthr = input(''Enter extent threshold for saving clusters: '');';
            str3 = 'dthr = input(''Enter display threshold(s): '');';
            str = 'robust_results3(EXPT, hthr, sthr, ''display thresholds'', dthr);' ;
            str = [str1 str2 str3 str];
            
            str = robust_toolbox_gui('ExpandString', str);
            uicontrol(F, 'String', 'Robust Results', ...
                'Position', [110 topbutton-(3-1)*butspace fullbutxy], ...
                'CallBack', str, ...
                'Interruptible', 'on', ...
                'ForegroundColor', 'k', 'FontWeight', 'b');

            % -------------------------------------------
            % Next section - Robust Seed Analysis
            uicontrol(F, 'Style', 'Text', ...
                'String', 'Robust Seed Analysis', 'FontSize', 14, ...
                'HorizontalAlignment', 'Center', ...
                'Position', [115 topbutton-(4.2-1)*butspace fullbutxy], ...
                'ForegroundColor', 'y', 'FontWeight', 'b');
            % -------------------------------------------


            % Select Seeds - opens up scn_roi_gui to: 1) set ROIs, 2) extract data,
            % 3) save seeds
            str1 = 'disp(''You need to get clusters of interest and extract data from them to use as seeds.''); ';
            str = [str1 'scn_roi_gui;'];                                         % callback function
            str = robust_toolbox_gui('ExpandString', str);                 % display then execute
            uicontrol(F, 'String', 'Cluster Menu (Setup)', ...
                'Position', [110 topbutton-(5-1)*butspace fullbutxy], ...
                'CallBack', str, ...
                'Interruptible', 'on', ...
                'ForegroundColor', 'k', 'FontWeight', 'b');

            % Robust Seed
            %str1 = 'disp(EXPT.SNPM.connames), whc = input(''Which contrasts? (e.g., [1:5]) : ''); ';
            str1 = 'doglobal = input(''Use global contrast covariate? (1/0) ''); ';
            str = [str1 'EXPT = robseed(EXPT, cl, 1:length(EXPT.SNPM.P), 0, EXPT.mask, doglobal);'];
            str = robust_toolbox_gui('ExpandString', str);                 % display then execute
            uicontrol(F, 'String', 'Robust Seed Analysis', ...
                'Position', [110 topbutton-(6-1)*butspace fullbutxy], ...
                'CallBack', str, ...
                'Interruptible', 'on', ...
                'ForegroundColor', 'k', 'FontWeight', 'b');

            % Robust Seed Results 
            
            str1 = 'hthr = input(''Enter p-value threshold for saving clusters: '');';
            str2 = 'sthr = input(''Enter extent threshold for saving clusters: '');';
            str3 = 'dthr = input(''Enter display threshold(s): '');';
            str = 'robust_results3(EXPT, hthr, sthr, ''display thresholds'', dthr);' ;
            str = [str1 str2 str3 str];
            
            %str = 'robust_results2(EXPT, .005, 10);' ;
            str = robust_toolbox_gui('ExpandString', str);
            uicontrol(F, 'String', 'Robust Seed Results', ...
                'Position', [110 topbutton-(7-1)*butspace fullbutxy], ...
                'CallBack', str, ...
                'Interruptible', 'on', ...
                'ForegroundColor', 'k', 'FontWeight', 'b');


            % -------------------------------------------
            % Next section - display
            uicontrol(F, 'Style', 'Text', ...
                'String', 'Display results map', 'FontSize', 14, ...
                'HorizontalAlignment', 'Center', ...
                'Position', [115 topbutton-(8-1)*butspace fullbutxy], ...
                'ForegroundColor', 'y', 'FontWeight', 'b');
            % -------------------------------------------

            buttontext = {'Display Menu' ...
                'Threshold/Display' ...
                'Activation + correlation' ...
                'Cluster tools menu' ...
                'Intercept' ...
                'Cov 1' ...
                'Cov 2' ...
                'Cov 3' ...
                };
            str = 'robust_toolbox_gui(''results'');';                          % callback function
            pop1 = uicontrol(F, 'Style', 'popupmenu', 'Tag', 'ResultsPop', ...
                'String', buttontext, ...
                'Position', [110 topbutton-(9-1)*butspace+15 fullbutxy], ...
                'CallBack', str, ...
                'Interruptible', 'on', ...
                'ForegroundColor', 'k', 'FontWeight', 'b');




            set(F, 'Pointer', 'Arrow', 'Visible', 'on')




        case lower('Color')
            %=======================================================================
            % robust_toolbox_gui('Color')
            %-----------------------------------------------------------------------
            % %-Developmental livery
            % varargout = {[0.7, 1.0, 0.7], 'Lime Green'};
            %-Distribution livery
            varargout = {[.8 0.7 1.0], 'Purple'};


        case lower('ExpandString')
            %=======================================================================
            % robust_toolbox_gui('ExpandString')
            % Expand an action button callback string (a command string to be
            % evaluated)
            % so that it first displays the command, and then executes it
            %-----------------------------------------------------------------------
            str = varargin{1}; str2 = [];
            for i = 1:length(str)
                if str(i) == ''''
                    str2(end+1) = '''';
                    str2(end+1) = '''';
                else
                    str2(end+1) = str(i);
                end
            end

            str = ['disp(''' char(str2) '''), ' str ];     % display then execute

            varargout = {str};


        case lower('initFigure')
            %=======================================================================
            % [F, winwid, winh] = robust_toolbox_gui('initFigure')
            %-----------------------------------------------------------------------
            % Get the position of the main BrainVowager menu, or if
            % not available, default screen pos.

            % default sizes, etc.
            S = get(0, 'ScreenSize');

            winwid = 300;               % window width
            winh = 400;                 % window height
            pos = [S(3)/2+150, S(4)/2-140, winwid, winh];  % default

            % look for figures to place this one next to, in this order (last is preferred):
            h = [];
            h = findobj('Tag', 'BrainVowager_gui Menu');
            h = findobj('Tag', 'scn_fir_results_gui Menu');

            if ~isempty(h),
                pos = get(h, 'Position');
                winwid = pos(3);
                winh = pos(4);
                pos(1) = pos(1) + winwid;   % put next to main figure
            end

            %-Open robust_toolbox_gui menu window
            %----------------------------------------------------------------------
            F = create_figure('robust_toolbox_gui Menu');

            set(F, 'Color', [1 1 1]*.8, ...
                'Name', robust_toolbox_gui('Ver'), ...
                'NumberTitle', 'off', ...
                'Position', pos, ...
                'Resize', 'off', ...
                'Tag', 'robust_toolbox_gui Menu', ...
                'Pointer', 'Watch', ...
                'MenuBar', 'none', ...
                'Visible', 'off');

            set(get(F, 'CurrentAxes'), 'Visible', 'off');

            varargout{1} = F; varargout{2} = winwid; varargout{3} = winh;


        case lower('results')
            %=======================================================================
            % robust_toolbox_gui('results')
            %----------------------------------------------------------------------

            if isempty(EXPT), warning('EXPT is empty.  not created, or not declared global in base workspace?'); end

            % get overlay image, if we have the ROI gui open
            han = findobj('Tag', 'scn_roi_gui Menu');
            if ~isempty(han)
                scn_roi_gui('overlay');
                P = guidata(han);
            else
                P = []; % default overlay
            end

            % callbacks
            % buttontext = {'Display Menu' 'Threshold/Display' 'Thresh/Disp with Scatterplots' 'Intercept' 'cov1' 'Cluster tools menu'};
            pop1 = findobj('Tag', 'ResultsPop');
            indx = get(pop1, 'Value');

            % callback functions

            callbk = {@noop_callback ...
                @robust_results_threshold_callback ...
                @active_plus_corr_scatterplot_plugin ...
                @scn_roi_gui ...
                @intercept_callback ...
                @cov1_callback ...
                @cov2_callback ...
                @cov3_callback ...
                };
            feval(callbk{indx});

            
        otherwise
            error('Unknown action string')
    end
    
    % return EXPT and cl to base workspace
    assignin('base', 'EXPT', EXPT);
    assignin('base', 'cl', cl);
    
    
end

%--------------------
% Callback functions
%--------------------
function noop_callback()
end

function robust_results_threshold_callback()
    robust_results_threshold_wrapper();
end

function intercept_callback()
    robust_results_threshold_wrapper('intercept');
end

function cov1_callback()
    robust_results_threshold_wrapper('cov1');
end

function cov2_callback()
    robust_results_threshold_wrapper('cov2');
end

function cov3_callback()
    robust_results_threshold_wrapper('cov3');
end

function robust_results_threshold_wrapper(cov_or_int_string)
    
    pthr = spm_input('Enter p-value threshold (uncorrected). Type "Inf" (without quotes) for FDR: ', 1);
    kthr = spm_input('Enter extent threshold (min cluster size) in voxels: ');
    display_string = spm_input('Display results? ', [], 'y/n', {''; 'nodisplay'}, 1);
    display_string = display_string{1};

    if(~exist('cov_or_int_string', 'var') || isempty(cov_or_int_string))
        pimg = scan_get_files(1, 'rob_p*img', 'Select a p-image file to get clusters from', pwd);
        
        [clpos, clneg] = robust_results_threshold(pthr, kthr, 'p', pimg, display_string);
    else
        [clpos, clneg] = robust_results_threshold(pthr, kthr, display_string, cov_or_int_string);
    end
    
    disp(' ');
    disp('Created variables clpos and clneg in the base workspace.')
    disp('These contain info (clusters) for positive-response ')
    disp('and negative-response regions, respectively.')
    disp('You can save these to disk, and/or display using:')
    disp('cluster_orthviews, cluster_surf, montage_clusters, or other functions');
    
    % load data
    if exist(fullfile(pwd, 'SETUP.mat'))
        load SETUP

        for i = 1:size(SETUP.files, 1)
            isfile(i) = exist(deblank(SETUP.files(i,:)), 'file');
        end
    
        if sum(isfile > 0) == size(SETUP.files, 1)
            % ok
            clpos = tor_extract_rois(SETUP.files, clpos);
            clneg = tor_extract_rois(SETUP.files, clneg);
            
            disp('clpos.timeseries and clneg.timeseries contain data ')
            disp('averaged over voxels within each suprathreshold region.')
            
        else
            disp('Warning: cannot find image files.  Update in SETUP.files');
        end
        
    end
           
    assignin('base', 'clpos', clpos);
    assignin('base', 'clneg', clneg);
    

end




function [EXPT, cl] = load_expt_variable(verbose)
    % load EXPT

    try
        EXPT = evalin('base', 'EXPT');
        if verbose && ~isempty(EXPT), disp('Using EXPT already in base memory.'); end
    catch
        EXPT = [];
    end

    try
        cl = evalin('base', 'cl');
        if verbose && ~isempty(cl), disp('Retrieved cl from base workspace.'); end
    catch
        cl = [];
    end

  
    if isempty(EXPT) && exist(fullfile(pwd, 'EXPT.mat'), 'file')
        if verbose, disp('loading EXPT.mat from file'); end
        load EXPT;

    elseif isempty(EXPT) && ~exist(fullfile(pwd, 'EXPT.mat'), 'file')
        if verbose
            disp('You need to go to the directory containing individual subject results dirs ');
            disp('this directory should have information about the experiment stored in ')
            disp('a variable called EXPT and saved in the file EXPT.mat')
            disp('No EXPT file in current directory.  Create manually or with Setup button.');

            fprintf(1, '\n')
        end
        EXPT = [];

    end

    if verbose
        fprintf(1, '\n')

        disp('==================================================')
        disp('EXPT is a structure that contains information about')
        disp('the images and design.');

        disp('Current values of EXPT:');
        if isempty(EXPT), disp('Empty.'); else disp(EXPT); end

        %assignin('base', 'EXPT', EXPT);
        disp('==================================================')

        fprintf(1, '\n')
        disp('==================================================')
        disp('cl is a structure that contains information about ')
        disp('results ''blobs'' and is used for ROI stats and') 
        disp('visualization.');

        disp('Current values of cl:');
        if isempty(cl), disp('Empty.'); else disp(cl); end
        %assignin('base', 'cl', cl);
        disp('==================================================')
        fprintf(1, '\n')

        disp('==================================================')
        disp('More notes on EXPT:');
        disp('After creating, you can manually edit EXPT in the');
        disp('workspace and save it in the file EXPT.mat');
        disp('Enter covariates in 2nd level robust analysis as');
        disp('columns in EXPT.cov');
        disp('If these are contrast coded, they will be treated');
        disp('as categorical predictors which can have unequal');
        disp('N in each group.  If they are anything else,');
        disp('They will be centered during analysis to make them');
        disp('independent of the intercept.');
        disp('==================================================')
        fprintf(1, '\n')
    end

end