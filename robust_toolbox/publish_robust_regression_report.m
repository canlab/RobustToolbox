% Runs batch analyses and publishes HTML report with figures and stats to
% results/published_output in local study-specific analysis directory.

% Run this from the main CANLab robust regression results directory (basedir)

close all

resultsdir = pwd;

fprintf('Creating HTML report for results in:\n%s\n', resultsdir);

% ------------------------------------------------------------------------
% Check that we have a valid mediation results dir
% ------------------------------------------------------------------------

is_robust_dir = exist(fullfile(resultsdir, 'SETUP.mat'));
is_run = exist(fullfile(resultsdir, 'rob_tmap_0001.nii')) | exist(fullfile(resultsdir, 'rob_tmap_0001.img'));

if ~is_robust_dir
    fprintf('%s\nis not a CANlab robust results directory because SETUP.mat is missing.\nSkipping report.\n', resultsdir);
    return
end

if ~is_run
    fprintf('Robust regression (robfit.m) does not appear to have run correctly because rob_tmap_0001.nii is missing.\nSkipping report.\n');
    return
end

% ------------------------------------------------------------------------
% Set HTML report filename and options
% ------------------------------------------------------------------------

pubdir = fullfile(resultsdir, 'published_output');
if ~exist(pubdir, 'dir'), mkdir(pubdir), end

pubfilename = ['robust_regression_report_' scn_get_datetime];

p = struct('useNewFigure', false, 'maxHeight', 800, 'maxWidth', 1600, ...
    'format', 'html', 'outputDir', fullfile(pubdir, pubfilename), 'showCode', false);

% ------------------------------------------------------------------------
% Run and report status
% ------------------------------------------------------------------------

publish('robust_results_batch', p)

myhtmlfile = fullfile(pubdir, pubfilename, 'robust_results_batch.html');

if exist(myhtmlfile, 'file')
    
    fprintf('Saved HTML report:\n%s\n', myhtmlfile);
    
    web(myhtmlfile);
    
else
    
    fprintf('Failed to create HTML report:\n%s\n', myhtmlfile);
    
end

