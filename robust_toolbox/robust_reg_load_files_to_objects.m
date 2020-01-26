function [trob, names, mask_obj, nsubjects, weights, SETUP] = robust_reg_load_files_to_objects(robust_reg_directory)

% Load files from a CANLab robust regression directory into objects
% Objects contain information needed for results display and tables.
%
% trob      statistic_image object with one image per contrast (the model intercept is the first image)
% names     cell array of names of each image (intercept, regressor 1, etc.)
% mask_obj  image defining the set of voxels analyzed
% nsubjects image summarizing number of subjects with valid data in each voxel
% weights   fmri_data object with weight maps for each subject (i.e., input image)
% SETUP     struct saved in directory, including design matrix SETUP.X

if ~exist(fullfile(robust_reg_directory, 'SETUP.mat'), 'file')
    error('robust_reg_directory must be a CANLab robust regression directory with SETUP.mat');
end


dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

printhdr('Robust regression');
fprintf('Loading files from %s\n', robust_reg_directory);
disp(' ');

load SETUP
disp(char({SETUP.dir SETUP.name}));


mask_obj = fmri_data(fullfile(robust_reg_directory, 'mask.nii'), 'noverbose');
nsubjects = fmri_data(fullfile(robust_reg_directory, 'nsubjects.nii'), 'noverbose');
weights = fmri_data(fullfile(robust_reg_directory, 'weights.nii'), 'noverbose');

[n, k] = size(SETUP.X);

names = {'Group average (Intercept)'};
for i = 2:k, names = [names sprintf('Regressor %d', i - 1)]; end

% Load all t images as statistic_image objects
d = dir('rob_tmap*');

trob = statistic_image('image_names', char(d.name), 'type', 't', 'dfe', n - k);

% Replace P-values with those from robust reg, as dfe varies
d = dir('rob_p_*');
prob = statistic_image('image_names', char(d.name));

trob = replace_empty(trob);
prob = replace_empty(prob);
trob.p = prob.dat;
trob = remove_empty(trob);

% Print some output
% --------------------------------------------------------------


if isfield(SETUP, 'name') && ~isempty(SETUP.name)
    
    printstr(SETUP.name)
    
end

if isfield(SETUP, 'dir') && ~isempty(SETUP.dir)
    
    printstr(SETUP.dir)
    
end

printhdr('Design Matrix')

ynstr = {'No' 'Yes'};
ok2string = @(isok) ynstr{double(isok) + 1} ;

isok = all(SETUP.X(:, 1) == mean(SETUP.X(:, 1), 1));
fprintf('First regressor (image) is intercept: %s\n', ok2string(isok));

isok = all(abs(mean(SETUP.X(:, 2:end))) < 1000 * eps)  ;
fprintf('Regressors mean-centered (intercept reflects group mean): %s\n', ok2string(isok));

fprintf('Number of regressors (including intercept): %d\n', k);

printstr(' ');


end
