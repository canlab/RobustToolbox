function R = robust_nonparam_results(R,myalpha)
% R = robust_nonparam_results(R,[corrected alpha level])
%
% Updated: Jan 2008, by Tor Wager, to add primary uncorrected thresholds

k = size(R.X,2);
v = size(R.correct.t,2);

if nargin < 2, myalpha = .05; end

% banner
fprintf(1,'\n* ======================================================================== *')
fprintf(1,'\n*                                                                          *')
fprintf(1,'\n*                 Running function: robust_nonparam_results.m              *')
fprintf(1,'\n* Robust regression with nonparametric correction for multiple comparisons *')
fprintf(1,'\n*                                                                          *')
fprintf(1,'\n*                        Tor Wager, version July 2006                      *')
fprintf(1,'\n*                                                                          *')
fprintf(1,'\n* ======================================================================== *\n')

fprintf(1,'\nGeneral information:\n')
if isfield(R,'savefile'), disp(['Associated file is ' R.savefile '.mat']), end
fprintf(1,'Mask of in-analysis voxels is: \n%s\n', R.maskname)
fprintf(1,'First data image name is: \n%s\n', R.image_names(1,:))
fprintf(1,'\n')
fprintf(1,'Intercept is param. %3.0f, called %s\n', R.wh_intercept,R.names{R.wh_intercept})
fprintf(1,'Tails for max T: %s\n', R.mytail)
fprintf(1,'Corrected alpha level is: %3.2f\n',myalpha)
fprintf(1,'\n')

if ~strcmp(R.volInfo.fname,R.maskname), error('R.volInfo and R.maskname do not match! Wrong mask info?'); end


% remove extra empty iterations
R = setup_iterations(R,k,0,R.volInfo.c);

fprintf(1,'Valid iterations for t-tests: %3.0f\n', size(R.maxt,2))
fprintf(1,'\n')

% -------------------------------------------------------------------
% Map-wise stats
% -------------------------------------------------------------------
if isfield(R,'mapwise'), R = rmfield(R,'mapwise'); end

R.mapwise.alpha = myalpha;
if ~isempty(R.maxf), dofstat = 1; else dofstat = 0; end

% critical values, 2-tailed
R.mapwise.tcrit = prctile(R.maxt',100*(1-myalpha));
if dofstat, R.mapwise.fcrit = prctile(R.maxf',100*(1-myalpha)); end

% max correct perm
ctmax = max(R.correct.t');
ctmin = min(R.correct.t');
if dofstat, cfmax = max(R.correct.f); end

% suprathreshold vox
R.mapwise.tsig = sign(R.correct.t') .* (abs(R.correct.t)' >= repmat(R.mapwise.tcrit,v,1));
R.mapwise.tpos = sum(R.mapwise.tsig > 0);
R.mapwise.tneg = sum(R.mapwise.tsig < 0);

if dofstat
    R.mapwise.fsig = R.correct.f >= R.mapwise.fcrit;
    R.mapwise.fcount = sum(R.mapwise.fsig);

    f = [cfmax; R.mapwise.fcrit; R.mapwise.fcount]';
end

t = [ctmax; ctmin; R.mapwise.tcrit; R.mapwise.tpos; R.mapwise.tneg]';

fprintf(1,'\nMapwise statistics\n-------------------------------\n');
if dofstat
    print_matrix(f,{'Fmax' 'Fcrit' 'Nsig'},R.names);
    fprintf(1,'\n');
end
print_matrix(t,{'Tmax' 'Tmin' 'Tcrit' 'Npos' 'Nneg'},R.names);
fprintf('\n');

fprintf('\nTmax\tmaximum t-value within search mask');
fprintf('\nTmin\tminimum t-value within search mask');
fprintf('\nTcrit\tT-value needed for map-wise correction');
fprintf('\nNpos\tNumber of voxels with positive t-values exceeding threshold');
fprintf('\nTneg\tNumber of voxels with negative t-values exceeding threshold');
fprintf('\n');

% -------------------------------------------------------------------
% Cluster-wise stats
% -------------------------------------------------------------------
if isfield(R,'svc'),R = rmfield(R,'svc'); end

numc = R.volInfo.c;


% critical values (2-tailed if "twotail" option run in analysis)
R.svc.desc = 'tcrit and fcrit are regions x regressors';

for c = 1:numc

    n(c,1) = sum(R.volInfo.clindx == c);

    R.svc.tcrit(c,:) = prctile(R.maxt_by_cl{c}',100*(1-myalpha));
    if dofstat, R.svc.fcrit(c,1) = prctile(R.maxf_by_cl{c}',100*(1-myalpha)); end

    % voxels in this cluster
    this_cl = R.volInfo.clindx == c;

    % index vectors of voxels only in this cluster, v x k or v x 1
    this_t = R.correct.t' .* repmat(this_cl,1,k);
    if dofstat, this_f = R.correct.f' .* this_cl; end

    % max correct perm
    ctmax(c,:) = max(this_t);
    ctmin(c,:) = min(this_t);
    if dofstat,cfmax(c,1) = max(this_f); end

    % suprathreshold vox
    R.svc.tsig{c} = sign(this_t) .* (abs(this_t) >= repmat(R.svc.tcrit(c,:),v,1));
    if dofstat, R.svc.fsig{c} = this_f >= R.svc.fcrit(c,1); end

    % clusters are columns, rows are voxels
    R.svc.tpos(c,:) = sum(R.svc.tsig{c} > 0);
    R.svc.tneg(c,:) = sum(R.svc.tsig{c} < 0);
    if dofstat, R.svc.fcount(c,1) = sum(R.svc.fsig{c}); end

    % cluster names
    R.cl_names{c} = ['Cl. ' num2str(c)];

    % clusters  **** need to do for each regressor, done in
    % robust_nonparam_displayregions
    try
        R.svc.t_clusters{c} = iimg_indx2clusters(max(R.svc.tsig{c}'),R.volInfo,R.svc.tcrit(c,1));
    catch
% %         disp('transpose issue')
% %         R.svc.t_clusters{c} = iimg_indx2clusters(max(R.svc.tsig{c}')',R.volInfo,R.svc.tcrit(c,1));
    end
    
    cl = iimg_indx2clusters(this_t(:,1),R.volInfo);
    xyz(c,:) = cl.mm_center;

end


if dofstat, f = [xyz n cfmax R.svc.fcrit  R.svc.fcount]; end

robust_nonpar_set_level_svc(R);

fprintf('\nClusterwise statistics (SVCs)\n-------------------------------\n');
fprintf('\nThis section shows results for each column in your design matrix')
fprintf('\nin order, and lists contiguous ROIs within your search mask,')
fprintf('\nindicating how many voxels within each show significant results ')
fprintf('\nwith small-volume correction within the ROI.');
fprintf('\n')

if dofstat
    print_matrix(f,{'x' 'y' 'z' 'nVox' 'Fmax' 'Fcrit' 'Nsig'},R.cl_names);
    fprintf(1,'\n');
end

for kk = 1:k
    t = [xyz n ctmax(:,kk) ctmin(:,kk) R.svc.tcrit(:,kk)  R.svc.tpos(:,kk) R.svc.tneg(:,kk)];
    fprintf(1,'\nt-values: %s\n',R.names{kk});
    print_matrix(t,{'x' 'y' 'z' 'nVox' 'Tmax' 'Tmin' 'Tcrit' 'Npos' 'Nneg'},R.cl_names);
    fprintf(1,'\n');
end

fprintf(1,'\n* ======================================================================== *\n')
fprintf(1,  '*                        End Results table output                          *')
fprintf(1,'\n* ======================================================================== *\n')





plot_thresh_by_iterations(R,myalpha);

return




% -------------------------------------------------------------------
% -------------------------------------------------------------------


% Subfunctions


% -------------------------------------------------------------------
% -------------------------------------------------------------------





% -------------------------------------------------------------------
% SETUP: iterations, remove empty
% -------------------------------------------------------------------

function [R,permindx] = setup_iterations(R,k,niter,c)

if isfield(R,'maxt') && isfield(R,'maxt_by_cl')
    % find empty
    whomit = find(any(R.maxt == 0));
    R.maxt(:,whomit) = [];
    for i = 1:c
        whomit = find(any(R.maxt_by_cl{c} == 0));
        try R.maxt_by_cl{c}(:,whomit) = []; catch disp('Warning: maxt_by_cl is wrong length! Results are invalid for this cluster'), end
    end
else
    R.maxt = [];
    R.maxt_by_cl = cell(1,c);
end

if isfield(R,'maxf') && isfield(R,'maxf_by_cl')
    % find empty
    whomit = find(R.maxf == 0);
    R.maxf(whomit) = [];
    for i = 1:c
        whomit = find(R.maxf_by_cl{c} == 0);
        try R.maxf_by_cl{c}(whomit) = []; catch disp('Warning: maxf_by_cl is wrong length. Results are invalid for this cluster'),end
    end
else
    R.maxf = [];
    R.maxf_by_cl = cell(1,c);
end

startat = size(R.maxt,2) + 1;

newt = zeros(k,niter);
newf = zeros(1,niter);
R.maxt = [R.maxt newt];     % max abs t (2-tailed);
R.maxf = [R.maxf newf];

% by-cluster
for i = 1:c
    R.maxt_by_cl{c} = [R.maxt_by_cl{c} newt];
    R.maxf_by_cl{c} = [R.maxf_by_cl{c} newf];
end

endat = size(R.maxt,2);
permindx = startat:endat;   % which indices in output arrays for each permutation

% the code below is different from that in robust_reg_nonparam:
texist = size(R.maxt,2) + 1;
fexist = size(R.maxf,2) + 1;
fprintf(1,'\nExisting iterations: f: %3.0f  t: %3.0f.\n',fexist,texist);

if min(fexist,texist) < 4000
    fprintf(1,'\n**************************************************\n')
    fprintf(1,'\nWarning:\n')
    if fexist < 4000
        fprintf(1,'Results for f-maps have < 4000 iterations, and are not likely to be stable\n')
    end
    if texist < 4000
        fprintf(1,'Results for t-maps have < 4000 iterations, and are not likely to be stable\n')
    end

    fprintf(1,'\n**************************************************\n')
end

return


% -------------------------------------------------------------------
% DATA: check for R.Y and load if necessary.  Check image names
% -------------------------------------------------------------------

function [R,Y] = load_data_subfunction(R)
sprintf('Loading data. ');

if isfield(R,'Y') && ~isempty(R.Y)
    disp('R.Y data field found.  Using it.');
    Y = R.Y;

    if ~exist(deblank(R.maskname),'file')
        fprintf(1,'\n* Error: R.maskname is not a valid file.\nImage: %s\nPlease select.\n',R.maskname)
        R.maskname = spm_get(1,'*','Select new mask image');
    end

else
    % check image names and get data

    if ~exist(deblank(R.image_names(1,:)),'file')
        fprintf(1,'\n* Error: R.image_names are not valid files.  Assuming names are correct but path is wrong!\n')
        R.image_names = filename_get_new_root_dir(R.image_names, ...
            spm_get(-1,'*','Select root dir for files'),2);
    end

    if ~exist(deblank(R.maskname),'file')
        fprintf(1,'\n* Error: R.maskname is not a valid file.\nImage: %s\nPlease select.\n',R.maskname)
        R.maskname = spm_get(1,'*','Select new mask image');
    end

    % load data and save in R
    Y = iimg_get_data(R.maskname,R.image_names);
    R.Y = Y;
end
if ~isfield(R,'volInfo') || ~strcmp(R.volInfo.fname,R.maskname)
    R.volInfo = iimg_read_img(R.maskname,1);
else
    fprintf(1,'\n* Using existing volInfo structure. Remove before running if mask info has changed.');
end

return




% -------------------------------------------------------------------
% PLOTS: not used yet
% -------------------------------------------------------------------

function plot_thresh_by_iterations(R,myalpha)

n = size(R.maxt,2);
k = size(R.X,2);
stepsize = round(n ./ 100);
myperc = 100* (1 - myalpha);

x = 10:stepsize:n;
tor_fig;
%colors = {'r' 'b' 'g' 'k' 'm' 'c' 'y'};

thr = zeros(length(x),k);
for i = 1:length(x)
    thr(i,:) = prctile(R.maxt(:,1:x(i))',myperc);
end
plot(x,thr,'LineWidth',2);


legend(R.names);
title('Mapwise t-threshold by number of iterations');
xlabel('Iterations'), ylabel('Critical t');

return


