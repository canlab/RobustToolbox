function cl = robust_nonparam_interactive_scatterplot(meth,X,cl,wh_interest,image_names)
% cl = robust_nonparam_interactive_scatterplot(meth,X,cl,wh_interest,image_names)
%
% meth: Method: 
%       'max' plots max correlation, uses cl.Ymaxdata field (see
%       robust_nonparam_get_sigregions.m)
%       'mean' plots mean correlation; you must enter valid image_names
%
% wh_interest: which column of X is of interest
% Works on existing orthviews figure
%
% tor wager, aug. 06
%
% Examples:
% partial correlations of first column of R.X, for cl_pos clusters
% robust_nonparam_interactive_scatterplot('max',R.X,cl_pos,1);
%
% partial correlations of 1st col. with input images, mean of cluster
% robust_nonparam_interactive_scatterplot('mean',R.X,[A.cl_pos
% A.cl_neg],1,R.image_names);
%
%


%wh_interest = 1;

% ---------------------------------------------------
% setup model
% ---------------------------------------------------
nuis = X;
nuis(:,wh_interest) = [];
for i = 1:length(cl)
    cl(i).xdat = X(:,wh_interest); 
    cl(i).nuisance = nuis; 
end


% ---------------------------------------------------
% setup data
% ---------------------------------------------------
switch meth
    case 'max'
        % max effect
        for i = 1:length(cl), cl(i).timeseries = cl(i).Ymaxeffect(:,1); end
    case 'mean'
        cl = extract_contrast_data({image_names},cl);
        for i = 1:length(cl), cl(i).timeseries = cl(i).CONTRAST.data(:,1); end

    otherwise
        error('Unknown method')
end


%R.image_names = filename_get_new_root_dir(R.image_names,spm_get(-1),2);
% cl = extract_contrast_data({R.image_names},cl);
% for i = 1:length(cl)
% cl(i).xdat = X(:,1);
% cl(i).nuisance = X(:,2);
% cl(i).timeseries = cl(i).CONTRAST.data(:,1);
% end

% ---------------------------------------------------
% interactive scatterplot
% ---------------------------------------------------
[xdat,ydat] = cluster_interactive_scatterplot(cl);

doauto = input('Save pngs of clusters automatically? ');

if doauto
    
    % This automatically executes the windowbtnup function
    cl = cluster_export_pngs(cl,1,[],1);

end

return






% ---------------------------------------------------
% subfunctions needed for 'auto' mode only
% ---------------------------------------------------

function cluster_interactive_callback(varargin)

% get data stored in figure
% should contain: clusterfield, xlab, ylab
% if doesn't exist, create as default: 'timeseries'
fh = findobj('Tag','Graphics'); figure(fh);
data = guidata(fh);

N = fieldnames(data);
for i=1:length(N)
    eval([N{i} ' = data.' N{i} ';']);
end

% check for data in specified field of cl
if ~isfield(cl,clusterfield), display_error('nodata',clusterfield); end
if ~isfield(cl,'xdat'), display_error('noxdata',clusterfield); end

spm_orthviews_showposition;

% activate scatterplot axis
axes(axish);

% find closest cluster and return index
% ---------------------------------------------------
pos = spm_orthviews('Pos')';


% check to see if we're in a cluster
wh = 0; 
centers = cat(1,cl.mm_center);

% find closest cluster, based on center
d = distance(pos,centers); wh = find(d == min(d)); wh = wh(1);

% only accept if cursor is w/i 2 mm of a voxel in the cluster
d = distance(pos,cl(wh).XYZmm'); d = min(d); if d > 15, wh = 0;,end



if wh
    %cluster_table(cl(wh));
    fprintf(1,'Cl. %3.0f, Voxels: %3.0f, Coords: %3.0f, %3.0f, %3.0f\n',wh,cl(wh).numVox,cl(wh).mm_center(1), cl(wh).mm_center(2),cl(wh).mm_center(3));

    %spm_orthviews('Reposition',cl(wh).mm_center);
    
    % get data from structure
    ydat = cl(wh).(clusterfield);
    xdat = cl(wh).xdat;
    
    % check for data
    if isempty(ydat),display_error('nodata',clusterfield); end
    if isempty(xdat), display_error('noxdata',clusterfield), end
        
    cla;
  
    % get design matrix, with 1st column as regressor of interest
    dorobust = 1;
    if isfield(cl, 'nuisance') && ~isempty(cl(wh).nuisance)
        X = [xdat cl(wh).nuisance];
    else
        X = xdat;
    end
    if no_intercept(X), X(:,end+1) = 1; end
    
    %[tmp,tmp,r] = partialcor(X,ydat,1,1,dorobust);
    
    plot_correlation(X,ydat,'robust','xlabel',data.xlab,'ylabel',data.ylab);
    
%     if isfield(cl, 'nuisance') && ~isempty(cl(wh).nuisance)
%         [r,str,sig,ry,rx,h,rr] = prplot(ydat,[xdat cl(wh).nuisance],1,1,{'ko'});
%     else
%         plot_correlation_samefig(xdat,ydat,[],'ko',0,1);
%     end
    %axis auto;
%     xlabel(data.xlab);
%     ylabel(data.ylab);

    drawnow
else
    cla;
    str = sprintf('No nearby cluster. Distance is %3.0f',d);
    fprintf(1,str); pause(.5); fprintf(repmat('\b',1,length(str)));
end

return



function display_error(meth,clusterfield)

switch meth
    case 'nodata'
        disp(['You have asked to get data from cl.' clusterfield ', but there is no data there.']);
        error('For a fix, try cl = tor_extract_rois(my_image_name_matrix_here,cl);');
        
    case 'noxdat'
        error('You must assign the x-axis data to plot to cl(*).xdat);');
end

return


function val = no_intercept(X)

val = all(any(diff(X)));

return